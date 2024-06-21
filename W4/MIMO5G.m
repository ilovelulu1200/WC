clear
clc
close all

Tx =2;         %傳送端個數
Rx =2;         %接收端個數


%% setting
data_num = 200;               % data量(200frames/SNR)
SNR_in_dB = 0:5:40;         % 自己設訂雜訊大小

BER_SNR_ZF=zeros(1,length(SNR_in_dB));

QAM = 16;
Eavg = 10;
NF = 1/sqrt(Eavg);
q_bit = log2(QAM);        % 一個symbol可以傳幾個bit
N = Tx;

%% start SNR 0:5:40
for a=1:length(SNR_in_dB)
    SNR = 10^(SNR_in_dB(a)/10);
    No = 10^(-SNR_in_dB(a)/10);
    Rx1_No = No; %之後做雜訊正規化需要使用
    Rx2_No = No;
    Es = 1;
    BER_ZF = 0;                        

    %% start 200 frames / SNR
    for seperate = 1:data_num
        data = randi([0 QAM-1],12,14,Tx);          % 隨機產生0~15 for 16QAM 矩陣大小12*14*2
        data1=repmat(data,137,40,1); %擴增矩陣為1644*560*2
        bin_data = de2bi(data1,q_bit,'left-msb');      % 0~15 to binary
        X = qammod(data1,QAM,'gray','UnitAveragePower',true);             % 0~15 to complex (Modulation); remember to normalize
        X_G=[zeros(202,560,Tx);X(1:822,:,:);zeros(1,560,Tx);X(823:1644,:,:);zeros(201,560,Tx)]; %add guard band
        x=ifft(ifftshift(X_G,1))*sqrt(2048); % phase to time

        %% add CP
        x_CP = zeros(1,1228800,Tx);
        index=1;
        for symbol=1:14*4*10
            if	mod(symbol,28)-1
                x_CP(1, index:index+2048+144-1,:)=[ x(2048-144+1:2048,symbol,:) ; x(:,symbol,:)];
                index = index+2048+144;
            else
                x_CP(1, index:index+2048+208-1,:)=[ x(2048-208+1:2048,symbol,:) ; x(:,symbol,:)];
                index = index+2048+208;
            end
        end

        %% channel
        Power_dB=[-2,-8,-10,-12,-15,-18].';
        Power_dB_MIMO = repmat(Power_dB,1,Rx,Tx);
        H_Channel=sqrt(10.^(Power_dB_MIMO./10)).*(sqrt(1/(2*Tx)).*(randn(6,Rx,Tx)+1i*randn(6,Rx,Tx)));
        h=[H_Channel;zeros(2042,Rx,Tx)];
        H=fftshift(fft(h,[],1),1);
        H_data=[H(203:1024,:,:);H(1026:1847,:,:)];
        H_frame = permute(repmat( H_data(:,:,:),1,1,1,560),[1 4 2 3]);

        %% conv with channel & add noise
        n = randn(Tx,1228800)*sqrt(No/2)+randn(Tx,1228800)*1i*sqrt(No/2) ;             % randn產生noise variance=No
        pure_y		= zeros(Rx,1228800+5);
        for Tx_n = 1:Tx
            for Rx_n = 1:Rx
                pure_y(Rx_n,:) = pure_y(Rx_n,:) + conv( x_CP(1,:,Tx_n) , H_Channel(:,Rx_n,Tx_n) );
            end
        end
        y=pure_y(:,1:1228800)+n;

        %% receive & remove CP
        y_rmCP		= zeros(2048 , 560,Tx);
        index  = 1;
        for symbol = 1:560
            if mod(symbol,28)-1;
                for Rx_n = 1:Rx
                    y_rmCP(:,symbol,Rx_n) = y(Rx_n,index+144:index+144+2048-1);
                end
                index  = index+144+2048;
            else
                for Rx_n = 1:Rx
                    y_rmCP(:,symbol,Rx_n) = y(Rx_n,index+208:index+208+2048-1);
                end
                index  = index +208+2048;
            end
        end

        %% time to phase
        Y_rmcp=fftshift(fft(y_rmCP/sqrt(2048)),1);

        %% remove guard band
        Y=[Y_rmcp(202+1:202+822,:,:);Y_rmcp(202+822+1+1:202+822+1+822,:,:)];

        %% noise normalize
        X_hat_ZF=zeros(1644,560,Tx);
        norm_Y = zeros(1644,560,2);
		norm_H = zeros(1644,560,2,2);
		norm_Y(:,:,1)   = Y(:,:,1)          ./ sqrt(Rx1_No);
		norm_Y(:,:,2)   = Y(:,:,2)          ./ sqrt(Rx2_No);
		norm_H(:,:,1,:) = H_frame(:,:,1,:)  ./ sqrt(Rx1_No);
		norm_H(:,:,2,:) = H_frame(:,:,2,:)  ./ sqrt(Rx2_No);
        
        % ZF detector (MATLAB code)
 		for SC = 1:1644
 			for slot = 1:560
 				unit_Y 				= reshape(Y(SC,slot,:),Rx,1);
 				unit_H 				= reshape(H_frame(SC,slot,:,:),Tx,Rx);
 				X_hat_ZF(SC,slot,:)	= inv( unit_H'* unit_H ) * unit_H' * unit_Y;
 			end
        end

        % %% ZF detector (C code)
        % X_hat_ZF = ZF_detector(norm_Y,norm_H,Tx);
        
        %% QAM demodulation & dec to bin
        data_hat_ZF = qamdemod(X_hat_ZF,QAM,'gray','UnitAveragePower',true);           %complex to 0~15  ; remember to inverse-normalize
        bin_data_hat = de2bi(data_hat_ZF,q_bit,'left-msb');          %0~15 to binary
        
        %% 算BER
        BER_ZF = BER_ZF+sum(bin_data~=bin_data_hat,'all');

    end

    %% 按照SNR把算好的BER存在矩陣裡
    BER_SNR_ZF(1,a) = BER_ZF/(data_num*q_bit*Tx*1644*560);

end



%輸入SNR_in_dB和BER
figure(1)
semilogy(SNR_in_dB,BER_SNR_ZF(1,:),'r-o','LineWidth',2)
hold on
grid on
title('BER of MIMO')
xlabel('SNR (dB)')
ylabel('BER')
legend('16QAM ZF')

