clear
clc
close all

%% setting
MIMO_range= 1:6; %天線數量
frame_num = 50;
SNR_in_dB = 0:5:40;        % 自己設訂雜訊大小
C_in_SNR  = zeros(length(MIMO_range),length(SNR_in_dB));

%% antenna 1x1 to 6x6
for smui=1:length(MIMO_range)
	Tx	= MIMO_range(smui);
	Rx	= Tx;

    %% SNR 0:5:40
	for a=1:length(SNR_in_dB)
		SNR = 10^( SNR_in_dB(a)/10);
		No  = 10^(-SNR_in_dB(a)/10);
		C 	= 0;

        %% frames / SNR
		for frame = 1:frame_num	
			fprintf("MIMO : %d/%d\t SNR : %d/%d \t frame : %d/%d\n",smui,length(MIMO_range),a,length(SNR_in_dB),frame,frame_num);

			%% 通道
			PowerdB 	= [ -2 -8 -10 -12 -15 -18].';
			PowerdB_MIMO= repmat(PowerdB,1,Rx,Tx);
			H_Channel 	= sqrt(10.^(PowerdB_MIMO./10));
			Ntap		= length(PowerdB);
			H_Channel   = H_Channel .* ( sqrt( 1/(2*Tx) ) .* ( randn(Ntap,Rx,Tx) + 1i*randn(Ntap,Rx,Tx) ) );
			h = [H_Channel ; zeros(2042,Rx,Tx)];
			H = fftshift(fft(h,[],1),1);
			H_Data 	= [H(203:1024,:,:);H(1026:1847,:,:)];	%rmGB
			H_frame = permute(repmat( H_Data(:,:,:),1,1,1,560),[1 4 2 3]  );

            %% calculate capacity
			for SC = 1:1644
				for slot = 1:560
					unit_H 	= reshape(H_frame(SC,slot,:,:),Rx,Tx);
					C		= C + abs( log2( det( eye(Tx) + (1/Tx) .* SNR .* (unit_H*unit_H'))));
				end
			end
		end
		C_in_SNR(smui,a) = C / (1644 * 560 * frame_num);
    end
end

%% pic
figure(1)
plot(SNR_in_dB,C_in_SNR(1,:),'r-','LineWidth',2);
hold on;
plot(SNR_in_dB,C_in_SNR(2,:),'g-','LineWidth',2);
hold on;
plot(SNR_in_dB,C_in_SNR(3,:),'b-','LineWidth',2);
hold on;
plot(SNR_in_dB,C_in_SNR(4,:),'c-','LineWidth',2);
hold on;
plot(SNR_in_dB,C_in_SNR(5,:),'m-','LineWidth',2);
hold on;
plot(SNR_in_dB,C_in_SNR(6,:),'k-','LineWidth',2);
hold on;
grid on;
title('Capacity of MIMO OFDM');
xlabel('SNR (dB)');
ylabel('Capacity(bits/sec/Hz)');
legend('1x1','2x2','3x3','4x4','5x5','6x6');