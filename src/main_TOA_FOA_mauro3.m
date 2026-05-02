clear all
N=100 %number of montecarlo trials
address = '4840D6';
message_hex = '202CC371C32CE0';

Fs=500e6;             % Sampling frequency
Fm=20e6;              %signal frequency
SNR = [10 15,20,30,40]
%SNR=20
numbit =4;
adsamlingrate=[100,50, 25,10]; %downsampling rates
%adsamlingrate=10
flagplot=0;

for i=1:length(SNR)
    for j=1:length(adsamlingrate)
    
[rmsTOA(i,j),stdTOA(i,j),meanTOA(i,j),rmsTOA2(i,j),stdTOA2(i,j),meanTOA2(i,j), rmsTOARife(i,j),stdTOARife(i,j),meanTOARife(i,j),...
    rmsTOAth(i,j),stdTOAth(i,j),meanTOAth(i,j),rmsTOAthint(i,j),stdTOAthint(i,j),meanTOAthint(i,j),...
    rmsFOAFFT(i,j),stdFOAFFT(i,j),meanFOAFFT(i,j),rmsFOARife(i,j),stdFOARife(i,j),meanFOARife(i,j),rmsFOARifeM(i,j),stdFOARifeM(i,j),meanFOARifeM(i,j)] =...
    MontecarloTOAFOA3(N,address, message_hex,Fs,Fm, numbit,SNR(i),adsamlingrate(j),flagplot);
    end
i
end


%%   Fa le figure
xaxis=SNR;
yaxis=Fs./adsamlingrate*0.4/10^6;
figure
surf(xaxis,yaxis,stdTOA')
title('std error TOA (Correlation Algorithm)')
xlabel('SNR (dB)') 
ylabel('Receiver Band (Mhz)') 
zlabel('std error (sec)')
set(gca, 'ZScale', 'log')
fontsize(gca,24,"points")


figure
surf(xaxis,yaxis,stdTOA2')
title('std error TOA (Correlation + Interpolation)')
xlabel('SNR (dB)') 
ylabel('Receiver Band (Mhz)') 
zlabel('std (sec)')
set(gca, 'ZScale', 'log')
fontsize(gca,24,"points")


figure
surf(xaxis,yaxis,stdTOARife')
title('std error TOA (Correlation  + Rife)')
xlabel('SNR (dB)') 
ylabel('Receiver Band (Mhz)') 
zlabel('std (sec)')
set(gca, 'ZScale', 'log')
fontsize(gca,24,"points")




figure
surf(xaxis,yaxis,stdTOAth')
title('std error TOA (Rising Edge Detection)')
xlabel('SNR (dB)') 
ylabel('Receiver Band (Mhz)') 
zlabel('std (sec)')
set(gca, 'ZScale', 'log')
fontsize(gca,24,"points")




figure
surf(xaxis,yaxis,stdTOAthint')
title('std error TOA (Rising Edge detection + Interpolation)')
xlabel('SNR (dB)') 
ylabel('Receiver Band (Mhz)') 
zlabel('std (sec)')
set(gca, 'ZScale', 'log')
fontsize(gca,24,"points")



figure
surf(xaxis,yaxis,stdFOAFFT')
title('std error FOA (FFT)')
xlabel('SNR (dB)') 
ylabel('Receiver Band (Mhz)') 
zlabel('std (sec)')
set(gca, 'ZScale', 'log')
fontsize(gca,24,"points")


figure
surf(xaxis,yaxis,stdFOARife')
title('std error FOA (FFT+Rife)')
xlabel('SNR (dB)') 
ylabel('Receiver Band (Mhz)') 
zlabel('std (sec)')
set(gca, 'ZScale', 'log')
fontsize(gca,24,"points")

figure
surf(xaxis,yaxis,stdFOARifeM')
title('std error FOA (FFT+ Modiefied Rife)')
xlabel('SNR (dB)') 
ylabel('Receiver Band (Mhz)') 
zlabel('std (sec)')
set(gca, 'ZScale', 'log')
fontsize(gca,24,"points")

% % % % scrive le std
% % % writematrix(stdFOAFFT, 'stdFOAFFT_4bit.xls');
% % % writematrix(stdFOARife, 'stdFOARife_4bit.xls');
% % % writematrix(stdFOARifeM, 'stdFOAmRife_4bit.xls');
% % % 
% % % 
% % % writematrix(stdTOAth, 'stdTOA_th_4bit.xls');
% % % writematrix(stdTOAthint, 'stdTOA_th_int_4bit.xls');
% % % writematrix(stdTOARife, 'stdTOA_corr_Rife_4bit.xls');
% % % writematrix(stdTOA, 'stdTOA_corr_4bit.xls');
% % % writematrix(stdTOA2, 'stdTOA_corr_int_4bit.xls');
% % % 
% % % %scrive le rms
% % % 
% % % writematrix(rmsFOAFFT, 'rmsFOAFFT_4bit.xls');
% % % writematrix(rmsFOARife, 'rmsFOARife_4bit.xls');
% % % writematrix(rmsFOARifeM, 'rmsFOAmRife_4bit.xls');
% % % 
% % % 
% % % writematrix(rmsTOAth, 'rmsTOA_th_4bit.xls');
% % % writematrix(rmsTOAthint, 'rmsTOA_th_int_4bit.xls');
% % % writematrix(rmsTOARife, 'rmsTOA_corr_Rife_4bit.xls');
% % % writematrix(rmsTOA, 'rmsTOA_corr_4bit.xls');
% % % writematrix(rmsTOA2, 'rmsTOA_corr_int_4bit.xls');
% % % 
% % % % scrive le mean
% % % 
% % % writematrix(meanFOAFFT, 'meanFOAFFT_4bit.xls');
% % % writematrix(meanFOARife, 'meanFOARife_4bit.xls');
% % % writematrix(meanFOARifeM, 'meanFOAmRife_4bit.xls');
% % % 
% % % 
% % % writematrix(meanTOAth, 'meanTOA_th_4bit.xls');
% % % writematrix(meanTOAthint, 'meanTOA_th_int_4bit.xls');
% % % writematrix(meanTOARife, 'meanTOA_corr_Rife_4bit.xls');
% % % writematrix(meanTOA, 'meanTOA_corr_4bit.xls');
% % % writematrix(meanTOA2, 'meanTOA_corr_int_4bit.xls');


%%
save("results\nuovi_results_4bit_100prove")
