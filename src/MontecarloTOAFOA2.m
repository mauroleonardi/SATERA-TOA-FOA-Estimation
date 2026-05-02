function [rmsTOA,stdTOA,meanTOA,rmsTOA2p,stdTOA2p,meanTOA2p, rmsTOARife,stdTOARife,meanTOARife,...
          rmsTOAth,stdTOAth,meanTOAth,rmsTOAthint,stdTOAthint,meanTOAthint,...
          rmsFOAFFT,stdFOAFFT,meanFOAFFT,rmsFOARife,stdFOARife,meanFOARife,rmsFOARifeM,stdFOARifeM,meanFOARifeM] = MontecarloTOAFOA2(N,address, message_hex,Fs,Fm, numbit,SNR,adsamlingrate,flagplot)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%N = 100;             % MONTE CARLO
FoEstRMSEfreq = zeros(1,N); % RMSE Initialization
FoEstSTDfreq = ones(1,N); % Standard deviation Initialization
FoEstMeanfreq = zeros(1,N); % BIAS Initialization
FoEstRMSEfreqrife = zeros(1,N);
FoEstMeanfreqrife = zeros(1,N);
FoEstSTDfreqrife = ones(1,N);
FoEstRMSEfreqmrife = zeros(1,N);
FoEstMeanfreqmrife = zeros(1,N);
FoEstSTDfreqmrife = ones(1,N);
FoEstRMSETOA = zeros(1,N);
FoEstMeanTOA = zeros(1,N);
FoEstSTDTOA = ones(1,N);

rng("shuffle")
for j= 1:N
    %ADS-B message
    %TYPE and ADS-B emitter category
    df = 17 ;
    ca = 5 ;
    %address = '4840D6';
    %message_hex = '202CC371C32CE0';
    %M=5;
    adsb.df = df ;
    adsb.ca = ca ;
    adsb.address = address;
    adsb.message_hex = message_hex;

    % Hex to Binary
    adsb.message_bin = [
        int32(dec2bin(hex2dec(adsb.message_hex(1:7)),28))-'0' ...
        int32(dec2bin(hex2dec(adsb.message_hex(8:end)),28))-'0'
        ] ;
    adsb.df = int32(dec2bin(adsb.df,5))-'0';
    adsb.ca = int32(dec2bin(adsb.ca,3))-'0';
    adsb.address = int32(dec2bin(hex2dec(adsb.address),24))-'0';

    % Assemble
    adsb.payload = [ adsb.df adsb.ca adsb.address adsb.message_bin ];

    % CRC
    adsb.crc24 = comm.CRCGenerator([1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 0 1]) ;
    adsb.payload = step(adsb.crc24, logical(adsb.payload)' )';

    % PPM modulation
    adsb.ppm.zero = [ 0 1 ] ;
    adsb.ppm.one = [ 1 0 ] ;
    for x=0:length(adsb.payload)-1
        if adsb.payload(x+1) == 0
            adsb.encoded(2*x+1:2*x+2) = adsb.ppm.zero ;
        else
            adsb.encoded(2*x+1:2*x+2) = adsb.ppm.one ;
        end
    end
    adsb.encoded; % 224

    % add preamble
    adsb.preamble = [ 1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 ] ;
    adsb.tx = [ adsb.preamble adsb.encoded ] ; % 240


    %%
    %Fs=300e6;             % Sampling frequency
    %Fm=20e6;              % carrier frequency
    P0 = 0:1/(Fs/(10e6)):(1-1/(Fs/(10e6))); % Pulse Rise points
    P1 = ones(1,4*Fs/(10e6));               % pulse mid points
    P2 = (1-1/(Fs/(10e6))):-1/(Fs/(10e6)):0;% Pulse Decay points
    Pulse1 = [P0,P1,P2];
    %Pulse0 = zeros(1,Fs*0.4*10^6);
    adsb.txup =zeros(1,round(120*10^(-6)*Fs));
    for i = 1:length(adsb.tx)
        if adsb.tx(i)==1
           adsb.txup(round(1+(i-1)*0.5*10^(-6)*Fs):round((i-1)*0.5*10^(-6)*Fs+length(Pulse1)))= Pulse1;
        %else
        %    adsb.txup = [adsb.txup,Pulse0];
        end
    end

    for i =  8.75*10^-6*Fs:Fs*10^-6*0.5:Fs*10^-6*0.5*239
        if adsb.txup(i)==1 & adsb.txup(i-Fs*10^-6*0.5)==1
           adsb.txup((i-1*Fs*10^-6*0.5):(i))= 1;
        %else
        %    adsb.txup = [adsb.txup,Pulse0];
        end
    end

    expreambleFs=adsb.txup(1:Fs*10^-6*8+1);
    %%


    %SNR = 10;
    %numbit =10; % quantization bits
    %adsamlingrate=50;

    %%

    %%



    FoRange = -1e6:1:1e6;
    Fo= FoRange(randi(length(FoRange))); %frequency offset

   %% mdificato da mauro per miglirare la granularità dello shift in
   % frequenza
   Fo=2*rand(1)*10^6-10^6;

   %% fine modifica 


    p= rand; %random phase
    
    %tau1=randi(Fs/10e4/adsamlingrate)+30;
    %tau = adsamlingrate*tau1;
    %tau2= 10*Fs/10e5-tau;
    %tau= 1500;

    %%% Migliorato da mauro per aggiungere granularità al TOA
    tau=round(10/10^6*Fs+rand(1)*1*10^(-6)*Fs); %AGGIUNGE 10 us + rand di 1 us
   % tau=round(10/10^6*Fs/adsamlingrate)*adsamlingrate;
    
    
    
    tau2= round(10/10^6*Fs); %aggiunge un us alla fine

    %%% fine miglioramento

    datalength = length(adsb.txup);
    t = (0:1/Fs:(datalength-1)/Fs);
    adsb.txup = adsb.txup.*cos(2*pi*(Fm+Fo)*t+2*pi*p);
    adsb.txupzeros=[zeros(1,tau) adsb.txup zeros(1,tau2)];
    adsb_noise = awgn(double(adsb.txupzeros),SNR) ;
    %%

    xr = adsb_noise;
    datalength = length(xr);
    t = 0:1/Fs:(datalength-1)/(Fs);

    carrierI=2*cos(2*pi*Fm*t);
    carrierQ=-2*sin(2*pi*Fm*t);

    receiver_I=carrierI.*xr;
    ffs = 2e6; % cut-off frequecy
    receiver_I_Lowpass = decimate(receiver_I,adsamlingrate);
    %figure
    %plot(t*10e5,receiver_I_Lowpass);title('in-phase component');
    %ylim([0 1.5])
    %xlim([0,20])

    receiver_Q=carrierQ.*xr;
    receiver_Q_Lowpass = decimate(receiver_Q,adsamlingrate);
    %figure
    %plot(t*10e5,receiver_Q_Lowpass);title('quadrature component');hold on
    %ylim([0 1])
    %xlim([0,20])


    A_s = abs(receiver_I_Lowpass+receiver_Q_Lowpass*1i); % envelope




    %%
    %figure
    %plot(t*10e5,abs(receiver_I_Lowpass+receiver_Q_Lowpass*1i));title("ADS-B signal");hold on
    %ylim([0 1.5])
    %xlim([0,20])

    ad_Fs= Fs/adsamlingrate; % sampling frequency
    xad=A_s; % input
    datalength = length(xad);
    t = 0:1/ad_Fs:(datalength-1)/ad_Fs;
    %figure
    %plot(t*10e5,xad);title('base band signal');grid on;        % base band signal
    %ylim([0 1.5])
    %xlim([0,20])


    maxnoise=rms(1000*xad(1:round((tau-300)/(Fs/ad_Fs))))/1000;
    ad_signal= xad/maxnoise;

    quan_signal_env=round(ad_signal);
    for i = 1:length(quan_signal_env)
        if quan_signal_env(i) > 2^numbit
            quan_signal_env(i) = 2^numbit;
        end
    end
    %figure
    %plot(t*10e5,quan_signal_env);title('quantized signal');hold on%
    %ylim([0 10])
    %xlim([0,20])

    %%
    xadi=receiver_I_Lowpass; % sampling
    datalength = length(xadi);
    t = 0:1/ad_Fs:(datalength-1)/ad_Fs;


    %maxnoise=max(1000*xadi(1:tau/(Fs/ad_Fs)))/1000;
    ad_signal= xadi/maxnoise;

    quan_signal_I=round(ad_signal);
    for i = 1:length(quan_signal_I)
        if quan_signal_I(i) > 2^(numbit-1)
            quan_signal_I(i) = 2^(numbit-1);
        elseif quan_signal_I(i) < -2^(numbit-1)
            quan_signal_I(i) = -2^(numbit-1);
        end
    end
    %figure
    %plot(t*10e5,quan_signal_I);title('I quantized signal');hold on%
    %ylim([0 10])
    %xlim([0,20])

    %%

    xadq=receiver_Q_Lowpass; % sampling
    datalength = length(xadq);
    t = 0:1/ad_Fs:(datalength-1)/ad_Fs;

    %maxnoise=max(xadq(1:tau/(Fs/ad_Fs)));
    ad_signal= xadq/maxnoise;

    quan_signal_Q=round(ad_signal);
    for i = 1:length(quan_signal_Q)
        if quan_signal_Q(i) > 2^(numbit)-1
            quan_signal_Q(i) = 2^(numbit)-1;
        elseif quan_signal_Q(i) < -2^(numbit)-1
            quan_signal_Q(i) = -2^(numbit)-1;
        end
    end
    %figure
    %plot(t*10e5,quan_signal_Q);title('Q quantized signal');hold on%
    %ylim([0 10])
    %xlim([0,20])

    % xr = quan_signal_Q;
    % datalength = length(xr);
    % t = 0:10/Fs:(datalength-1)/(Fs/10);
    % plot(t*10e5,quan_signal_env)
    % hold on
    % ylim([0 8])
    % xlim([0,130])
    % xlabel('Time(\mu s)')
    % ylabel('Amplitude')
    % title('Envelope of the signal after A/D@3bits,30MHz.')

    %%
    Ff = ad_Fs/3;


    %%%%%% qui ci può essere un errore 
    x0 = quan_signal_I.*cos(2*pi*Ff*t)-quan_signal_Q.*sin(2*pi*Ff*t); % bandpass signal
    x03 = quan_signal_I.*cos(2*pi*Ff*t);%-quan_signal_Q.*sin(2*pi*Ff*t); % bandpass signal
    %x0 = receiver_I_Lowpass.*cos(2*pi*Ff*t)-receiver_Q_Lowpass.*sin(2*pi*Ff*t);
    x02 = quan_signal_I+1i*quan_signal_Q;

    dt = 1/ad_Fs;
    absfftx0=abs(fft(x0));
    absfftx02=abs(fft(x02));
   
    %absfftx03=abs(fft(x03));
%     figure
%     plot(absfftx0)
%     hold on
%     plot(absfftx02)
    %plot(absfftx03)
   
    %Y = fftshift(absfftx0);
    %fhift = (-length(absfftx0)/2:length(absfftx0)/2-1)*(ad_Fs/length(absfftx0));
    [b1,n1]=max(absfftx0);
    [b2,n02]=max(absfftx02);
    if n02 > length(absfftx0)/2
         n022=n02-length(absfftx02);
         f02=(n022-1)*ad_Fs/length(absfftx02); 
    else
        f02=(n02-1)*ad_Fs/length(absfftx02); 
    end
    % if n1 > length(absfftx0)/2
    %     n1=n1-length(absfftx0)-1;
    % end
    f=(n1-1)*ad_Fs/length(absfftx0); % FFT-based estimated frequency
  %  f02=(n02-1)*ad_Fs/length(absfftx02); % FFT-based estimated frequency
    


    e= f - Fo-Ff;% error
    e2= f02 - Fo;% error
    E(j)=e;
   % eeee(j)=e;
   % eeee2(j)=e2;

%%% RIFE Tecnhiques
    
    if absfftx0(n1+1) < absfftx0(n1-1)
        n2= n1 - 1;
        r = -1;
    else
        n2= n1 + 1;
        r = 1;
    end

    delta_n =  r*absfftx0(n2)/(absfftx0(n1)+absfftx0(n2));

    k= n1+delta_n-1;
    fe=k*ad_Fs/length(absfftx0);
    ee= (fe - Fo-Ff); % Rife error
    EE(j)=ee;

%     %%% RIFE su BB da mauro
%  
%     if absfftx02(n02+1) < absfftx02(n02-1)
%         n2= n02 - 1;
%         r = -1;
%     else
%         n2= n02 + 1;
%         r = 1;
%     end
% 
%     delta_n2 =  r*absfftx02(n2)/(absfftx02(n02)+absfftx02(n2));
%     if n02 > length(absfftx0)/2
%         k2= n022+delta_n2-1;
%     else
%         k2= n02+delta_n2-1;
%     end
%     fe2=k2*ad_Fs/length(absfftx02);
%     ee2= (fe2 - Fo); % Rife error Mauro
%     EE2(j)=ee2;



    %%% Modfied rife
    freq = M_Rifenew(x0,1/ad_Fs);
    eee= (freq-Fo-Ff); % M-Rife error
    1;
    % %%
    EEE(j)=eee;

    adsbSig=quan_signal_env;
    adsbSig = adsbSig/max(adsbSig);

    %figure
    %plot(adsbSig);
    % standard preamble
    adsbParam.SyncSequence = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];
    adsbParam.SyncSequenceLength = length(adsbParam.SyncSequence);
    adsbParam.SyncSequenceHighIndices = find(adsbParam.SyncSequence);
    adsbParam.SyncSequenceNumHighValues = length(adsbParam.SyncSequenceHighIndices);
    adsbParam.SyncSequenceLowIndices = find(~adsbParam.SyncSequence);
    adsbParam.SyncSequenceNumLowValues = length(adsbParam.SyncSequenceLowIndices);
    % with sampling frequency Fs，number of samples of each pulse: Fs/(2e6)
    adsbParam.SamplesPerSymbol = ad_Fs/(2e6);
    % extend the preamble sequence
   % syncSignal = reshape(ones(adsbParam.SamplesPerSymbol,1)...
    %    *adsbParam.SyncSequence, 16*adsbParam.SamplesPerSymbol, 1);
    % turn 0 to -1: [1 -1 1 -1 -1 -1 -1 1 -1 1 -1 -1 -1 -1 -1 -1]
    %adsbParam.SyncFilter = single(2*syncSignal-1);



 %% aggiunto da mauro per migliorare la funzione attesa per la correlazione
    expreamble = [ 1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 ]; 
    %expreamble2 = [ 1 0 1 0];% 0 0 0 1 0 1 0 ]; 
     %expreambleFs=[];

  %   expreambleFs=adsb.txup(1:Fs*10^-6*8+1)
%     for i = 1:length(expreamble)
%         if expreamble(i)==1
%             expreambleFs = [expreambleFs,Pulse1];
%         else
%             expreambleFs = [expreambleFs,Pulse0];
%         end
%     end
%    expreambleFs(expreambleFs==0)=-1;
    
    expreambleFad=decimate(double(expreambleFs),adsamlingrate);
    expreambleFad(expreambleFad<0.5)=-1;
    expreambleFad(expreambleFad>0.5)=1;
%%

    % Cross-correlation calculation
    %xFilt = xcorr(adsbSig,adsbParam.SyncFilter);

    xFilt= xcorr(abs(x02)/max(abs(x02)),expreambleFad);

    xFilt2= xcorr(abs(x02)/max(abs(x02)),expreambleFad);

    %xFilt3= xcorr(resample(abs(x02)/max(abs(x02)),100,1),resample(expreambleFad,100,1));
    %xFilt3= xcorr(resample(adsbSig/max(adsbSig),100,1),resample(double(adsbParam.SyncFilter),100,1));


%% TOA con soglia
    
    %ciccio=sign(dif(abs(x02)))
    Edge=diff(sign(abs(x02)-4));
    RiseEdge= find(Edge>0);
    FallEdge= find(Edge<0);
    TOA_tH=RiseEdge(1);
    TOA_tH_err=(TOA_tH-tau/adsamlingrate)/ad_Fs;
%    TOA_tH2=RiseEdge(1)+FallEdge(1))/2-0.25*10;

%% TOA con soglia su salita e discesa
    
    %ciccio=sign(dif(abs(x02)))
    %Edge=diff(sign(abs(x02)-4));
    %RiseEdge= find(Edge>1);
    %FallEdge= find(Edge<-1);
    TOA_tHa=(RiseEdge(1)+FallEdge(1))/2-0.26*10^-6*ad_Fs;
    TOA_tHa_err=(TOA_tHa-tau/adsamlingrate)/ad_Fs;
%    TOA_tH2=RiseEdge(1)+FallEdge(1))/2-0.25*10;


%% TOA con soglia su segnale interpolato
    uprate=10;
    Edge2=diff(sign(abs(interp(x02,uprate))-4));
    RiseEdge2= find(Edge2>1);
    FallEdge2= find(Edge2<-1);
    TOA_tH2=RiseEdge2(1)/uprate;
    TOA_tH2_err=(TOA_tH2-tau/adsamlingrate)/ad_Fs;



    %% TOA con correlatore 


    %figure
    %plot(xFilt);
    xxFilt=xFilt(1:round(length(adsbSig)*1.5));
    [a,b] = max(xxFilt);

    
    xxFilt2=xFilt2(1:round(length(adsbSig)*1.5));
    [a2,b2] = max(xxFilt2);
%    b3=b-101+b2
 %   b2=b3


    %xxFilt3=xFilt3(1:round(length(adsbSig)*1.5*100));
    %[a3,b3] = max(xxFilt3);
    %b3=b3/100;

    %datatip(h,b,xFilt(b))
    % TOA
    siganlStart = b-length(adsbSig);
    TOAe=(siganlStart-tau/adsamlingrate)/ad_Fs;


    %siganlStart3p = b3-length(adsbSig)+1;
    %TOAe3p=(siganlStart3p-tau/adsamlingrate-1)/ad_Fs;

    siganlStart2p = b2-length(adsbSig);
    TOAe2p=(siganlStart2p-tau/adsamlingrate)/ad_Fs;


%    Seleziono su quale fare Rifle
    b=b2;
    xxFilt=xxFilt2;

%%% TOA con correlatore + Rifle

%%% Prova Rife sui tempi aggiunto da mauro
 if xxFilt(b+1) < xxFilt(b-1)
        bbb2= b - 1;
        r = -1;
    else
        bbb2= b + 1;
        r = 1;
    end



    delta_nT =  r*xxFilt(bbb2)/(xxFilt(b)+xxFilt(bbb2));
%% Toa correlatore + interpolatore

%%% prova interpolazione + derivatore
diff1= xxFilt(b+1)- xxFilt(b);
diff2=  xxFilt(b)- xxFilt(b-1);
delta=diff2/(diff1-diff2);
kTM=b+delta;
 siganlStartM = kTM-length(adsbSig);
    TOAeM=(siganlStartM-tau/adsamlingrate)/ad_Fs;

    kT= b+delta_nT;
    siganlStartRife = kT-length(adsbSig);
    TOAeRIFE=(siganlStartRife-tau/adsamlingrate)/ad_Fs;     % TOA Rife error
%%%%% fine prova
    TTRife(j)=TOAeRIFE; %Rifle
    TT2(j)=TOAeM;       %interpolatore
    TT(j)=TOAe2p;   %secco
    TTth(j)=TOA_tH_err;
    TTth_int(j)=TOA_tH2_err;
    %title("TOA"+num2str(siganlStart))
    %%
%     FoEstRMSEfreq(j) = (e).^2;
%     FoEstMeanfreq(j) = (e);
%     FoEstRMSEfreqrife(j) = (ee).^2;
%     FoEstMeanfreqrife(j) = (ee);
%     FoEstRMSEfreqmrife(j) = (eee).^2;
%     FoEstMeanfreqmrife(j) = (eee);
%     FoEstRMSETOA(j) = (TOAe).^2;
%     FoEstMeanTOA(j) = (TOAe);

%% other measurement
Edge2=diff(sign(abs(interp(x02,uprate))-4));
RiseEdge2= find(Edge2>1);
FallEdge2= find(Edge2<-1);

sig_ampl=mean(abs(x02(find(abs(x02)>4))))*maxnoise;
sig_ampl_e(j)=sig_ampl-1+maxnoise;
a=a;

% RMSEfreqfft = sqrt(sum(FoEstRMSEfreq,2)/N);
% Meanfreqfft = sum(FoEstMeanfreq,2)/N;
% for i=1:length(FoEstMeanfreq)
%     FoEstSTDfreq(i)  = (FoEstMeanfreq(i)-Meanfreqfft).^2;
% end
% STDfreqfft = sqrt(sum(FoEstSTDfreq,2)/N);
% Meanfreqerror = sum(abs(FoEstMeanfreq),2)/N;
% 
% RMSEfreqrife = sqrt(sum(FoEstRMSEfreqrife,2)/N);
% Meanfreqrife = sum(FoEstMeanfreqrife,2)/N;
% for i=1:length(FoEstMeanfreqrife)
%     FoEstSTDfreqrife(i)  = (FoEstMeanfreqrife(i)-Meanfreqrife).^2;
% end
% STDfreqrife = sqrt(sum(FoEstSTDfreqrife,2)/N);
% 
% RMSEfreqmrife = sqrt(sum(FoEstRMSEfreqmrife,2)/N);
% Meanfreqmrife = sum(FoEstMeanfreqmrife,2)/N;
% for i=1:length(FoEstMeanfreqmrife)
%     FoEstSTDfreqmrife(i)  = (FoEstMeanfreqmrife(i)-Meanfreqmrife).^2;
% end
% STDfreqmrife = sqrt(sum(FoEstSTDfreqmrife,2)/N);
% 
% 
% 
% 
% RMSETOA = sqrt(sum(FoEstRMSETOA,2)/N);
% MeanTOA = sum(FoEstMeanTOA,2)/N;
% for i=1:length(FoEstMeanTOA)
%     FoEstSTDTOA(i)  = (FoEstMeanTOA(i)-MeanTOA).^2;
% end
% STDTOA = sqrt(sum(FoEstSTDTOA,2)/N);
% MeanTOAerror = sum(abs(FoEstMeanTOA),2)/N;
% 

% if flagplot
%     figure
%     histogram('Categories',{'RMSE','Mean','STD'},'BinCounts',[RMSEfreqfft Meanfreqerror STDfreqfft],'BarWidth',0.3)
%     figure
%     histogram('Categories',{'RMSE','Mean','STD'},'BinCounts',[RMSETOA MeanTOAerror STDTOA],'BarWidth',0.3)
% 
%     hist(FoEstMeanTOA,10);title('Errors of TOA in 1000 Monte Carlo.');hold on
%     xlabel('second')
% end

% rmsTOA=RMSETOA;
% stdTOA=STDTOA;
% meanTOA=MeanTOAerror;
% 

end

rmsTOARife=rms(TTRife);
stdTOARife=std(TTRife);
meanTOARife=mean(TTRife);

rmsTOA=rms(TT);
stdTOA=std(TT);
meanTOA=mean(TT);


rmsTOA2p=rms(TT2);
stdTOA2p=std(TT2);
meanTOA2p=mean(TT2);

rmsTOAth=rms(TTth);
stdTOAth=std(TTth);
meanTOAth=mean(TTth);

rmsTOAthint=rms(TTth_int);
stdTOAthint=std(TTth_int);
meanTOAthint=mean(TTth_int);


rmsFOAFFT=rms(E);
stdFOAFFT=std(E);
meanFOAFFT=mean(E);


rmsFOARife=rms(EE);
stdFOARife=std(EE);
meanFOARife=mean(EE);


rmsFOARifeM=rms(EEE);
stdFOARifeM=std(EEE);
meanFOARifeM=mean(EEE);



end