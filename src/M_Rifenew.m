function [freq] = M_Rifenew(x0,dt)

    sample_time = dt;
    m = length(x0);
    t=1:m;
    flag = 0;
    absfftx0=abs(fft(x0));
    [b1,n1]=max(absfftx0);          % max spectrum line

%-------------------------------
% If n1 is very small, shift right by 1/3 of the sampling frequency
    if n1 < 3                       
        flag = 1;
        x0 = x0.*exp(j*2*pi*(0:m-1)/3);
        g1=exp(j*2*pi*t/(3*m)) ; x1=x0.*g1;
        g2=exp(-j*2*pi*t/(3*m)); x2=x0.*g2;
        absfftx0=abs(fft(x0));%absfftx0(1)=0;absfftx0(m)=0;
        absfftx1=abs(fft(x1));%absfftx1(1)=0;absfftx1(m)=0;
        absfftx2=abs(fft(x2));%absfftx2(1)=0;absfftx2(m)=0;
        %Rife·¨¹À¼Æ²»×ãf
        [b1,n1]=max(absfftx0);
        if n1+1 < m
            b2=absfftx0(n1+1);
        else
            b2=absfftx0(n1);
        end
        if n1-1>1
            b3=absfftx0(n1-1);
        else
            b3=absfftx0(n1);
        end
%If n1 is very large, shift left by 1/3 of the sampling frequency
    elseif n1 > m - 3               
        flag = 2;
        x0 = x0.*exp(-1i*2*pi*(0:m-1)/3);
        g1=exp(1i*2*pi*t/(3*m)) ; x1=x0.*g1;
        g2=exp(-1i*2*pi*t/(3*m)); x2=x0.*g2;
        absfftx0=abs(fft(x0));%absfftx0(1)=0;absfftx0(m)=0;
        absfftx1=abs(fft(x1));%absfftx1(1)=0;absfftx1(m)=0;
        absfftx2=abs(fft(x2));%absfftx2(1)=0;absfftx2(m)=0;

        [b1,n1]=max(absfftx0);
        if n1+1 < m
            b2=absfftx0(n1+1);
        else
            b2=absfftx0(n1);
        end
        if n1-1>1
            b3=absfftx0(n1-1);
        else
            b3=absfftx0(n1);
        end
    else                        
        g1=exp(1i*2*pi*t/(3*m)) ; x1=x0.*g1; % right shift the signal
        g2=exp(-1i*2*pi*t/(3*m)); x2=x0.*g2; % left shift the signal
        absfftx0(1)=0;absfftx0(m)=0;
        absfftx1=abs(fft(x1));%absfftx1(1)=0;absfftx1(m)=0;
        absfftx2=abs(fft(x2));%absfftx2(1)=0;absfftx2(m)=0;

    if n1+1 < m
            b2=absfftx0(n1+1);
        else
            b2=absfftx0(n1);
        end
        if n1-1>1
            b3=absfftx0(n1-1);
        else
            b3=absfftx0(n1);
        end
    end

    if   b2>=b3
        fx0=(n1-1+1*b2/(b1+b2))/(m*sample_time);
    else
        fx0=(n1-1-1*b3/(b1+b3))/(m*sample_time);
    end

    [b1,n1]=max(absfftx1);
    if n1+1 < m
        b2=absfftx1(n1+1);
    else
        b2=absfftx1(n1);
    end


    if n1-1>1
        b3=absfftx1(n1-1);
    else
        b3=absfftx1(n1);
    end
    if   b2>=b3
        fx1=(n1-1+1*b2/(b1+b2))/(m*sample_time);
    else
        fx1=(n1-1-1*b3/(b1+b3))/(m*sample_time);
    end

    [b1,n1]=max(absfftx2);
    if n1+1 < m
        b2=absfftx2(n1+1);
    else
        b2=absfftx2(n1);
    end
    if n1-1>1
        b3=absfftx2(n1-1);
    else
        b3=absfftx2(n1);
    end
    if   b2>=b3
        fx2=(n1-1+1*b2/(b1+b2))/(m*sample_time);
    else
        fx2=(n1-1-1*b3/(b1+b3))/(m*sample_time);
    end

    divide=1/(m*sample_time);%
    distx0=abs(rem(fx0,divide)-0.5*divide);
    distx1=abs(rem(fx1,divide)-0.5*divide);
    distx2=abs(rem(fx2,divide)-0.5*divide);
    temp=[distx0,distx1,distx2];


    [~,pos]=min(temp);
    if pos==1
        freq=fx0;
    elseif pos==2
        freq=fx1-1/(3*m*dt);
    else
        freq=fx2+1/(3*m*dt) ;
    end

    if flag ==1
        freq = freq - 1/dt/3;
    elseif flag == 2
        freq = freq + 1/dt/3;
    end

Fs = 1/dt;
    if freq>1/2*Fs
       freq = Fs- freq;
    end
    return;



