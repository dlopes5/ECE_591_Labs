% ECE 591 
% Software Defined Radio
% Midterm Project

clc;clear;close all

%% Part 1
user_input=input('Enter the string you would like to send: ', 's');
user_input_ascii= int2bit(double(user_input),7);
user_input_ascii= reshape(user_input_ascii,[1,length(user_input)*7]);



fs = 96000; % sampling frequency
Ts = 1/fs; % sampling duration
symbolrate = 1000; % transmitted pulses/second should be an integer divisor of fs
sps = fs/symbolrate; % number of samples in one symbol
tau=1/sps;
delay=zeros(1,round(tau*fs));

%% 2 PAM
% Transmitter setup
% define the basic rectangular pulse shape

pulse = ones(1,sps); % 1 pulse
pulseham= hamming(sps);%creates hamming pulse
pulsehann=hann(sps);

%delayedPulse=[delay, pulse];
data = user_input_ascii;
%data = [1 0 1 1 0 0 1 1 1 0 0 0 1 1 0 1 1 0];

time=linspace(-6,6,sps);
truncate_sinc=sinc(time);


[number, length0] = size(data);
symbol = 2*data-1;
for i = 1: length0
    for j=1:sps
        signal((i-1)*sps + j) = symbol(i) * pulse(j);
        signalHamm((i-1)*sps + j) = symbol(i) * pulseham(j);
        signalHann((i-1)*sps + j) = symbol(i) * pulsehann(j);
        %signalSinc((i-1)*sps + j) = symbol(i) * truncate_sinc(j);
    end
end

%use the truncated sinc to mod the signal
signalSinc= sinc_mod(symbol,sps);
plot(signalSinc)

%create delayed signal
delayedSignal=[delay, signal];
delayindex= 0:Ts:(length0*sps-1)*Ts+tau;

index = 0: Ts: (length0*sps-1)*Ts;
freq=fft(signal);
freq=fftshift(freq);

freqHamm=fft(signalHamm);
freqHamm=fftshift(freqHamm);

freqHann=fft(signalHann);
freqHann=fftshift(freqHann);

freqSinc=fft(signalSinc);
freqSinc=fftshift(freqSinc);

freqSincTrunk=fft(truncate_sinc);
freqSincTrunk=fftshift(freqSincTrunk);

%
figure();
plot(index, signal);
axis([0 1.5*length0*sps*Ts -2 2]);
xlabel('time');
ylabel('amplitude');
title('basband transmission signal of 2PAM');
figure();
plot(abs(freq));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of 2PAM signal');
figure()
plot(delayindex, delayedSignal)
title('2 Pam Signal with delay')
xlabel('time(sec)')
ylabel('symbol')
figure()
plot(index, signalHamm);
axis([0 1.5*length0*sps*Ts -2 2]);
xlabel('time');
ylabel('amplitude');
title('basband transmission signal of 2PAM (Hamming window)');
figure()
plot(index, signalHann);
axis([0 1.5*length0*sps*Ts -2 2]);
xlabel('time');
ylabel('amplitude');
title('basband transmission signal of 2PAM (Hanning window)');

figure();
plot(abs(freqHamm));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of 2PAM signal (hamming window)');

figure();
plot(abs(freqHann));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of 2PAM signal (hanning window)');

% figure();
% plot(index, signalSinc);
% axis([0 1.5*length0*sps*Ts -2 2]);
% xlabel('time');
% ylabel('amplitude');
% title('basband transmission signal of 2PAM (Sinc Pulse)');

figure();
plot(abs(freqHann));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of 2PAM signal (sinc window)');

figure();
plot(abs(freqSincTrunk));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of truncated sinc (sinc window)');

%% ON/OFF Keying
% Transmitter setup
% define the basic rectangular pulse shape
pulse = ones(1,sps); % 1 pulse
pulseHammOnOff= hamming(sps);%creates hamming pulse
pulseHannOnOff=hann(sps);
pulseSincOnOff=sinc(time);

data = uint8(user_input_ascii);
%data = [1 0 1 1 0 0 1 1 1 0 0 0 1 1 0 1 1 0];

[number, length0] = size(data);
symbol = 2*data-1;
for i = 1: length0
    for j=1:sps
    signalOnOff((i-1)*sps + j) = symbol(i) * pulse(j);
    signalOnOffHamm((i-1)*sps + j) = symbol(i) * pulseHammOnOff(j);
    signalOnOffHann((i-1)*sps + j) = symbol(i) * pulseHannOnOff(j);
    signalOnOffSinc((i-1)*sps + j) = symbol(i) * pulseSincOnOff(j);
    end
end

index = 0: Ts: (length0*sps-1)*Ts;

delayedSignalOnOff=[delay, signalOnOff];
delayindex= 0:Ts:(length0*sps-1)*Ts+tau;

freq=fft(signalOnOff);
freq=fftshift(freq);

freqOnOffHamm=fft(signalOnOffHamm);
freqOnOffHamm=fftshift(freqOnOffHamm);

freqOnOffHann=fft(signalOnOffHann);
freqOnOffHann=fftshift(freqOnOffHann);

freqOnOffSinc=fft(signalOnOffSinc);
freqOnOffSinc=fftshift(freqOnOffSinc);
%
figure();
plot(index, signalOnOff);
axis([0 1.5*length0*sps*Ts -2 2]);
xlabel('time');
ylabel('amplitude');
title('basband transmission signal of On/off keying');
figure();
plot(abs(freq));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of On/Off Keying signal');
figure()
plot(delayindex, delayedSignalOnOff)
title('On/Off keying Signal with delay')
xlabel('time(sec)')
ylabel('symbol')
figure()
plot(abs(freqOnOffHamm));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of On/Off Keying signal (hamming window)');
figure()
plot(abs(freqOnOffHann));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of On/Off Keying signal (hanning window)');
figure()
plot(abs(freqOnOffSinc));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of On/Off Keying signal (sinc window)');

%% 4PAM

if length(user_input_ascii)/2~=1
    pad_zero=[0];
    padd_user_ascii=[user_input_ascii,pad_zero ];
else
    padd_user_ascii=user_input_ascii;
end




pulseHamm4= hamming(sps);%creates hamming pulse
pulseHann4=hann(sps);
pulseSinc4=sinc(time);

FourPamMod=[]; 
cnt=1;

%code that parses and 4PAM mods the input signal
for i=1:2:length(padd_user_ascii)-1

    twoBitTmp=[padd_user_ascii(i), padd_user_ascii(i+1)];
        if twoBitTmp== [0,0]
            FourPamMod(cnt)=-3;
        end
        if twoBitTmp== [0,1]
            FourPamMod(cnt)=-1;
        end
        if twoBitTmp== [1,0]
            FourPamMod(cnt)=3;
        end
        if twoBitTmp== [1,1]
            FourPamMod(cnt)=1;
        end
        cnt=cnt+1;
end

%Dans Transmitter
for i = 1: length(FourPamMod)
    for j=1:sps
    signal1((i-1)*sps + j) = FourPamMod(i) * pulse(j);
    signal4Hamm((i-1)*sps + j) = FourPamMod(i) * pulseHamm4(j);
    signal4Hann((i-1)*sps + j) = FourPamMod(i) * pulseHann4(j);
    %signal4Sinc((i-1)*sps + j) = FourPamMod(i) * pulseSinc4(j);
    end
end

signal4Sinc= sinc_mod(FourPamMod,sps);
plot(signal4Sinc)


freq4Hamm=fft(signal4Hamm);
freq4Hamm=fftshift(freq4Hamm);

freq4Hann=fft(signal4Hann);
freq4Hann=fftshift(freq4Hann);

freq4Sinc=fft(signal4Sinc);
freq4Sinc=fftshift(freq4Sinc);

delayedSignalFourPam=[delay, signal1];
delayedIndexFourPam= 0:Ts:(length(FourPamMod)*sps-1)*Ts+tau;

index1=0:Ts:(length(FourPamMod)*sps-1)*Ts;
figure()
plot(index1,signal1)
title('4 Pam Signal')
xlabel('time(sec)')
ylabel('symbol')
figure()
plot(delayedIndexFourPam, delayedSignalFourPam)
title('4 Pam Signal with delay')
xlabel('time(sec)')
ylabel('symbol')

figure()
plot(abs(freq4Hamm));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of 4PAM signal (Hamm window)')

figure()
plot(abs(freq4Hann));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of 4PAM signal (hann window)')

figure()
plot(abs(freq4Sinc));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of 4PAM signal (sinc window)')

%% 8 PAM 
cnt=1;

if mod(length(user_input_ascii),3)==2
    pad_zero=[0];
    padd_user_ascii2=[user_input_ascii , pad_zero];
elseif mod(length(user_input_ascii),3)==1
    pad_zero=[0, 0];
    padd_user_ascii2=[user_input_ascii , pad_zero];
else
    padd_user_ascii2=user_input_ascii;
end

for i=1:3:length(padd_user_ascii2)

    twoBitTmp=[padd_user_ascii2(i), padd_user_ascii2(i+1), padd_user_ascii2(i+2)];
        if twoBitTmp== [0,0,0]
            EightPamMod(cnt)=-7;
        end
        if twoBitTmp== [0,0,1]
            EightPamMod(cnt)=-5;
        end
        if twoBitTmp== [0,1,0]
            EightPamMod(cnt)=-1;
        end
        if twoBitTmp== [0,1,1]
            EightPamMod(cnt)=-3;
        end
        if twoBitTmp== [1,0,0]
            EightPamMod(cnt)=7;
        end
        if twoBitTmp== [1,0,1]
            EightPamMod(cnt)=5;
        end
        if twoBitTmp== [1,1,0]
            EightPamMod(cnt)=1;
        end
        if twoBitTmp== [1,1,1]
            EightPamMod(cnt)=3;
        end
        cnt=cnt+1;
end

pulseHamm8= hamming(sps);%creates hamming pulse
pulseHann8=hann(sps);
pulseSinc8=sinc(time);


for i = 1: length(EightPamMod)
    for j=1:sps
    signal8((i-1)*sps + j) = EightPamMod(i) * pulse(j);
    signal8Hamm((i-1)*sps + j) = EightPamMod(i) * pulseHamm8(j);
    signal8Hann((i-1)*sps + j) = EightPamMod(i) * pulseHann8(j);
    %signal8Sinc((i-1)*sps + j) = EightPamMod(i) * pulseSinc8(j);
    end
end


signal8Sinc= sinc_mod(EightPamMod,sps);
figure()
plot(signal8Sinc)


total_delay=0.0068;
signal8Sinc=[zeros(1,round(total_delay*fs)), signal8Sinc];
figure()
plot(signal8Sinc)

freq8Hamm=fft(signal8Hamm);
freq8Hamm=fftshift(freq8Hamm);

freq8Hann=fft(signal8Hann);
freq8Hann=fftshift(freq8Hann);

freq8Sinc=fft(signal8Sinc);
freq8Sinc=fftshift(freq8Sinc);

delayedSignalEightPam=[delay, signal8];
delayedIndexEightPam= 0:Ts:(length(EightPamMod)*sps-1)*Ts+tau;

index8=0:Ts:(length(EightPamMod)*sps-1)*Ts;
figure()
plot(index8,signal8)
title('8 Pam Signal')
xlabel('time(sec)')
ylabel('symbol')
figure()
plot(delayedIndexEightPam,delayedSignalEightPam)
title('8 Pam Signal with delay')
xlabel('time(sec)')
ylabel('symbol')

figure()
plot(abs(freq8Hamm));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of 8PAM signal (Hamm window)')

figure()
plot(abs(freq8Hann));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of 8PAM signal ((Hann window)')

figure()
plot(abs(freq8Sinc));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of 8PAM signal (sinc window)')

figure()
plot(signal8Sinc)




%% Part 2 Receiver

%%%%%%%%%%%%sample the recieved signal%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rectangular
sampled2PAM=samplerRect(delayedSignal, tau,fs,sps);
sampled4PAM=samplerRect(delayedSignalFourPam, tau, fs, sps);
sampled8PAM=samplerRect(delayedSignalEightPam, tau, fs, sps);
sampledOOK=samplerRect(delayedSignalOnOff, tau, fs, sps);
% %hamming
% sampled2PAMHamm=samplerRect(signalHamm, 0,fs,sps);
% sampled4PAMHamm=samplerRect(signal4Hamm, 0,fs,sps);
% sampled8PAMHamm=samplerRect(signal8Hamm, 0,fs,sps);
% %hanning
% sampled2PAMHann=samplerRect(signalHann, 0,fs,sps);
% sampled4PAMHann=samplerRect(signal4Hann, 0,fs,sps);
% sampled8PAMHann=samplerRect(signal8Hann, 0,fs,sps);
% %sinc
% sampled2PAMSinc=samplerRect(signalSinc, 0,fs,sps);
% sampled4PAMSinc=samplerRect(signal4Sinc,0,fs,sps);
sampled8PAMSinc=samplerSinc(signal8Sinc,total_delay,fs,sps);

%%%%%%%%%%%%%%%%% Quantize the signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rectangular
quantized2PAM=quantizer2PAM(sampled2PAM);
quantized4PAM=quantizer4PAM(sampled4PAM);
quantized8PAM=quantizer8PAM(sampled8PAM);
% %hamming
% quantized2PAMHamm=quantizer2PAM(sampled2PAMHamm);
% quantized4PAMHamm=quantizer4PAM(sampled4PAMHamm);
% quantized8PAMHamm=quantizer8PAM(sampled8PAMHamm);
% %hanning
% quantized2PAMHann=quantizer2PAM(sampled2PAMHann);
% quantized4PAMHann=quantizer4PAM(sampled4PAMHann);
% quantized8PAMHann=quantizer8PAM(sampled8PAMHann);
% %sinc
% quantized2PAMSinc=quantizer2PAM(sampled2PAMSinc);
% quantized4PAMSinc=quantizer4PAM(sampled4PAMSinc);
quantized8PAMSinc=quantizer8PAM(sampled8PAMSinc);


%%%%%%%%%%%%%%%%%%%% decode the signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
decoded2PAM=decoder2PAM(quantized2PAM);
decoded4PAM=decoder4PAM(quantized4PAM);
decoded8PAM=decoder8PAM(quantized8PAM);

decoded8PAMSinc=decoder8PAM(quantized8PAMSinc);

%%%%%%%%%%%%%%%%%%%%% Convert to ascii %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mess2PAM=char(bin2dec(reshape(decoded2PAM,7,[]).')).';
% disp(mess2PAM)
% mess4PAM=char(bin2dec(reshape(decoded4PAM,7,[]).')).';
% disp(mess4PAM)
% mess8PAM=char(bin2dec(reshape(decoded8PAM,7,[]).')).';
% disp(mess8PAM)


mess8PAMSinc=char(bin2dec(reshape(decoded8PAMSinc,7,[]).')).';
disp(mess8PAMSinc)
% mess2=char(bin2dec(reshape(sampled4PAM,7,[]).')).';
% disp(mess2)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sinc_mod=sinc_mod(signal,sps)
    timeSinc=linspace(-6,6,12*sps);
    slidingSinc=sinc(timeSinc);
    sum_sinc=slidingSinc*signal(1);
    for i=2:length(signal)
        current_sinc=slidingSinc*signal(i);
        current_sinc=[zeros(1,(i-1)*sps), current_sinc];
        sum_sinc=[sum_sinc, zeros(1,sps)];
        sum_sinc=current_sinc+sum_sinc;
    end
    sinc_mod=sum_sinc;
end

%%%%%%%%%%%%%%%%%%%% SAMPLERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sampler for rectangular pulse
function sampled_sig=samplerRect(signal, tau, fs, sps )
cnt=1;
    for i=round(tau*fs)+sps/2:sps:length(signal)
        sampled_sig(cnt)=signal(i);
        cnt=cnt+1;
    end
end

function sampled_sig=samplerSinc(signal, tau, fs, sps )
cnt=1;
    for i=round(tau*fs)+6*sps:sps:length(signal)-6*sps
        sampled_sig(cnt)=signal(i);
        cnt=cnt+1;
    end
end


%%%%%%%%%%%%%%%%%%%% QUANTIZERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quantizer
function quantized_sig=quantizer2PAM(signal)
    for i=1:length(signal)
        if signal(i)<0
            quantized_sig(i)=-1;
        else
            quantized_sig(i)=1;
        end
    end
end

function quantized_sig=quantizer4PAM(signal)
    for i=1:length(signal)
        if signal(i)<-2
            quantized_sig(i)=-3;
        end
        if signal(i)>2
            quantized_sig(i)=3;
        end
        if signal(i)>=0 && signal(i)<=2 
            quantized_sig(i)=1;
        end
        if signal(i)>=-2 && signal(i)<0 
            quantized_sig(i)=-1;            
        end
    end
end

function quantized_sig=quantizer8PAM(signal)
    for i=1:length(signal)
        if signal(i)<-6
            quantized_sig(i)=-7;
        end
        if signal(i)>6
            quantized_sig(i)=7;
        end
        if signal(i)>=-6 && signal(i)<-4 
            quantized_sig(i)=-5;
        end
        if signal(i)>=-4 && signal(i)<=-2 
            quantized_sig(i)=-3;            
        end
        if signal(i)>=-2 && signal(i)<0 
            quantized_sig(i)=-1;            
        end     
        if signal(i)>=0 && signal(i)<=2 
            quantized_sig(i)=1;            
        end  
         if signal(i)>2 && signal(i)<=4 
            quantized_sig(i)=3;
         end
         if signal(i)>4 && signal(i)<=6 
           quantized_sig(i)=5; 
         end        
    end
end

function quantized_sig=quantizerOOK(signal)
    for i=1:length(signal)
        if signal(i)<=0.5
            quantized_sig(i)=0;
        else
            quantized_sig(i)=1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% DECODERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function decode_sig=decoder2PAM(signal)
    for i=1:length(signal)
        if signal(i)==-1
            decode_sig(i)='0';
        else
            decode_sig(i)='1';
        end
    end
end

function decode_sig=decoder4PAM(signal)
    cnt=1;
    for i=1:length(signal)
        if signal(i)==-1
            decode_sig(cnt)=0;
            decode_sig(cnt+1)=1;
        elseif signal(i)==-3
            decode_sig(cnt)=0;
            decode_sig(cnt+1)=0;
        elseif signal(i)==1
            decode_sig(cnt)=1;
            decode_sig(cnt+1)=1;
        else
            decode_sig(cnt)=1;
            decode_sig(cnt+1)=0;            
        end
        cnt=cnt+2;
    end
    num_pad=0;
    len_decode=length(decode_sig);
    while mod(len_decode,7)~=0
        num_pad=num_pad+1;
        len_decode=len_decode-1;
    end
    if num_pad ~= 0
        decode_sig=decode_sig(:,1:length(decode_sig)-num_pad);
    end
    decode_sig = sprintf('%d', decode_sig);  % If a is a double.
end

function decode_sig=decoder8PAM(signal)
    cnt=1;
    for i=1:length(signal)
        if signal(i)==-7
            decode_sig(cnt)=0;
            decode_sig(cnt+1)=0;
            decode_sig(cnt+2)=0;            
        elseif signal(i)==-5
            decode_sig(cnt)=0;
            decode_sig(cnt+1)=0;
            decode_sig(cnt+2)=1;
        elseif signal(i)==-3
            decode_sig(cnt)=0;
            decode_sig(cnt+1)=1;
            decode_sig(cnt+2)=1;
        elseif signal(i)==-1
            decode_sig(cnt)=0;
            decode_sig(cnt+1)=1;
            decode_sig(cnt+2)=0;
        elseif signal(i)==1
            decode_sig(cnt)=1;
            decode_sig(cnt+1)=1;
            decode_sig(cnt+2)=0;
        elseif signal(i)==3
            decode_sig(cnt)=1;
            decode_sig(cnt+1)=1;
            decode_sig(cnt+2)=1;
        elseif signal(i)==5
            decode_sig(cnt)=1;
            decode_sig(cnt+1)=0;
            decode_sig(cnt+2)=1;
        elseif signal(i)==7
            decode_sig(cnt)=1;
            decode_sig(cnt+1)=0;
            decode_sig(cnt+2)=0;            
        end
        cnt=cnt+3;        
    end
        num_pad=0;
    len_decode=length(decode_sig);
    while mod(len_decode,7)~=0
        num_pad=num_pad+1;
        len_decode=len_decode-1;
    end
    if num_pad ~= 0
        decode_sig=decode_sig(:,1:length(decode_sig)-num_pad);
    end
    decode_sig = sprintf('%d', decode_sig);  % If a is a double.
end



