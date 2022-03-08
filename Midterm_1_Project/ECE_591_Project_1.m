% ECE 591 
% Software Defined Radio
% Midterm Project

clc;clear;close all

%% Part 1
user_input=input('Enter the string you would like to send: ', 's');
user_input_ascii= int2bit(double(user_input),7);
user_input_ascii= reshape(user_input_ascii,[1,length(user_input)*7]);



fs = 48000; % sampling frequency
Ts = 1/fs; % sampling duration
symbolrate = 1000; % transmitted pulses/second should be an integer divisor of fs
sps = fs/symbolrate; % number of samples in one symbol
tau=1/sps;
delay=zeros(1,round(tau*fs));

%% 2 PAM
% Transmitter setup
% define the basic rectangular pulse shape

pulse = ones(1,sps); % 1 pulse
%delayedPulse=[delay, pulse];
data = user_input_ascii;
%data = [1 0 1 1 0 0 1 1 1 0 0 0 1 1 0 1 1 0];

[number, length0] = size(data);
symbol = 2*data-1;
for i = 1: length0
    for j=1:sps
        signal((i-1)*sps + j) = symbol(i) * pulse(j);
    end
end

%create delayed signal
delayedSignal=[delay, signal];
delayindex= 0:Ts:(length0*sps-1)*Ts+tau;

index = 0: Ts: (length0*sps-1)*Ts;
freq=fft(signal);
freq=fftshift(freq);
%
figure(1);
plot(index, signal);
axis([0 1.5*length0*sps*Ts -2 2]);
xlabel('time');
ylabel('amplitude');
title('basband transmission signal of 2PAM');
figure(2);
plot(abs(freq));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of 2PAM signal');
figure(3)
plot(delayindex, delayedSignal)
title('2 Pam Signal with delay')
xlabel('time(sec)')
ylabel('symbol')

%% ON/OFF Keying
% Transmitter setup
% define the basic rectangular pulse shape
pulse = ones(1,sps); % 1 pulse
data = uint8(user_input_ascii);
%data = [1 0 1 1 0 0 1 1 1 0 0 0 1 1 0 1 1 0];

[number, length0] = size(data);
symbol = 2*data-1;
for i = 1: length0
    for j=1:sps
    signalOnOff((i-1)*sps + j) = symbol(i) * pulse(j);
    end
end
index = 0: Ts: (length0*sps-1)*Ts;

delayedSignalOnOff=[delay, signalOnOff];
delayindex= 0:Ts:(length0*sps-1)*Ts+tau;

freq=fft(signalOnOff);
freq=fftshift(freq);
%
figure(4);
plot(index, signalOnOff);
axis([0 1.5*length0*sps*Ts -2 2]);
xlabel('time');
ylabel('amplitude');
title('basband transmission signal of On/off keying');
figure(5);
plot(abs(freq));
xlabel('frequency');
ylabel('magnitude');
title('spectrum of On/Off Keying signal');
figure(6)
plot(delayindex, delayedSignalOnOff)
title('On/Off keying Signal with delay')
xlabel('time(sec)')
ylabel('symbol')

%% 4PAM
if length(user_input_ascii)/2~=1
    pad_zero=[0];
    padd_user_ascii=[pad_zero, user_input_ascii ];
end

FourPamMod=[]; 
cnt=1;

%code that parses and 4PAM mods the input signal
for i=1:2:length(padd_user_ascii)

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
    end
end

delayedSignalFourPam=[delay, signal1];
delayedIndexFourPam= 0:Ts:(length(FourPamMod)*sps-1)*Ts+tau;

index1=0:Ts:(length(FourPamMod)*sps-1)*Ts;
figure(7)
plot(index1,signal1)
title('4 Pam Signal')
xlabel('time(sec)')
ylabel('symbol')
figure(8)
plot(delayedIndexFourPam, delayedSignalFourPam)
title('4 Pam Signal with delay')
xlabel('time(sec)')
ylabel('symbol')

%% 8 PAM 
cnt=1;

if mod(length(user_input_ascii),3)==2
    pad_zero=[0];
    padd_user_ascii2=[pad_zero, user_input_ascii ];
elseif mod(length(user_input_ascii),3)==1
    pad_zero=[0, 0];
    padd_user_ascii2=[pad_zero, user_input_ascii ];
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

for i = 1: length(EightPamMod)
    for j=1:sps
    signal8((i-1)*sps + j) = EightPamMod(i) * pulse(j);
    end
end

delayedSignalEightPam=[delay, signal8];
delayedIndexEightPam= 0:Ts:(length(EightPamMod)*sps-1)*Ts+tau;

index8=0:Ts:(length(EightPamMod)*sps-1)*Ts;
figure(9)
plot(index8,signal8)
title('8 Pam Signal')
xlabel('time(sec)')
ylabel('symbol')
figure(10)
plot(delayedIndexEightPam,delayedSignalEightPam)
title('8 Pam Signal with delay')
xlabel('time(sec)')
ylabel('symbol')

