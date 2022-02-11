%% SDR Lab 2
%% PART 2.2

% Signal Characteristics
fc=1e3;
fm=100;
fs=80e3;
B=[0.1, 0.3, 0.5, 1, 2, 5, 10, 30];

%% Modulated Signal

t=[0:1/fs:.01];
for b=1:8
    s=0;
    S=0;
    f=[-5e3:5e3];
    for n=-100:100
        s=s+(besselj(n,B(b))*cos((2*pi*(fc+n*fm).*t))); % Time Domain Singal-Tone FM Signal 
    end
    figure(b)
    plot(t,s)
    
    % Frequency Analysis
    figure(b+8)
    for n=-100:100
        impulse1= f== fc + n*fm;
        impulse2= f== -fc-n*fm;
        S=S+(besselj(n,B(b))*(impulse1+impulse2)); % Frequency Domain Singal-Tone FM Signal 
    end
    f=f+fc;
    plot(f,abs(.5*S))
end

%% Estimate of Transmission BW
Transmission_BW = 2*fc*(1+(1/B(4))); % ~176 kHz => Actual is 180 kHz