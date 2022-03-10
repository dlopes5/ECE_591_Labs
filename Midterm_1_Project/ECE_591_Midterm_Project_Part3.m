% ECE 591 Software Defined Radio
% Midterm Project 1
% Part 2
% Daniel Lopes, Ryan Ferreira, Kevin Ventura


%% Part 3: Decode the Given Data

clc;clear;close all;
load trans_signal.mat;
who
Ts=1/fs;
sps=fs/symbolrate;
delaysamples=round(fs*tau);
index=0:Ts:(length(signal)-1)*Ts;
plot(index,signal)
real_delay=tau+delta;

sampled_sig=samplerSinc(signal, real_delay, fs,sps);
quantized_sig=quantizer8PAM(sampled_sig, gain);
decoded_sig=decoder8PAM(quantized_sig);

mess8PAMSinc=char(bin2dec(reshape(decoded_sig,7,[]).')).';
disp(mess8PAMSinc)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%S SAMPLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sampled_sig=samplerSinc(signal, tau, fs, sps )
cnt=1;
    for i=round(tau*fs)+6*sps:sps:length(signal)-6*sps
        sampled_sig(cnt)=signal(i);
        cnt=cnt+1;
    end
end

%%%%%%%%%%%%%%%%%%%% QUANTIZERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function quantized_sig=quantizer8PAM(signal, g)
    signal=signal/g;
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

%%%%%%%%%%%%%%%%%%%%%%%% DECODERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

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

