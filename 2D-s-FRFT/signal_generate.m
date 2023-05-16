function [ data_r, t_fast, ta, Nr, Na] = signal_generate( fc, B, fs, pulse_duration, PRF, pulse_number, amplitude, R0, v, aa, flag)
%   LFM pulse signal generation function.
%   Input:
%   fc: Carrier frequency.  B: Bandwidth of signal. fs: Sampling rate.
%   pulse_duration: Pulse duration.  PRF: Pulse repetition frequency.
%   pulse_number: Pulse number to generate. 
%   amplitude, R0, v, aa: Target motion parameter vectors. 
%   flag: A flag to check the first n targets.
%   Output:
%   data_r: Signal after pulse compression.
%   t_fast: Fast time axis.  ta: Slow time axis.
%   Nr: Point number in a pulse.
%   Na: Pulse number

if (nargin < 8)
    error('Not enough input parameters!');    
end
if (nargin == 8)
    if length(amplitude)~=length(R0)
       error('Check the length of target motion vectors!'); 
    end
    target_num=length(R0);
    v=zeros(1,target_num);
    aa=zeros(1,target_num);
    flag=target_num;
    warning('We set v=0 and aa=0.');     
end
if (nargin == 9)
    if length(amplitude)~=length(R0)||length(R0)~=length(v)||length(amplitude)~=length(v)
       error('Check the length of target motion vectors!'); 
    end
    target_num=length(R0);
    aa=zeros(1,target_num);
    flag=target_num;
    warning('We set v=0.');     
end
if (nargin == 10)
    if length(amplitude)~=length(R0)||length(R0)~=length(v)||length(amplitude)~=length(v)...
       ||length(amplitude)~=length(aa)||length(R0)~=length(aa)||length(v)~=length(aa)
       error('Check the length of target motion vectors!'); 
    end
    flag=target_num;
    warning('We will get the result of all targets.');     
end
if (nargin == 11)
    if length(amplitude)~=length(R0)||length(R0)~=length(v)||length(amplitude)~=length(v)...
       ||length(amplitude)~=length(aa)||length(R0)~=length(aa)||length(v)~=length(aa)
       error('Check the length of target motion vectors!'); 
    end  
end

%% Fast time parameters.
c=3e8;
lamda = c/fc;
if fs<B
   error('fs<B. Signal is undersampling!');
end
ts=1/fs;     % sampling interval
pulse_num=floor(pulse_duration/ts);    % Samples in a pulse.
kr = B/pulse_duration;   %chirp rate

%% Slow time parameters.
pri = 1/PRF;
Na = pulse_number;    % pulse number
ta = (-Na/2:Na/2-1)*pri;   % slow time axis

%% Target trajectory

target_num=length(R0);
if flag>target_num
   flag=target_num;
end
if flag<1
   flag=1; 
end
for i=1:flag  % We can choose the number of targets to be used.
    R(i,:)=R0(i)+v(i)*ta+0.5*aa(i)*ta.^2;   
end

%% Signal generate
delay_max=2*max(max(R))/c;    % max delay
delay_min=2*min(min(R))/c;    % min delay
t_fast = (delay_min-pulse_num*ts):ts:(delay_max+pulse_num*ts);   % fast time axis

t_start=t_fast(1);  % The start point of the fast time axis.
r_axis_real=t_fast*c/2;     % Real range axis.
Nr = length(t_fast);    % Length of fast time axis.

receiver_sig=zeros(Na,Nr);
for tt=1:flag
    for i = 1:Na
        delay = 2*R(tt,i)/c;
        receiver_sig(i,:) = receiver_sig(i,:)+amplitude(tt)*rectpuls(t_fast-pulse_duration/2-delay,pulse_duration).*exp(1i*2*pi*fc*(-delay)).*exp(1i*pi*kr*(t_fast-pulse_duration/2-delay).^2);   
    end
end
st=rectpuls(t_fast-t_start-pulse_duration/2,pulse_duration).*exp(1i*pi*kr*(t_fast-t_start-pulse_duration/2).^2);   % Reference signal

fft_st=fft(st);
Nr_fft=fft(receiver_sig,[],2);   
for na=1:Na
    data_r(na,:)=ifft((Nr_fft(na,:)).*conj(fft_st));  % Pulse compression.
end
% figure  
% imagesc(r_axis_real,ta,abs(data_r)); 
% xlabel('Range (m)');ylabel('Slow time (s)');
% title('Range Walk');
end

