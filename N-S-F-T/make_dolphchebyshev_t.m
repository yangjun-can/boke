%% 构造切比雪夫滤波器
function [ filter ,w] = make_dolphchebyshev_t(lobefrac,tolerance)
%% 参数
% lobefrac 通带截频
% tolerance 泄露容限
%% 窗长
w =fix((1 / pi) * (1/lobefrac) * acosh(1/tolerance)); %  窗长  fix():向0取整
if ~mod(w,2)==1                                       % ~：非
    w=w-1;						                      % 保证w是奇数，滤波器是对称的
end
%% 滤波器
% filter = zeros(1,w);
t0 = cosh(acosh(1/tolerance) / (w-1));

% for i=0:(w-1)
%     filter(i+1) = Cheb(w-1, t0*cos(pi*i/w))*tolerance; % 切比雪夫滤波器的时域
% end
 filter = Cheb(w-1, t0*cos(pi*(0:w-1)/w))*tolerance;

%figure;plot(abs(filter));title('win');
temp=fft(filter, w);			                     % 切比雪夫滤波器的频域
%figure;plot(real(temp));title('FREwin');
filter=fftshift(temp);                               % 将零频点移到频谱的中间（将fft处理之后的pi-2pi部分搬移至-pi-0，从而使零频分量居于频谱的中心位置。）
%figure;plot(real(temp));title('FREwin_half');

% for i=1:w
%     filter(i) = real(filter(i));                	 % 返回的是实部
% end
filter = real(filter);
end