%% 利用切比雪夫窗函数生成的平坦窗函数
function [out_time,out_sizet,out_freq]=make_multiple_t(filtert, w, n, b)
%% 参数
% filtert 切比雪夫滤波器（频域）
% w 滤波器的窗长
% n 信号长度
% b 滤波器个数
% out_time  平坦窗滤波器（时域）
% out_sizet 平坦窗滤波器窗长
% out_freq  平坦窗滤波器（频域）
%%  滤波器的窗长和数量均不超过n
if b >= n || w >= n
    return;
end
%% 
g=zeros(1,n);
h=zeros(1,n);
w1=fix(w/2);
g(1:w-w1)=filtert((w1+1):w);
g(n-w1+1:n)=filtert(1:w1);
% figure;plot(abs(g));title('win');
g=fft(g,n);

% s = 0;
% for i=1:b
%     s = g(i)+s ;									% 滤波器能量
% end
s=sum(g(1:b));
offset = fix(b/2);
for i=0:(n-1)
    h(mod((i+offset),n)+1) = s;
    s = s + (g(mod((i+ b),n)+1) - g(i+1));
end

h_max=max(abs(h));
% for i=1:n
%     h(i) =h(i)/ h_max;
% end
h = h/ h_max;
step=exp(-2j*pi*w1/n);

% offsetc = 1;
% for i=1:n
%     h(i)= h(i)*offsetc;
%     offsetc= offsetc *step;							% 在频域上对每个幅值进行移位
% end
h=h.*step.^(0:n-1);
g=ifft(h,n);
filtert(1:w)=g(1:w);

out_time = filtert;								    % 时域信号
out_sizet = w;										% 滤波器频宽
out_freq = h;										% 频域信号
end