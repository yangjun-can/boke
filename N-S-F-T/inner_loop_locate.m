%% 定位内循环
function [x_samp,J ]=inner_loop_locate(origx, n, cur_filter_time,cur_filter_sizet, ...
    num,cur_B, ai,b)
%% 参数
% origx 原始信号
% n 信号长度
% cur_filter_time 平坦窗滤波器（时域）
% cur_filter_sizet 平坦窗滤波器窗长
% cur_filter_freq 平坦窗滤波器（频域）
% num 大桶个数
% cur_B 分桶数
% a ai的模逆
% ai 重排时域下标伸缩参数（相当于课本上的sigma）
% b 重排时域下标平移参数 （相当于课本上的-a*sigma）
% x_samp ：桶值
% J ：num个大桶坐标集
%%  要求信号长度是分桶数的整数倍，否则输出Warning
if mod(n,cur_B)       
   fprintf('Warning: n is not divisible by cur_B, which algorithm expects.\n'); 
end
%% 分桶操作
x_sampt =zeros(1,n,'single');  
index=b;
for i=0:(cur_filter_sizet-1)   
    x_sampt(mod(i,cur_B)+1) = origx(index+1) * cur_filter_time(i+1) + x_sampt(mod(i,cur_B)+1); 
    index =rem((index+ai),n);	
    % origx(index+1)：随机重排，重排坐标index = b+ai*i,故 ai=sigma,b=-a*sigma
    % origx(index+1) * cur_filter_time(i+1)：滤波（时域相乘，频域卷积）
    % x_sampt(mod(i,cur_B)+1) =....： 频率降采样（对应的时域值）
end
x_samp=fft(x_sampt,cur_B);                    % 降采样频域上的采样值

%% 寻找大桶位置集
samples = zeros(1,cur_B);                     % 降采样后的功率谱
for i=1:cur_B
    samples(i) = (abs(x_samp(i)))^2;							
end
% figure;plot(samples);
J=find_largest_indices( num, samples ,cur_B);  % num个大桶的坐标
end