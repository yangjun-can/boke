%% 定位内循环 (分桶)
function [x_samp,J]=nusfft_inner_loop_locate(origx, n, flat_filter_time, ...
    num,cur_B,ai,b,pe_major,inv_pe_major, sa_major,f)

% function [x_samp,J]=nusfft_inner_loop_locate(origx, n, flat_filter_time, ...
%     num,cur_B,sa_major,permute,f)

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
% x2=zeros(n,1);
% index=b;
%% 随机重排
x1=inv_pe_major*origx.';
% for i=0:(n-1)
%     x2(i+1)=x1(index+1);
%     index =rem((index+ai),n);
% end
% x2=x1(rem(b+(0:n-1)*ai,n)+1);
% x3=permute*origx.';
x3=pe_major*x1(rem(b+(0:n-1)*ai,n)+1);
%  figure;
%  subplot(2,1,1);plot(abs(x3));title('(a)');ylabel('amplitude');xlabel('Sampling point in time domion');
%  subplot(2,1,2);plot(f,abs(ndftld(x3.', n, f)));title('(b)');ylabel('amplitude');xlabel('Sampling point in frequency domion');
 %% 加平坦窗滤波
x4=x3.* flat_filter_time;
%% 频域下采样
% x5=sqrt(n/cur_B)*sa_major*x4.';
x5=n/cur_B*sa_major*x4;
% x5=x5.';
%% NUDFT
f2=Nonuniform_sampling_point(0:cur_B-1,cur_B);
f2=f2(:);
x_samp =chebfun.nufft(x5,f2/cur_B,2);
% x_samp = ndftld(x5,cur_B,f2).';
% x_samp = fft(x5);
%  figure;
%  subplot(2,1,1);plot(abs(x5));title('(a)');ylabel('amplitude');xlabel('Sampling point in time domion');
%  subplot(2,1,2);stem(f2,abs(x_samp));title('(b)');ylabel('amplitude');xlabel('Sampling point in frequency domion');
%% 寻找大桶位置集
% samples = zeros(1,cur_B);                     % 降采样后的功率谱
% for i=1:cur_B
%     samples(i) = (abs(x_samp(i)))^2;
% end
samples = (abs(x_samp)).^2;
J=find_largest_indices( num, samples ,cur_B);  % num个大桶的坐标
end