function [out_Large,BB_loc,b_loc]=nusfft_iter(x,k,f)
% 迭代的 NUSFT算法
%% 参数初始化
% profile on;
x=single(x);
f=single(f);
n=length(x);           % 信号长度
%k                     % 信号稀疏度
repetitions = 1;	   % 迭代次数
Bcst_loc=16;            % loc分桶常量
tolerance_loc = 1e-1;     % loc噪声（泄露）容限 test中n=64,取1e-4； test中n>64,取1e-6；
BB_loc = fix(Bcst_loc*sqrt(n*k/(log2(n))));
B_loc = floor_to_pow2(BB_loc);                 % loc所用B
B_thresh = 2*k;                                % 大桶的个数
lobefrac_loc = 0.5 / BB_loc;				   % loc通带区域[-lobefrac_loc,lobefrac_loc]
b_loc = fix(1.2*1.1*(n/BB_loc));
% 预计算重排矩阵
sequence1=gpuArray(single(exp(-2i*pi*(0:n-1).*f(1:n)/n))); 
pe_major=n*ifft(sequence1);
inv_pe_major=inv(pe_major);
% 预计算下采样矩阵
sequence2=gpuArray(single(exp(-2i*pi*(0:B_loc-1).*(f((1:B_loc)*n/B_loc)-(1:B_loc)'*(n/B_loc-1))/B_loc)));  % 预计算下采样矩阵
sa_major1=B_loc*ifft(sequence2); 
sa_major=repmat(sa_major1,[1,n/B_loc]);
sa_major=sa_major/pe_major;
%% NUSFFT估计频率
%% 大桶数不超过桶数
if B_thresh > B_loc 
    return;
end
%% 构造平坦窗滤波器
[filtert,w_loc ]=make_dolphchebyshev_t(lobefrac_loc, tolerance_loc);                 % 构造 filtert 切比雪夫滤波器，w_loc：窗长
[filter_time, filter_sizet, filter_freq]=make_multiple_t(filtert, w_loc, n, b_loc);    % 构造平坦窗滤波器，filter_time：滤波器时域 filter_sizet：滤波器窗长 filter_freq：滤波器频域
filter_time=single(filter_time);
filter_freq=single(filter_freq);
% aaaaa=fft(filtert);
% figure;plot(abs(aaaaa));title('dolphchebyshev filter time domain - window function');
% figure;
% subplot(2,2,1);plot(abs(filter_time));title('(a)');ylabel('amplitude');
% subplot(2,2,3);plot(abs(filter_freq));title('(c)');ylabel('amplitude');
% subplot(2,2,2);semilogy(abs(filter_time));title('(b)');ylabel('amplitude');
% subplot(2,2,4);semilogy(abs(filter_freq));title('(d)');ylabel('amplitude');
%% 迭代
% for ii=1:repetitions
    % 外循环，每次得到heavy值的位置和大小
        [I,out_I,out]= iter(x, n,filter_time,filter_sizet, filter_freq, ...
            B_thresh, B_loc,pe_major,inv_pe_major,sa_major,f);    
    % I     大值的位置（升序）
    % out_I 对应的估值
    % out   估值关于位置的函数
% end
%% 用估计结果构造信号的频域
num_candidates= length(out_I);      % 所估计的大值个数
out_Large = zeros(1,n);
candidates=zeros(2,num_candidates);
% 将算法估计的结果保存在二维矩阵，第二行是位置，第一行是相应的幅值
candidates(1,:)=abs(out_I);  % 大值的幅度序列
candidates(2,:)=I;           % 大值的位置序列
temp=sortrows((candidates).').';          % 大值按幅值升序排序
% 取K个最大的幅值，组成信号频域 out_Large
key=temp(2,num_candidates - k+ (1:k));
out_Large(key+1) = out(key+1)*n;
% profile viewer;
end