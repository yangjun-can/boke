function [Sig,Sig_freq,f,A,u] = genSigOffGrid(N,k,SNR_dB)
% Generate 1-D signals with off-grid randomly distributed frequencies
% Input:
%   N0: signal length
%   K: sparsity level
%   a_min: the minimum sigal amplitude
%   sigma_n: noise variance
% Output:
%   Sig: signal (time)
%   A: amplitude of each frequency (K-length vector)
%   u: frequency locations
%%
% 随机生成频率的K个位置（可不在格点）
u = randi(N,[k,1])-1;
uInd = unique(u);
while size(uInd,1)<k
    u = randi(N,[k,1])-1;
    uInd = [uInd;u];
    uInd = unique(uInd);
end
rdI = randperm(size(uInd,1));%将一列序号随机打乱，序号必须是整数。
uInd = uInd(rdI(1:k));
u  = sort(uInd);%随机生成的K个格点位置

% u_off = mod(u+0.5*rand(k,1),N);%位置移动，使不在格点
u_off = mod( u+0.01, N);
%取模 mod(a,b) 与取余 rem（a,b)操作的区别
%当 a 和 b 都是正数的时候, 二者结果一样: a 除以 b 后的余数
%当任何一个位置出现负数的时候, 先按正数算出结果的绝对值. 然后对于结果 mod 取和b 一样的符号,rem 取和 a 一样的符号
 fu = u_off/N;
% fv = v/N1;
% A/a和A./a ：矩阵A中元素都除以a，输出矩阵
% A.\a：a除以矩阵A中各元素，输出矩阵
% A/B：相当于A乘以B的逆
% A./B：矩阵点右除：要求两矩阵维数相等，即MxN维矩阵除以MxN维矩阵。矩阵对应位置元素相除输出，A矩阵对应元素除以B矩阵对应元素
% A\B：相当于A的逆乘以B
% A.\B:矩阵点左除：要求两矩阵维数相等，即MxN维矩阵除以MxN维矩阵。矩阵对应位置元素相除输出，B矩阵对应元素除以A矩阵对应元素
%%  构造含离点频率的频域信号
% A = a_min*exp(1i*2*pi*rand(k,1));%随机生成信号的振幅
A = N/256*exp(1i*2*pi*rand(k,1));%随机生成信号的振幅
index=1;
Sig_freq=zeros(1,N);
f=zeros(N,1);
for i=0:N-1
    if( ismember(i,u) )
        Sig_freq(i+1)=A(index);
%         f(i+1)=u_off(index); 
        index=index+1;
    end
end
for i=1:N
    f(i)=Nonuniform_sampling_point(i-1,N);
end
%% 构造含离点频率的时域信号
Sig = indft1d(Sig_freq, N, f);
%  Sig = awgn(Sig,SNR_dB,'measured');%加信噪比为SNR_dB的高斯白噪声

% NOISE=randn(N);
% NOISE=NOISE-mean(NOISE);
% signal_power = 1/length(x)*sum(x.*x);
% noise_variance = signal_power / ( 10^(SNR/10) );
% NOISE=sqrt(noise_variance)/std(NOISE)*NOISE;
% NOISE=normrnd(0,sigma_n,N);
% Sig = Sig+NOISE;%含躁信号
end