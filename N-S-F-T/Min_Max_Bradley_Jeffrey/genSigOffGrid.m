function [Sig,A,u] = genSigOffGrid(N,k,a_min,sigma_n)
% Generate 1-D signals with off-grid randomly distributed frequencies
% Input:
%   N: signal length 
%   k: sparsity level
%   a_min: the minimum sigal amplitude
%   sigma_n: noise variance length 
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
%% 位置移动，使不在格点
% u = mod(u+0.1*rand(K,1),N0);
u = mod(u+0.05,N);
%取模 mod(a,b) 与取余 rem（a,b)操作的区别
%当 a 和 b 都是正数的时候, 二者结果一样: a 除以 b 后的余数
%当任何一个位置出现负数的时候, 先按正数算出结果的绝对值. 然后对于结果 mod 取和b 一样的符号,rem 取和 a 一样的符号

% A/a和A./a ：矩阵A中元素都除以a，输出矩阵
% A.\a：a除以矩阵A中各元素，输出矩阵
% A/B：相当于A乘以B的逆
% A./B：矩阵点右除：要求两矩阵维数相等，即MxN维矩阵除以MxN维矩阵。矩阵对应位置元素相除输出，A矩阵对应元素除以B矩阵对应元素
% A\B：相当于A的逆乘以B
% A.\B:矩阵点左除：要求两矩阵维数相等，即MxN维矩阵除以MxN维矩阵。矩阵对应位置元素相除输出，B矩阵对应元素除以A矩阵对应元素
%%
%构造含离点频率的时域信号
A = a_min*exp(1i*2*pi*rand(k,1));%随机生成信号的振幅
fu = u/N; % 真实频率值
% w = 2*pi*fu;% 角频率
Sig=zeros(1,N);
% UP_Sig=zeros(1,Len);
for ii=1:k
    %每个频率对应的时域信号叠加
     sig=A(ii)*exp(1i*2*pi*fu(ii)*(0:N-1));
%      sig=A(ii)*sin(w(ii)*(0:N-1));
     Sig = Sig+sig;
%      sigK=A(ii)*exp(1i*2*pi*fu(ii)*(0:Len-1));
%      sigK=A(ii)*sin(w(ii)*(0:Len-1));
%      UP_Sig = UP_Sig+sigK;
end
noise = sigma_n*sqrt(2)/2*(randn(1,N)+1i*randn(1,N));%噪声
Sig = Sig+noise;%含躁信号
% noise2 = sigma_n*sqrt(2)/2*(randn(1,Len)+1i*randn(1,Len));%噪声
% UP_Sig = UP_Sig+noise2;%含躁信号

end