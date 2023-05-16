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
% �������Ƶ�ʵ�K��λ�ã��ɲ��ڸ�㣩
u = randi(N,[k,1])-1;
uInd = unique(u);
while size(uInd,1)<k
    u = randi(N,[k,1])-1;
    uInd = [uInd;u];
    uInd = unique(uInd);
end
rdI = randperm(size(uInd,1));%��һ�����������ң���ű�����������
uInd = uInd(rdI(1:k));
u  = sort(uInd);%������ɵ�K�����λ��

% u_off = mod(u+0.5*rand(k,1),N);%λ���ƶ���ʹ���ڸ��
u_off = mod( u+0.01, N);
%ȡģ mod(a,b) ��ȡ�� rem��a,b)����������
%�� a �� b ����������ʱ��, ���߽��һ��: a ���� b �������
%���κ�һ��λ�ó��ָ�����ʱ��, �Ȱ������������ľ���ֵ. Ȼ����ڽ�� mod ȡ��b һ���ķ���,rem ȡ�� a һ���ķ���
 fu = u_off/N;
% fv = v/N1;
% A/a��A./a ������A��Ԫ�ض�����a���������
% A.\a��a���Ծ���A�и�Ԫ�أ��������
% A/B���൱��A����B����
% A./B��������ҳ���Ҫ��������ά����ȣ���MxNά�������MxNά���󡣾����Ӧλ��Ԫ����������A�����ӦԪ�س���B�����ӦԪ��
% A\B���൱��A�������B
% A.\B:����������Ҫ��������ά����ȣ���MxNά�������MxNά���󡣾����Ӧλ��Ԫ����������B�����ӦԪ�س���A�����ӦԪ��
%%  ���캬���Ƶ�ʵ�Ƶ���ź�
% A = a_min*exp(1i*2*pi*rand(k,1));%��������źŵ����
A = N/256*exp(1i*2*pi*rand(k,1));%��������źŵ����
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
%% ���캬���Ƶ�ʵ�ʱ���ź�
Sig = indft1d(Sig_freq, N, f);
%  Sig = awgn(Sig,SNR_dB,'measured');%�������ΪSNR_dB�ĸ�˹������

% NOISE=randn(N);
% NOISE=NOISE-mean(NOISE);
% signal_power = 1/length(x)*sum(x.*x);
% noise_variance = signal_power / ( 10^(SNR/10) );
% NOISE=sqrt(noise_variance)/std(NOISE)*NOISE;
% NOISE=normrnd(0,sigma_n,N);
% Sig = Sig+NOISE;%�����ź�
end