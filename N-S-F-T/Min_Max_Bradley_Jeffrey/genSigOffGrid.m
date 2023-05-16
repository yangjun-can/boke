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
%% λ���ƶ���ʹ���ڸ��
% u = mod(u+0.1*rand(K,1),N0);
u = mod(u+0.05,N);
%ȡģ mod(a,b) ��ȡ�� rem��a,b)����������
%�� a �� b ����������ʱ��, ���߽��һ��: a ���� b �������
%���κ�һ��λ�ó��ָ�����ʱ��, �Ȱ������������ľ���ֵ. Ȼ����ڽ�� mod ȡ��b һ���ķ���,rem ȡ�� a һ���ķ���

% A/a��A./a ������A��Ԫ�ض�����a���������
% A.\a��a���Ծ���A�и�Ԫ�أ��������
% A/B���൱��A����B����
% A./B��������ҳ���Ҫ��������ά����ȣ���MxNά�������MxNά���󡣾����Ӧλ��Ԫ����������A�����ӦԪ�س���B�����ӦԪ��
% A\B���൱��A�������B
% A.\B:����������Ҫ��������ά����ȣ���MxNά�������MxNά���󡣾����Ӧλ��Ԫ����������B�����ӦԪ�س���A�����ӦԪ��
%%
%���캬���Ƶ�ʵ�ʱ���ź�
A = a_min*exp(1i*2*pi*rand(k,1));%��������źŵ����
fu = u/N; % ��ʵƵ��ֵ
% w = 2*pi*fu;% ��Ƶ��
Sig=zeros(1,N);
% UP_Sig=zeros(1,Len);
for ii=1:k
    %ÿ��Ƶ�ʶ�Ӧ��ʱ���źŵ���
     sig=A(ii)*exp(1i*2*pi*fu(ii)*(0:N-1));
%      sig=A(ii)*sin(w(ii)*(0:N-1));
     Sig = Sig+sig;
%      sigK=A(ii)*exp(1i*2*pi*fu(ii)*(0:Len-1));
%      sigK=A(ii)*sin(w(ii)*(0:Len-1));
%      UP_Sig = UP_Sig+sigK;
end
noise = sigma_n*sqrt(2)/2*(randn(1,N)+1i*randn(1,N));%����
Sig = Sig+noise;%�����ź�
% noise2 = sigma_n*sqrt(2)/2*(randn(1,Len)+1i*randn(1,Len));%����
% UP_Sig = UP_Sig+noise2;%�����ź�

end