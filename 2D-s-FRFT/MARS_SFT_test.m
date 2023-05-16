%% Test of Robust MARS-SFT
% clear;
close all;
N0=2048;
N1=2048;
N = N0*N1;
L = lcm(N0,N1); % line length

K = 10;  % signal sparsity
gama = 74;  % ��λ��������
T = 1;   % ��������

% n_d = 2;  % ��ѭ������
% n_s = 3;  % ��ͶƱ����
%% Generate signal
sigma_n = 1; % ������ʱ���׼�� set this as 0 for noiseless signals, otherwise set as 1
sigma_f = sqrt(N)*sigma_n;  % Ƶ�������ı�׼��
if sigma_n==1
%  A = sqrt(sigma_n^2*SNR);     % chirp�ź�ʱ��ĵ�Ƶ��ֵ
  SNR_dB = 60;  % ��СƵ�� �����
  SNR = 10^(SNR_dB/10);
  a_min = sqrt(sigma_f^2*SNR); % ��С����
%     Pd=0.99;  % ��Ͱ������
%     Pfa=Pd^(L*SNR+1); % ��Ͱ�龯����
     Pfa=1e-2;
     epsilon = -2*L*log(Pfa)*sigma_n^2; %Threshold for frequency detectoin in a slice; set as 1e-10 for noiseless case
     gamma = 1*1e4;   % Threshold for 1-sparse detection; set as 1e-10 for noiseless case
else
    a_min = 10^(SNR_dB/10);
    epsilon = 1e-10;% %Threshold for frequency detectoin in a slice; set as 1e-10 for noiseless case
    gamma = 1e-10; % Threshold for 1-sparse detection; set as 1e-10 for noiseless case
end

% Generate on-grid frequencies
 [Sig, A, u, v] = genSigOnGrid(N0,N1,K,a_min,sigma_n);

% Generate off-grid frequencies
%[Sig, A, u, v] = genSigOffGrid(N0,N1,K,a_min,sigma_n);

% Generate on-grid frequency clusters
%[Sig, A, u, v] = genSigOnGridClusters(N0,N1,K,2,a_min,sigma_n); K = K*25;

%% the groud truth
Omega_gt = [u,v];  % ground truth frequency locations
X_true=zeros(N0,N1);
X_true(sub2ind([N0,N1],u+1,v+1)) = A; % ��ʵƵ��
% figure;mesh(0:N0-1,0:N1-1,abs(X_true));
%%  fft
% tic
% X2=fft2(Sig);
% disp([' fft2��ʱ��Ϊ ',  num2str(toc)]);
% figure;mesh(0:N0-1,0:N1-1,abs(X2));title("FFT2���")
% error2=norm(X2-X_true)
%% Generate window; genearate rect window for on-grid cases
%% Chebyshev window
% att = 70; % PSR in dB
% rho_w = 10^(att/10);
% win0 = chebwin(N0,att);
% win1 = chebwin(N1,att);
%  chebwin(n,r)����n��Ĵ��������丵��Ҷ�任����԰겨�Ʊ������rdB���԰��ǵȲ��Ƶġ�ע�⣺��nΪż��ʱ���������ĳ���Ϊ(n+1)��

%% Rectangular window
win0 = ones(N0,1);
win1 = ones(N1,1);
%% MARS_SFT
Win = win0*win1.';
%A�Ĺ���ת�þ���ΪA'��ת��ΪA.'
% [Omega, hA, P] = MARS_ISFT(Sig, Win, N0, N1, T, epsilon, gamma, n_d, n_s);
tic
[Omega, hA] = MARS_SFT_corr(Sig, Win, N0, N1, T, epsilon, gamma, gama);
% [Omega, hA] = MARS_SFT_vote(Sig, Win, N0, N1, T, epsilon, gamma, n_d, n_s);
% [Omega, hA] = MARS_SFT_nonoise(Sig, Win, N0, N1, T);
disp([' MARS-SFT��ʱ��Ϊ ',  num2str(toc)]);

%% reconstructed results
[Omega, ind] = sortrows(Omega); % ���Ƶ�λ�ã�����
%[B,index] = sortrows(A) Ĭ������A��һ�е���ֵ�������ƶ�ÿһ�У������һ�е���ֵ����ͬ�ģ��������ұȽϡ�
hA = hA(ind);   % ���Ƶ�Ƶ��ֵ
X=zeros(N0,N1);
u = Omega(:,1);%����λ�ã�
v = Omega(:,2);%����λ�ã�
X(sub2ind([N0,N1],u+1,v+1)) = hA;
% figure;mesh(0:N0-1,0:N1-1,abs(X));

%% Visualize
Omega_gt = round(Omega_gt); % ��������
visual_localization(N0,N1,Omega_gt,Omega);
% figure;spy(X);%ϡ�������ӻ�
error=norm(abs(X-X_true))/N
