%clc
%clear
%close
%% 2D chirp �ź�
K = 10; %ϡ���
a_max = 18.3;                     
sigma_n = 1; % ��������
delta = 0.1; 
% SNR_dB = 80; % �����
% SNR = 10^(SNR_dB/20);
% a_max = sqrt(sigma_n^2*SNR); % ����
fs=256;  %������Hz
T=1;   %�����������ʱ��
ts=1/fs;  % ������� 
n=round(T*fs); %�����������round��ȡ��
% x=linspace(0,T-1,n); %ʱ������
% y=linspace(0,T-1,n); %ʱ������
x=(0:n-1)*ts; X=ones(n,1)*x;% ʱ������
y=(0:n-1)*ts; Y=y'*ones(1,n);% ʱ������
A1=a_max;   A2=a_max;  A3=a_max;   % ����
M1= 197.9; M2= 28.7;  M3=218.8;   % �г�ʼƵ�� ��0��fs��֮�䣬�Һ�T�ĳ˻�������ʱ������������У���ȫϡ�裩��������й¶ 
N1=0.8;  N2=N1;   N3=N1; % �е�Ƶ��
R1= 224.2;  R2= 184.1; R3= 54.2;  % �г�ʼƵ��
S1=0.24;  S2=S1;   S3=S1;%-7.8 % �е�Ƶ��
F1=A1*exp(1j*pi*(2*M1*X+N1*X.^2+2*R1*Y+S1*Y.^2)); %�ź�1
F2=A2*exp(1j*pi*(2*M2*X+N2*X.^2+2*R2*Y+S2*Y.^2)); %�ź�2
F3=A3*exp(1j*pi*(2*M3*X+N3*X.^2+2*R3*Y+S3*Y.^2)); %�ź�2
F=F1+F2+F3;
% F = zeros(n,n);
% A = zeros(K,1);
% M = zeros(K,1);
% N = zeros(K,1);
% R = zeros(K,1);
% S = zeros(K,1);
% for ii = 1:K
%     A(ii) = a_max;% ����
%     M(ii) = (randi([0 n-1])+2*delta*rand(1)-delta)*fs/n;% �г�ʼƵ��
%     N(ii) = 0.8;% �е�Ƶ��
%     R(ii) = (randi([0 n-1])+rand(1)-0.5)*fs/n;% �г�ʼƵ��
%     S(ii) = 0.24;% �е�Ƶ��
%     F = F + A(ii)*exp(1j*pi*(2*M(ii)*X+N(ii)*X.^2+2*R(ii)*Y+S(ii)*Y.^2)); %�ź�1
% end
noise = sigma_n*sqrt(2)/2*(randn(n,n)+1i*randn(n,n));  % ʱ��Ӹ���˹����
SNR_t=snr(F,noise)
F = F+noise;%����ʱ���ź�
% F = awgn(F,15.587,'measured');
%   figure;imagesc(x,y,real(F));colormap(gray(256));colorbar% ��ӽ���ɫ��title('2D chirp�ź�');
%  figure; colormap("hsv");mesh((0:n-1)*T/n,(0:n-1)*T/n,real(F));zlabel('Amplitude');%axis([0 T 0 T]);title("Our methed");
%  colormap("winter");

%% 2D DFRFT
% �����ת��
% alfa=mod(acot(-N(1)),pi); p1=2*alfa/pi;
% beta=mod(acot(-S(1)),pi); p2=2*beta/pi;
alfa=mod(acot(-N1),pi); p1=2*alfa/pi;
beta=mod(acot(-S1),pi); p2=2*beta/pi;

% 2��1D DFRFT
% tic
% X_frft = DFRFT_2D_fft(F,p1,p2,ts,ts);
% disp([' ��fft���С��зֱ�DFRFT��ʱ��Ϊ ',  num2str(toc)]);
% figure;colormap("hot");mesh((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_frft)/n);zlabel('Amplitude');axis([0 fs 0 fs]);%title("Decomposition method");
% figure;colormap("colorcube");imagesc((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_frft)/n);  %��άͼ����ɫ��ǳ�������ȴ�С
  % 2D SFRFT ������Ƶ��
tic
[X_sfrft,Omega,us1,us2] = DFRFT_2D_sft(F,p1,p2,n,n,ts,ts,sigma_n);
disp(['sDFRFT��ʱ��Ϊ ',  num2str(toc)]);
%  figure;colormap("summer");mesh((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_sfrft)/n);zlabel('Amplitude');axis([0 fs 0 fs]);%title("Our methed");
 % figure;colormap("summer");imagesc((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_sfrft)/n);  %��άͼ����ɫ��ǳ�������ȴ�С
% colormap("hot");colormap("summer");colormap("colorcube");
%% λ�ÿ��ӻ� ������ʵλ�öԱȣ�
Omega_gt = [R1 M1
    R2 M2
    R3 M3]; 
Omega_gt = round(Omega_gt); % ��������
% visual_localization(n,n,Omega_gt,Omega*fs/n);
visual_localization(fs,fs,Omega_gt*n/fs,Omega);
