clc
clear
%% ����ʱ���ϵ�1D���Ե�Ƶ�ź�
% n=256;
% fs=n;          % ����Ƶ��
% ts=1/fs;       % ������� 
% t1=(0:n-1)*ts; % ʱ������
% fd=[100 150 200 200];     % ��ʼƵ��
% mu=[-1.7 -1.4 -1.2 -2.5]; % ��Ƶ��
% x1=exp(1j*pi*mu(1)*(t1.^2)).*exp(1j*2*pi*fd(1)*t1);
% x2=exp(1j*pi*mu(2)*(t1.^2)).*exp(1j*2*pi*fd(2)*t1);
% x3=exp(1j*pi*mu(3)*(t1.^2)).*exp(1j*2*pi*fd(3)*t1);
% x4=exp(1j*pi*mu(4)*(t1.^2)).*exp(1j*2*pi*fd(4)*t1);
% x=x1+x2/2+x3/3+x4/4;      % ���Ե�Ƶ�źţ�ʱ��
%X=fftshift(fft(x)); X1=X/fs;   % ���Ե�Ƶ�źŷ�ֵ��Ƶ��

%% 2D chirp �ź�
fs=2048;  %������Hz
T=2^-4;   %�����������ʱ��
A1=0.5;  A2=0.5; A3=0.5;   % ����
M1= 1088; M2=240; M3=1824; %640 % �г�ʼƵ�� ��0��fs��֮�䣬�Һ�T�ĳ˻�������ʱ������������У���ȫϡ�裩��������й¶ 
N1=-0.8; N2=N1;  N3=N1;%-16.4% �е�Ƶ��
R1= 416;  R2=1600;R3= 1632;%896 % �г�ʼƵ��
S1=0.24; S2=S1;  S3=S1;%-7.8 % �е�Ƶ��
ts=1/fs;       % ������� 
n=round(T*fs); %�����������round��ȡ��
% x=linspace(0,T-1,n); %ʱ������
% y=linspace(0,T-1,n); %ʱ������
x=(0:n-1)*ts; % ʱ������
y=(0:n-1)*ts; % ʱ������
Lx=length(x);
Ly=length(y);
X=ones(Ly,1)*x;
Y=y'*ones(1,Lx);
F1=A1*exp(1j*pi*(2*M1*X+N1*X.^2+2*R1*Y+S1*Y.^2)); %�ź�1
F2=A2*exp(1j*pi*(2*M2*X+N2*X.^2+2*R2*Y+S2*Y.^2)); %�ź�2
F3=A3*exp(1j*pi*(2*M3*X+N3*X.^2+2*R3*Y+S3*Y.^2)); %�ź�2
F=F1+F2+F3;
% F = awgn(F,4.587,'measured');
%  figure;imagesc(x,y,real(F));title('2D chirp�ź�');colormap(gray(256));colorbar% ��ӽ���ɫ��
% figure; colormap("hsv");mesh((0:n-1)*T/n,(0:n-1)*T/n,real(F));zlabel('Amplitude');%axis([0 T 0 T]);title("Our methed");
% colormap("winter");

%% 2D DFRFT
% �����ת��
alfa=mod(acot(-N1),pi); p1=2*alfa/pi;
beta=mod(acot(-S1),pi); p2=2*beta/pi;
  % 2��1D DFRFT
tic
X_frft = DFRFT_2D_fft(F,p1,p2,ts,ts);
disp([' ��fft���С��зֱ�DFRFT��ʱ��Ϊ ',  num2str(toc)]);
figure;colormap("hot");mesh((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_frft)/n);zlabel('Amplitude');axis([0 fs 0 fs]);%title("Decomposition method");
figure;colormap("colorcube");imagesc((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_frft)/n);  %��άͼ����ɫ��ǳ�������ȴ�С
  % 2D SFRFT ������Ƶ��
% tic
% [X_sfrft,us1,us2] = DFRFT_2D_sft(F,p1,p2,n,n,ts,ts,1);
% disp(['sDFRFT��ʱ��Ϊ ',  num2str(toc)]);
% figure;colormap("summer");mesh((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_sfrft)/n);zlabel('Amplitude');axis([0 fs 0 fs]);%title("Our methed");
% figure;colormap("colorcube");imagesc((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_sfrft)/n);  %��άͼ����ɫ��ǳ�������ȴ�С
% % colormap("hot");colormap("summer");
   % MARS_SFT ������Ƶ��
% tic
% T=3;
% win0 = ones(n,1);
% win1 = ones(n,1);
% Win = win0*win1.';
% epsilon = 2e-6; %Threshold for frequency detectoin in a slice; set as 1e-10 for noiseless case
% gamma = 2e-6;
% [Omega, hA] = MARS_SFT_vote(F, Win, n, n, T, epsilon,gamma, 3, 2);
% [Omega, ind] = sortrows(Omega); % ���Ƶ�λ�ã�����
% hA = hA(ind)/n;   % ���Ƶ�Ƶ��ֵ
% X_SFT=zeros(n,n);
% u = Omega(:,1);%����λ�ã�
% v = Omega(:,2);%����λ�ã�
% X_SFT(sub2ind([n,n],u+1,v+1)) = hA;
% disp(['MARS_SFT��ʱ��Ϊ ',  num2str(toc)]);
% figure;colormap("summer");mesh((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_SFT)/n);zlabel('Amplitude');axis([0 fs 0 fs]);%title("Our methed");

