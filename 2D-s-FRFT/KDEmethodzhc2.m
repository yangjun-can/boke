clear
clc
%% ��������
N0=4096; % �źŴ�С
k=1;   % ϡ���
B=256;  % ��Ͱ��
D=N0/B; % ��������
SNR_dB = -10; % ʱ�� signal SNR in dB
SNR = 10^(SNR_dB/10);
sigma_t = 1; % ʱ�������ı�׼��
sigma_f = sqrt(N0)*sigma_t;  % Ƶ�������ı�׼��
A = sqrt(sigma_t^2*SNR);     % chirp�ź�ʱ��ķ���
A_f = N0*A;                  %�ź���Ƶ��ķ���

%A_con = 25*sqrt(2);
il=24; % Ƶ������
ql=3*il; %��������
jl=round(ql/D); % ��Ͱ����
lamda=10000;% ���������
%% ��λ��������
fai_e1=zeros(lamda,1);
for sample=1:lamda
    noise = sigma_f/sqrt(2)*(randn(N0,1)+1i*randn(N0,1)); %���ɸ���˹����
    zi=A_f;
    mu=A_f;
    for m=1:D
        zi=zi+noise(m);
        %mu=mu+noise(m)*exp(-1i*2*pi/N0*(m+D*jl-ql));
        mu=mu+noise(m)*exp(-1i*2*pi/N0*(m-D/2-1));
    end
    fai_e1(sample) =angle(zi/mu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = -5*2*pi/N0:5*2*pi/N0/100:5*2*pi/N0;
%figure('windowst','d');
figure;
binranges = theta;
[bincounts] = histc(fai_e1,binranges); % ͳ��fai_e1�ڸ�������binranges�ڵ�ֵ�ĸ���
bar(binranges*N0/(2*pi),bincounts/lamda,'histc');  % ��Ƶ��ͼ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h = histogram(fai_e1); % ���ֲ�ֱ��ͼ

%% ȷ��gama
for gama = 0:N0 % round(L/2/k)
    gamal = (gama+0.5)*2*pi/N0; % ���ֵ
    xi = -gamal:gamal/100:gamal;   % �����������Χ
    pdf_dis = ksdensity(fai_e1,xi); % ��Ӧ�ĸ����ܶ�
    prob = trapz(xi,pdf_dis);  % ���� ʹ��trapzģ�����
    if prob>0.90
        break;
    end
end

% figure;  plot(xi,pdf_dis);
figure;  plot(N0*xi/pi/2,pdf_dis*2*pi/N0);
disp(['  gama = ',num2str(gama)]);