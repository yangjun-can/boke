clear
clc
%% 参数设置
N0=4096; % 信号大小
k=1;   % 稀疏度
B=256;  % 分桶数
D=N0/B; % 噪声点数
SNR_dB = -10; % 时域 signal SNR in dB
SNR = 10^(SNR_dB/10);
sigma_t = 1; % 时域噪声的标准差
sigma_f = sqrt(N0)*sigma_t;  % 频域噪声的标准差
A = sqrt(sigma_t^2*SNR);     % chirp信号时域的幅度
A_f = N0*A;                  %信号在频域的幅度

%A_con = 25*sqrt(2);
il=24; % 频率坐标
ql=3*il; %重排坐标
jl=round(ql/D); % 大桶坐标
lamda=10000;% 随机样本数
%% 相位误差的样本
fai_e1=zeros(lamda,1);
for sample=1:lamda
    noise = sigma_f/sqrt(2)*(randn(N0,1)+1i*randn(N0,1)); %生成复高斯噪声
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
[bincounts] = histc(fai_e1,binranges); % 统计fai_e1在给定区间binranges内的值的个数
bar(binranges*N0/(2*pi),bincounts/lamda,'histc');  % 画频率图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h = histogram(fai_e1); % 画分布直方图

%% 确定gama
for gama = 0:N0 % round(L/2/k)
    gamal = (gama+0.5)*2*pi/N0; % 最大值
    xi = -gamal:gamal/100:gamal;   % 随机变量容许范围
    pdf_dis = ksdensity(fai_e1,xi); % 相应的概率密度
    prob = trapz(xi,pdf_dis);  % 概率 使用trapz模拟积分
    if prob>0.90
        break;
    end
end

% figure;  plot(xi,pdf_dis);
figure;  plot(N0*xi/pi/2,pdf_dis*2*pi/N0);
disp(['  gama = ',num2str(gama)]);