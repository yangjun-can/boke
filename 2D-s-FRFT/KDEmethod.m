clear
%% 参数设置
N1_all=[64	128	256	512	1024 2048, 4096]; % 信号大小
N2_all=[64	128	256	512	1024 2048, 4096];
SNR_all= [71	69	68	67	66	65];

N1=2048;
N2=2048;
 %for index=1:length(SNR_all)

SNR_dB =43;       % 最小频率信噪比

N = N1*N2;
k=1;   % 稀疏度
kesail=76; % 大频率坐标
kesai2=75;
L = lcm(N1,N2); % line length
D=N/L; % 噪声点数
sigma_t = 1; % 时域噪声的标准差
sigma_f = sqrt(N)*sigma_t;  % 傅里叶域噪声的标准差

SNR = 10^(SNR_dB/10);
A_f = sqrt(sigma_f^2*SNR); % 最小幅度
A_f1 = sqrt(sigma_f^2*SNR)*exp(1i*2*pi*rand(1,1)); % 频率值
% A = sqrt(sigma_t^2*SNR);     % chirp信号时域的幅度
% A_f = N*A;
%% 构造随机切片的参数
tau1 = randi(N1,1)-1;
tau2 = randi(N2,1)-1;
alpha1 = randi(N1,1)-1;
alpha2 = randi(N2,1)-1;
% randi(imax,n)产生一个n*n矩阵，这个矩阵的元素都是小于等于imax的正整数。
while gcd(alpha1,alpha2) ~= 1 || gcd(alpha1, L/N2) ~= 1 || gcd(alpha2,L/N1) ~= 1
    % gcd 求两个整数的最大公约数，返回值是正整数
    alpha1 = randi(N1,1)-1;
    alpha2 = randi(N2,1)-1;
end
% h=mod(L*kesail*alpha1/N1+L*kesai2*alpha2/N2,L); % 大桶坐标
%% 构造样本 （计算的慢）
% l = 0:D-1;
% x = mod(kesail+l*alpha2*L/N2,N1)+1;
% y = mod(kesai2-l*alpha1*L/N1,N2)+1;
% ind = sub2ind([N1,N2], x, y);  % 返回位置（x，y）在大小[N0,N1]的矩阵中的按列索引号
% m1=0:N1-1;
% m2=0:N2-1;
% chirp_1 = exp(1i*2*pi*m1*tau1/N1);
% chirp_11 = exp(1i*2*pi*(m1-kesail)/N1);
% chirp_2 = exp(1i*2*pi*m2*tau2/N2);
% lamda=10000;% 随机样本数
% fai_e1=zeros(lamda,1); % 相位误差的样本
% for sample=1:lamda
%     noise = sigma_f/sqrt(2)*(randn(N1,N2)+1i*randn(N1,N2)); %生成复高斯噪声
%     n_mu=noise.*repmat(chirp_1.',1,N2).*repmat(chirp_2,N1,1);
%     n_zi=n_mu.*repmat(chirp_11.',1,N2);
%     zi=A_f*exp(1i*2*pi*(kesail*tau1/N1+kesai2*tau2/N2))+sum(n_zi(ind));
%     mu=A_f*exp(1i*2*pi*(kesail*tau1/N1+kesai2*tau2/N2))+sum(n_mu(ind));
%     fai_e1(sample) =angle(zi/mu);
% end
%% 构造样本
% l = 0:D-1;
% x = mod(kesail+l*alpha2*L/N2,N1);
% y = mod(kesai2-l*alpha1*L/N1,N2);
% chirp_1 = exp(1i*2*pi*x*tau1/N1);
% chirp_11 = exp(1i*2*pi*(x-kesail)/N1);
% chirp_12 = exp(1i*2*pi*(y-kesai2)/N1);
% chirp_2 = exp(1i*2*pi*y*tau2/N2);
% lamda=10000;% 随机样本数
% fai_e1=zeros(lamda,1); % 相位误差的样本
% for sample=1:lamda
%     noise = sigma_f/sqrt(2)*(randn(1,D)+1i*randn(1,D)); %生成复高斯噪声
%     n_mu=noise.*chirp_1.*chirp_2;
%     n_zi=n_mu.*chirp_12;
%     zi=A_f1*exp(1i*2*pi*(kesail*tau1/N1+kesai2*tau2/N2))+sum(n_zi);
%     mu=A_f1*exp(1i*2*pi*(kesail*tau1/N1+kesai2*tau2/N2))+sum(n_mu);
%     fai_e1(sample) =angle(zi/mu);
% end
%% 画频率图
% theta = -5*2*pi/N1:5*2*pi/N1/100:5*2*pi/N1;
% %figure('windowst','d');
% figure;
% binranges = theta;
% [bincounts] = histc(fai_e1,binranges); % 统计fai_e1在给定区间binranges内的值的个数
% bar(binranges*N1/(2*pi),bincounts/lamda,'histc');  % 画频率图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  h = histogram(fai_e1); % 画分布直方图

%% 确定gama
% for gama = 0:N1 % round(L/2/k)
%     gamal = (gama+0.5)*2*pi/N1; % 最大值
%     xi = -gamal:gamal/7000:gamal;   % 随机变量容许范围
%     pdf_dis = ksdensity(fai_e1,xi); % 相应的概率密度
%     prob = trapz(xi,pdf_dis);  % 概率 使用trapz模拟积分
%     if prob>0.999
%         break;
%     end
% end
% figure;  plot(xi,pdf_dis); xlabel('${\phi _1^{err}}$','Interpreter','latex');ylabel('Probability density');grid
% figure; plot(N1*xi/pi/2,pdf_dis/N1/pi);xlabel('m_1^{err}');ylabel('Probability density');grid
% disp(['  gama = ',num2str(gama)]);
 %end 
 %% 一稀疏检测
miu1 = mod(kesail+5*alpha2*L/N2,N1); %频率2位置
miu2 = mod(kesai2-5*alpha1*L/N1,N2);
A_f2 = sqrt(sigma_f^2*SNR)*exp(1i*2*pi*rand(1,1)); % 频率值
% 一稀疏桶
L10=A_f1*exp(1i*2*pi*(kesail*tau1/N1+kesai2*tau2/N2))*L/N;
L11=L10*exp(1i*2*pi*kesail/N1);
L12=L10*exp(1i*2*pi*kesai2/N2);
% 非一稀疏，另一大频率
L20=A_f2*exp(1i*2*pi*(miu1*tau1/N1+miu2*tau2/N2))*L/N;
L21=L20*exp(1i*2*pi*miu1/N1);
L22=L20*exp(1i*2*pi*miu2/N2);
% 噪声
 l=0:D-1;
 noisex = mod(kesail+l*alpha2*L/N2,N1); 
 noisey = mod(kesai2-l*alpha1*L/N1,N2);
 noise = sigma_f/sqrt(2)*(randn(1,D)+1i*randn(1,D)); %生成复高斯噪声
%  delta=4*L/N*sum(abs(noise))^2
% 噪声的投影
 noise0 = noise.*exp(1i*2*pi*noisex*tau1/N1).*exp(1i*2*pi*noisey*tau2/N2)*L/N; 
 noise1 = noise0.*exp(1i*2*pi*noisex/N1);
 noise2 = noise0.* exp(1i*2*pi*noisey/N2);
%% 无躁时
%  % 定位行位置
% true1_loc=round(wrapTo2Pi(angle(L11/L10))*N1/2/pi);
% fouse1_loc=round(wrapTo2Pi(angle( (L11+L21)/(L10+L20) ) *N1/2/pi));
% e1=fouse1_loc-true1_loc;
% % 定位列位置
% true2_loc=round( angle(L12/L10 )  *N2/2/pi);
% fouse2_loc=round( angle((L12+L22)/(L10+L20)) *N2/2/pi);
% e2=fouse2_loc-true2_loc;
% % 1稀疏检测方法
% detc1_ture=abs(L11-L10*exp(1i*2*pi*true1_loc/N1))^2
% detc1_foulse=abs(L11+L21-(L10+L20)*exp(1i*2*pi*fouse1_loc/N1))^2
% 
% detc2_ture=abs(L11-L10)^2
% detc2_foulse=abs(L21+L11-L20-L10)^2
% 
% 
% bin1_mod=[abs(L10),abs(L11),abs(L12)];
% detc3_ture=var(bin1_mod)
% bin2_mod=[abs(L10+L20),abs(L11+L21),abs(L12+L22)];
% detc3_foulse=var(bin2_mod)
% 
% detc4_ture=abs(bin1_mod(1)-bin1_mod(2))^2
% detc4_foulse=abs(bin2_mod(1)-bin2_mod(2))^2

%% 含躁时
% % % 定位行位置
% % true1_loc_noise=round(wrapTo2Pi(angle((L11+sum(noise1))/(L10+sum(noise0))))*N1/2/pi);
% % fouse1_loc_noise=round(wrapTo2Pi(angle( (L11+L21+sum(noise1))/(L10+L20+sum(noise0)) ) *N1/2/pi));
% % e1_noise=fouse1_loc_noise-true1_loc_noise;
% % % 定位列位置
% % true2_loc_noise=round( angle((L12+sum(noise2))/(L10+sum(noise0) ))  *N2/2/pi);
% % fouse2_loc_noise=round( angle((L12+L22+sum(noise2))/(L10+L20+sum(noise0))) *N2/2/pi);
% % e2_noise=fouse2_loc_noise-true2_loc_noise;
% % 
% % 1稀疏检测方法
% % detc1_ture_noise=abs(L11+sum(noise1)-(L10+sum(noise0))*exp(1i*2*pi*true1_loc_noise/N1))^2
% % detc1_foulse_noise=abs(L11+L21+sum(noise1)-(L10+L20+sum(noise0))*exp(1i*2*pi*fouse1_loc_noise/N1))^2
% % 
% % detc2_ture_noise=abs(L11+sum(noise1)-L10-sum(noise0))^2
% % detc2_foulse_noise=abs(L21+L11+sum(noise1)-L20-L10-sum(noise0))^2
% 
bin1_mod_noise=[abs(L10+sum(noise0)),abs(L11+sum(noise1)),abs(L12+sum(noise2))];
bin2_mod_noise=[abs(L10+L20+sum(noise0)),abs(L11+L21+sum(noise1)),abs(L12+L22+sum(noise2))];
detc3_ture_noise=var(bin1_mod_noise)
detc3_foulse_noise=var(bin2_mod_noise)
% detc4_ture_noise=(bin1_mod_noise(1)-bin1_mod_noise(2))^2
% detc4_foulse_noise=(bin2_mod_noise(1)-bin2_mod_noise(2))^2

%% 大桶检测
% 空桶
bin0=abs(sum(noise0))^2
% 1稀疏桶
bin1=abs(L10+sum(noise0))^2
% 阈值
% Pd=0.9;
% Pfa=Pd^(L*SNR+1);
Pfa=2e-3;
epsilon = -2*L*log(Pfa)
% 








%% KDE工具箱估计分布
%  pdf = kde(fai_e1,'rot',[],'E');figure; plot(pdf,'r-');
% lcv 根据最大似然估计自动确定带宽
% hall,hjsm 渐进性较好的插值估计器，根据MISE准则自动确定带宽
% rot,msp 基于标准偏差的快速方法，根据AMISE准则自动确定带宽
%% 按KDE定义估计分布
% [t, y_true, tt, y_KDE] = KernelDensityEstimation(fai_e1);
%% Kernel Density Estimation % 只能处理正半轴密度
function [t, y_true, tt, y_KDE] = KernelDensityEstimation(x)
%参数初始化
Max = round(max(x));           %数据中最大值
Min = round(min(x));           %数据中最小值
Ntotal = length(x);     %数据个数
tt = 0 : 0.1 : Max;     %精确x轴
t = 0 : Max;            %粗略x轴
y_KDE = zeros(10 * Max+1, 1);   %核密度估计值
sum1 = 0;                       %求和的中间变量
%计算带宽h
R = 1/(2*sqrt(pi));
m2 = 1;
% h = 3;
h = (R)^(1/5) / (m2^(2/5) * R^(1/5) * Ntotal^(1/5));
%%
%计算核密度估计
for i = 0 : 0.1 : Max
    for j = 1 : Ntotal
        sum1 = sum1 + normpdf(i-x(j));
    end
    y_KDE(round(i*10+1)) = sum1 / (h * Ntotal);
    sum1 = 0;
end

sum2 = sum(y_KDE)*0.1;  %归一化KDE密度
for i = 0 : 0.1 : Max
    y_KDE(round(i*10+1)) = y_KDE(round(i*10+1))/sum2;
end
%%
%计算真实密度的分布
y_true = zeros(Max+1,1);
for i = 0 : Max
    for j = 1 : Ntotal
        if (x(j) < i+1)&&(x(j) >= i)
            y_true(i+1) = y_true(i+1) + 1;
        end
    end
    y_true(i+1) = y_true(i+1) / Ntotal;
end
%% 绘图
figure;          %真实密度的分布图象
bar(t, y_true);
axis([Min Max+1 0 max(y_true)*1.1]);
figure;          %核密度估计的密度分布图象
plot(tt, y_KDE);
axis([Min Max 0 max(y_true)*1.1]);
end
