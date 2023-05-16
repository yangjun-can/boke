% QFM信号参数估计
clear 
close all;clc;
fs=128;      % 采样频率
ts=1/fs;     % 采样间隔
N_all=512;   % 采样点数
T = N_all*ts; %信号采样时长
n=(0:N_all-1)-floor(N_all/2);
% n=(0:N_all-1);
t=n*ts;      % 离散时间，从-T/2~T/2
a1=-5;b1=-15;% QFM信号的参数  abs(a1)<64
a2=12;b2=13; % abs(a2)<16
a3=2;b3=-8;  % abs(a2)<2.666
A=1.24;
S1=A*exp(-1i*2*pi*(a1*t+a2*t.^2+a3*t.^3)); % 分量1
% S2=exp(1i*2*pi*(b1*t+0.5*b2*t.^2+1/6*b3*t.^3)); % 分量2
sig=S1;%+S2; % QFM信号
% % SNR=40;
% % sig=awgn(sig,SNR);
% % sig_fft=fft(sig);
figure; plot(real(sig));ylabel('Amplitude');xlabel('Sampling points')
% % figure; plot(abs(sig_fft));
Sig_dechirp=A*exp(-1i*2*pi*a1*t);
Sig_dechirp_fft=fftshift(fft(fftshift(Sig_dechirp)));
figure; plot(abs(Sig_dechirp_fft));ylabel('Amplitude');xlabel('Sampling points')


%% 信号的固定时延自相关函数，定义为s(t+tao) s*(t-tao) [s(t+tao+tao0) s*(t-tao-tao0)]*
tao0=T/256;
tao0=floor(tao0/ts);  %这里有一点近似
m1=floor(-N_all/2 : N_all/2-1);   %tao1的离散范围
tao_m1=m1*ts;    %时延tao1轴
N_tao1=length(m1);    %时延tao1的长度
[row,col] = meshgrid(1:N_all, 1:N_tao1);
sig(N_all+1) = 0;         %信号最后增加一个0点，为了后面处理
right1=row+m1(col);  % s(t+tao)
left1=row-m1(col);   % s(t-tao)
right2=row+m1(col)+tao0;  % s(t+tao+tao0)
left2=row-m1(col)-tao0;   % s(t-tao-tao0)
right1(right1 >N_all) = N_all+1;   %如果大于N_all，那就到信号最后一个点上去，而最后一个点设为了0
right1(right1 <1) = N_all+1;
left1(left1 >N_all) = N_all+1;
left1(left1 <1) = N_all+1;
right2(right2 >N_all) = N_all+1;   %如果大于N_all，那就到信号最后一个点上去，而最后一个点设为了0
right2(right2 <1) = N_all+1;
left2(left2 >N_all) = N_all+1;
left2(left2 <1) = N_all+1;
Rx = sig(right1) .* conj(sig(left1)) .* conj( sig(right2).* conj(sig(left2)) );
clear row col right1 left1 right2 left2;
figure;imagesc(t,tao_m1,abs(Rx)); 
xlabel('Time $t/s$',Interpreter='latex',FontSize=13);
ylabel('Lag $\tau/s$',Interpreter='latex',FontSize=13);

% Rx=zeros(N_all,N_all);
% for t=-N_all/2:N_all/2-1
%     for tao=-N_all/2:N_all/2-1
%         Rx(t+N_all/2+1,tao+N_all/2+1)=exp(-1j*2*pi*(2*a1*tao0+2*a3*tao0^3))*exp(-1j*2*pi*(6*a3*tao0^2*tao+6*a3*tao0*tao^2))*exp(-1j*2*pi*(4*a2*tao0*t+6*a3*tao0*t^2));
%     end
% end

%% lvd
% [ gama_axis, frequency_axis, lvd ] = LVD_CZT( Rx(:,245) );
% figure;mesh(gama_axis*fs^2, frequency_axis*fs, abs(lvd) );
%%  FRFT
alfa= mod(acot(12*a3*(tao0)*ts),pi);
p1=2*alfa/pi;
p2=p1;
% X_frft=DFRFT_2D_fft(Rx, p1,0,ts, ts); 
% figure;mesh(abs(X_frft));title("DFRFT");
X_frft = DFRFT_2D_fft(fftshift(Rx),p1,p2,ts,ts);
% figure;imagesc(fftshift(abs(X_frft)));title("DFRFT");
figure;mesh(fftshift(abs(X_frft)));title("DFRFT");

% X_frft = DFRFT_2D_fft(fftshift(Rx),p1,p2,ts,ts);
% figure;imagesc(fftshift(abs(X_frft)));title("DFRFT");

% [X_sfrft,us1,us2] = DFRFT_2D_sft(Rx,p1,p2,N_all,N_all,ts,ts,1);
% figure;colormap("summer");mesh(abs(X_sfrft));
% zlabel('Amplitude');%axis([0 fs 0 fs]);title("Our methed");

X_sfrft=zeros(N_all,N_all);
X_sfrft(4*a2*tao0,6*a3*(tao0)^2)=A^4;
figure;colormap("summer");mesh(abs(X_sfrft.'));zlabel('Amplitude');
xlabel('Samples of $t$',Interpreter='latex')
ylabel('Samples of $\tau$',Interpreter='latex')

%%======================SRA moving target detection================
% clear;clc;
%%========================================================
% %%Parameter--constant
% C=3e8;                           %propagation speed传播速度/光速
% %%Parameter--radar characteristics雷达特性
% Fc=1e9;                          %carrier frequency 1GHz 载频 1GHz
% lambda=C/Fc;                     %wavelength 波长
% %%Parameter--target area目标参数
% Xmin=0;                          %target area in azimuth is within[Xmin,Xmax]方位向的目标区域在 [Xmin,Xmax] 范围内
% Xmax=50;                  
% Yc=10000;                        %center of imaged area成像区域的中心
% Y0=500;                          %target area in range is within[Yc-Y0,Yc+Y0]距离向的目标区域在[Yc-Y0,Yc+Y0]内
%                                  %imaged width 2*Y0 成像宽度 2*Y0
% %%Parameter--orbital information轨道信息
% V=100;                           %SAR velosity 100 m/s SAR 速度 100 m/s
% H=5000;                          %height 5000 m 高度 5000 米
% R0=sqrt(Yc^2+H^2);
% %%Parameter--antenna 天线参数
% D=4;                            %antenna length in azimuth direction 方位向的天线长度
% Lsar=lambda*R0/D;               %SAR integration length  SAR长度
% Tsar=Lsar/V;                     %SAR integration time  SAR积分时间
% %%Parameter--slow-time domain 方位慢时域参数
% Ka=-2*V^2/lambda/R0;            %doppler frequency modulation rate 多普勒频率调制率
% Ba=abs(Ka*Tsar);                %doppler frequency modulation bandwidth  多普勒频率调制带宽
% PRF=2*Ba;                        %pulse repitition frequency 脉冲重复频率
% PRT=1/PRF;                     %pulse repitition time 脉冲重复时间
% ds=PRT;                         %sample spacing in slow-time domain 慢时域中的样本间距
% Nslow=ceil((Xmax-Xmin+Lsar)/V/ds); %sample number in slow-time domain 慢时域中的样本间距
% % Nslow=2^nextpow2(Nslow);          %for fft
% sn=linspace((Xmin-Lsar/2)/V,(Xmax+Lsar/2)/V,Nslow);%discrete time array in slow-time domain 慢时域离散时间数组
% PRT=(Xmax-Xmin+Lsar)/V/Nslow;    %refresh 刷新脉冲重复时间
% PRF=1/PRT;                        %refresh 刷新脉冲重复频率
% ds=PRT;                         % 更新慢时域中的样本间距
% %%Parameter--fast-time domain  距离快时域参数
% Tr=5e-6;                         %pulse duration 10us 脉冲持续时间 10us
% Br=30e6;                        %chirp frequency modulation bandwidth 30MHz 啁啾调频带宽 30MHz
% Kr=Br/Tr;                        %chirp slope 啁啾斜率
% Fsr=3*Br;                        %sampling frequency in fast-time domain 快时域采样频率
% dt=1/Fsr;                         %sample spacing in fast-time domain 快时域中的样本间距
% Rmin=sqrt((Yc-Y0)^2+H^2); 
% Rmax=sqrt((Yc+Y0)^2+H^2+(Lsar/2)^2);
% Nfast=ceil(2*(Rmax-Rmin)/C/dt+Tr/dt);%sample number in fast-time domain 快时域样本数
% % Nfast=2^nextpow2(Nfast);                   %for fft
% tm=linspace(2*Rmin/C,2*Rmax/C+Tr,Nfast); %discrete time array in fast-time domain 快时域离散时间数组
% dt=(2*Rmax/C+Tr-2*Rmin/C)/Nfast;    %refresh 更新快时域中的样本间距
% Fsr=1/dt;                           %refresh 更新快时域中的采样频率
% %%Parameter--resolution
% DY=C/2/Br;                          %range resolution 距离分辨率
% DX=D/2;                            %cross-range resolution 交叉距离分辨率
% %%Parameter--point targets 点目标
% Ntarget=1;                           %number of targets  目标个数
% %format [x, y, reflectivity]
% Ptarget=[Xmin,Yc,1
%         Xmin,Yc+10*DY,1
%         Xmin+20*DX,Yc+50*DY,1];  
% disp('参数:')
% disp('Sampling Rate in fast-time domain');disp(Fsr/Br)
% disp('快时间域采样点数');disp(Nfast)
% disp('Sampling Rate in slow-time domain');disp(PRF/Ba)
% disp('慢时间域采样点数');disp(Nslow)
% disp('距离分辨率');disp(DY)
% disp('方位分辨率');disp(DX)     
% disp('合成孔径长度');disp(Lsar)     
% disp('目标位置【（0-50），（9500-10500）】');disp(Ptarget)
% %%========================================================
% %%Generate the raw signal data  生成原始信号数据
% K=Ntarget;                                %number of targets 目标数量
% N=Nslow;                                   %number of vector in slow-time domain 慢时域中的向量数
% M=Nfast;                                   %number of vector in fast-time domain 快时域中的向量数
% T=Ptarget;                                %position of targets 目标位置
% Srnm1=zeros(N,M);
% for k=1:1:K
%     sigma=T(k,3);
%     Dslow=sn*V-T(k,1);
%     R=sqrt(Dslow.^2+T(k,2)^2+H^2);
%     tau=2*R/C;
%     Dfast=ones(N,1)*tm-tau'*ones(1,M);
%     phase=pi*Kr*Dfast.^2-(4*pi/lambda)*(R'*ones(1,M));
%     Srnm1=Srnm1+sigma*exp(j*phase).*(0<Dfast&Dfast<Tr).*((abs(Dslow)<Lsar/2)'*ones(1,M));
% end
% % row=tm*C/2-2008;col=sn*V-26;
% % figure;colormap(gray);
% % imagesc(row,col,abs(Srnm1));  
% % axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
% % xlabel('\rightarrow\itRange in meters'),ylabel('\itAzimuth in meters\leftarrow'),
% % title('Stripmap SAR befor range and azimuth compression'),
% %%========================================================
% %%Range compression 距离压缩
% tr=tm-2*Rmin/C;
% Refr=exp(j*pi*Kr*tr.^2).*(0<tr&tr<Tr);
% Sr=ifty(fty(Srnm1).*(ones(N,1)*conj(fty(Refr))));
% Gr=abs(Sr);
% %%Azimuth compression 方位压缩
% ta=sn-Xmin/V;
% Refa=exp(j*pi*Ka*ta.^2).*(abs(ta)<Tsar/2);
% Sa=iftx(ftx(Sr).*(conj(ftx(Refa)).'*ones(1,M)));
% Ga=abs(Sa);
% %%========================================================
% %%graw the intensity image of signal 绘制信号的强度图像
% 
% figure;colormap(gray);
% subplot(211);
% row=tm*C/2-2008;col=sn*V-26;
% imagesc(row,col,Gr);           %intensity image of Sr
% axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
% xlabel('\rightarrow\itRange in meters'),ylabel('\itAzimuth in meters\leftarrow'),
% title('Stripmap SAR after range compression'),
% subplot(212);
% imagesc(row,col,Ga);          %intensity image of Sa
% axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
% xlabel('\rightarrow\itRange in meters'),ylabel('\itAzimuth in meters\leftarrow'),
% title('Stripmap SAR after range and azimuth compression')

%%=====================利用二维FRFT的相位的位移不变性进行目标检测===================================

%step1 2D FRFT
% tic
% p1=0.0005;p2=0.0005;
% p1=0.00001
% alfa=mod(acot(2*aa(1)/lamda),pi); p2=2*alfa/pi;
%ts1=Ba/PRF;%ds;
%ts2=Br/Fsr;%dt;
% X_frft = DFRFT_2D_fft(data_r,p1,p2,ts,pri);
% disp([' 按fft对行、列分别DFRFT的时间为 ',  num2str(toc)]);
% figure;mesh(abs(X_frft));title("DFRFT");

% %step2 phase
% P=X_frft/abs(X_frft);
% %step3 reconstruction
% alpha=p1*pi/2;                   % 旋转角度（旋转角以2pi为周期）
% beta=p2*pi/2;                   % 旋转角度（旋转角以2pi为周期）
% us1=2*pi*sign(sin(alpha))*sin(alpha)/(ts1*Nslow);
% us2=2*pi*sign(sin(beta))*sin(beta)/(ts2*Nfast);
% Recons = DFRFT_2D_fft(P,-p1,-p2,us1,us2);
% % 画重构图
% 
% figure;colormap(gray);
% imagesc(row,col,abs(Recons));  
% axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
% xlabel('\rightarrow\itRange in meters');ylabel('\itAzimuth in meters\leftarrow');
% title('Detect by 2dDFRFT');
