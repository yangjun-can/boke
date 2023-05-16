clc
clear
%% signal in time domain
% Analog signal [0,T]
T=30;         
Dt = 0.005;   % 栅格  很小，用来近似连续的模拟信号
t = 0 : Dt : T ;
[x_analog_time,~,~ ]= signal(t); % 模拟信号

%% Uniform Discrete-time signal
ts = 1;  % 采样间隔
uniform_sample_time= 0 : ts : T; % 采样序列
Fs = 1/ts;    % 采样频率
N=length(uniform_sample_time);% 采样点数      
% n=0:N-1; % uniform_sample_time = n*Ts;  % 采样序列
[x_Udis_time,HSx_Udis_time,WSx_Udis_time] = signal(uniform_sample_time); % 均匀采样的离散信号

%% NonUniform Discrete-time signal
deviation=0.01; % 采样扰动标准差
kesai=deviation*randn(1,N); % 采样扰动
Non_uniform_sample_time= uniform_sample_time + kesai; % 非均匀采样序列
[x_NUdis_time,HSx_NUdis_time,WSx_NUdis_time] = signal(Non_uniform_sample_time); % 非均匀的离散信号
% Half-sample symmetric
HSN=2*N;   % 半对称扩展的长度
sy_kesai=fliplr(kesai); %左右对称
HSkesai=cat(2,sy_kesai,kesai);  % 采样扰动的半对称扩展 按列将sy_kesai,kesai连接起来
HS_sample=HSkesai+(0:HSN-1);    % 半对称扩展后的采样点
% whole-sample symmetric
WSN=2*N-1; % 全对称扩展的长度
WSkesai=HSkesai;
WSkesai(N+1)=[];             % 采样扰动的全对称扩展
WS_sample=WSkesai+(0:WSN-1); % 全对称扩展后的采样点
%% signal in LCTD
%  a=0.1;b=0.2;c=0;d=10;  % LCT参数
a=cos(1.5394);  b=sin(1.5394); c=-sin(1.5394); d=cos(1.5394);
% a=3;b=1;c=1/2;d=1/2;
%% LCT of Uniform Discrete signal
[X_Udis_lct,us] = LCT(x_Udis_time,a,b,c,d, N, ts); % 均匀离散信号的LCT
[HSX_Udis_lct,HSus] = LCT(HSx_Udis_time,a,b,c,d, HSN, ts); % 半对称扩展的LCT
[WSX_Udis_lct,WSus] = LCT(WSx_Udis_time,a,b,c,d, WSN, ts); % 全对称扩展的LCT
% Rex_Udis_time=ILCT(X_Udis_lct,a,b,c,d, N, us);
% ReHSx_Udis_time=ILCT(HSX_Udis_lct,a,b,c,d, HSN, us);
% ReWSx_Udis_time=ILCT(WSX_Udis_lct,a,b,c,d, WSN, us);
% Rex_Udis_time=LCT(X_Udis_lct,d,-b,-c,a, N, us);
% ReHSx_Udis_time=LCT(HSX_Udis_lct,d,-b,-c,a, HSN, us);
% ReWSx_Udis_time=LCT(WSX_Udis_lct,d,-b,-c,a, WSN, us);
% %检验对称扩展的效果，不含非均匀采样
% figure; % LCT
% subplot(3,2,1); stem(0:N-1,real(X_Udis_lct)); title("Original real LCT") 
% subplot(3,2,2); stem(0:N-1,imag(X_Udis_lct)); title("Original imag LCT") 
% subplot(3,2,3); stem(0:HSN-1,real(HSX_Udis_lct)); title("HS real LCT") 
% subplot(3,2,4); stem(0:HSN-1,imag(HSX_Udis_lct)); title("HS imag LCT") 
% subplot(3,2,5); stem(0:WSN-1,real(WSX_Udis_lct)); title("WS real LCT") 
% subplot(3,2,6); stem(0:WSN-1,imag(WSX_Udis_lct)); title("WS imag LCT") 
% figure; % time
% subplot(3,2,1)
% stem(0:N-1,real(Rex_Udis_time),'.');hold on 
% stem(0:N-1,real(x_Udis_time));hold off 
% title("Original real") 
% subplot(3,2,2)
% stem(0:N-1,imag(Rex_Udis_time),'.');hold on 
% stem(0:N-1,imag(x_Udis_time));hold off 
% title("Original imag")
% subplot(3,2,3)
% stem(0:HSN-1,real(ReHSx_Udis_time),'.');hold on 
% stem(0:HSN-1,real(HSx_Udis_time));hold off 
% title("HS real")
% subplot(3,2,4)%(N+1:HSN)
% stem(0:HSN-1,imag(ReHSx_Udis_time),'.');hold on 
% stem(0:HSN-1,imag(HSx_Udis_time));hold off 
% title("HS imag")
% subplot(3,2,5)%(N:WSN)
% stem(0:WSN-1,real(ReWSx_Udis_time),'.');hold on 
% stem(0:WSN-1,real(WSx_Udis_time));hold off 
% title("WS real")
% subplot(3,2,6)
% stem(0:WSN-1,imag(ReWSx_Udis_time),'.');hold on 
% stem(0:WSN-1,imag(WSx_Udis_time));hold off 
% title("WS imag")

%% LCT of Non_Uniform Discrete signal
[X_NUdis_lct,~] = LCT(x_NUdis_time,a,b,c,d, N, ts); % 非均匀离散信号的LCT
[HSX_NUdis_lct,~] = LCT(HSx_NUdis_time,a,b,c,d, HSN, ts);% 半对称扩展的LCT
[WSX_NUdis_lct,~] = LCT(WSx_NUdis_time,a,b,c,d, WSN, ts);% 全对称扩展的LCT

%% reconstruction Uniform Discrete signal
%% in LCT
ReX_Udis_lct=reconstruction_Un_LCT(X_NUdis_lct,a,b,c,d,us,ts,kesai,Non_uniform_sample_time,N);
ReHSX_Udis_lct=reconstruction_Un_LCT(HSX_NUdis_lct,a,b,c,d,HSus,ts,HSkesai,HS_sample,HSN);
ReWSX_Udis_lct=reconstruction_Un_LCT(WSX_NUdis_lct,a,b,c,d,WSus,ts,WSkesai,WS_sample,WSN);
%% in time
Rex_Udis_time=LCT(ReX_Udis_lct,d,-b,-c,a, N, us);
ReHSx_Udis_time=LCT(ReHSX_Udis_lct,d,-b,-c,a, HSN, HSus);
ReWSx_Udis_time=LCT(ReWSX_Udis_lct,d,-b,-c,a, WSN, WSus);
%% plot 
% 无对称扩展
% figure; 
% % subplot(3,1,1); plot(t,real(x_analog_time)); title('Analog signal'); xlabel('t in msec'); ylabel('xa');
% subplot(2,2,1);
% stem(uniform_sample_time,real(Rex_Udis_time),'x');hold on 
% stem(uniform_sample_time,real(x_Udis_time));hold off
% title('时域实部');
% subplot(2,2,2);
% stem(uniform_sample_time,imag(Rex_Udis_time),'x');hold on 
% % stem(uniform_sample_time,imag(ReHSx_Udis_time(N+1:HSN)),'.');hold on 
% stem(uniform_sample_time,imag(x_Udis_time));hold off
% title('时域虚部');
% subplot(2,2,3);
% stem(0:N-1,real(ReX_Udis_lct),'x');hold on 
% % stem(0:N-1,real(ReHSX_Udis_lct(1:2:HSN)),'.');hold on 
% stem(0:N-1,real(X_Udis_lct));hold off
% title('LCT域虚部');
% subplot(2,2,4); 
% stem(0:N-1,imag(ReX_Udis_lct),'x');hold on
% % stem(0:N-1,imag(ReHSX_Udis_lct(1:2:HSN)),'.');hold on 
% stem(0:N-1,imag(X_Udis_lct));hold off
% title('LCT域实部');

% % 检验对称扩展效果 含非均匀    
% figure;  
% subplot(2,3,1);  
% stem(0:HSN-1,real(ReHSx_Udis_time),'.');hold on
% stem(0:HSN-1,real(HSx_Udis_time));hold off
% title('时域实部');
% subplot(2,3,2); 
% stem(0:HSN-1,imag(ReHSx_Udis_time),'.');hold on
% stem(0:HSN-1,imag(HSx_Udis_time));hold off
% title('时域虚部');
% subplot(2,3,3); 
% stem(0:HSN-1,abs(ReHSx_Udis_time),'.');hold on
% stem(0:HSN-1,abs(HSx_Udis_time));hold off
% title('时域幅度');
% 
% subplot(2,3,4);  
% stem(0:HSN-1,real(ReHSX_Udis_lct),'.');hold on
% stem(0:HSN-1,real(HSX_Udis_lct));hold off
% title('LCT域实部');
% subplot(2,3,5); 
% stem(0:HSN-1,imag(ReHSX_Udis_lct),'.');hold on
% stem(0:HSN-1,imag(HSX_Udis_lct));hold off
% title('LCT域虚部');
% subplot(2,3,6); 
% stem(0:HSN-1,abs(ReHSX_Udis_lct),'.');hold on
% stem(0:HSN-1,abs(HSX_Udis_lct));hold off
% title('LCT域幅度');

% 三种重构效果对比(时域) 
figure;  
subplot(3,1,1); 
stem(0:N-1,real(Rex_Udis_time),'.','filled');hold on
stem(0:N-1,real(ReHSx_Udis_time(N+1:HSN)),'*');hold on
stem(0:N-1,real(ReWSx_Udis_time(N:WSN)),'x','c');hold on
% plot(t,x_analog_time);hold on
stem(0:N-1,real(x_Udis_time));hold off
legend('Reconstruction without extensions','Reconstruction with HS','Reconstruction with WS','True uniform samples');
title('时域实部');
subplot(3,1,2); 
stem(0:N-1,imag(Rex_Udis_time),'.','filled');hold on
stem(0:N-1,imag(ReHSx_Udis_time(N+1:HSN)),'*');hold on
stem(0:N-1,imag(ReWSx_Udis_time(N:WSN)),'x','c');hold on
% plot(t,x_analog_time);hold on
stem(0:N-1,imag(x_Udis_time));hold off
legend('Reconstruction without extensions','Reconstruction with HS','Reconstruction with WS','True uniform samples');
title('时域虚部');
subplot(3,1,3); 
stem(0:N-1,abs(Rex_Udis_time),'.','filled');hold on
stem(0:N-1,abs(ReHSx_Udis_time(N+1:HSN)),'*');hold on
stem(0:N-1,abs(ReWSx_Udis_time(N:WSN)),'x','c');hold on
% plot(t,x_analog_time);hold on
stem(0:N-1,abs(x_Udis_time));hold off
legend('Reconstruction without extensions','Reconstruction with HS','Reconstruction with WS','True uniform samples');
title('时域幅度');

figure;% 三种重构效果对比(LCT域)
subplot(3,3,1);  
stem(0:N-1,real(ReX_Udis_lct),'.');hold on
stem(0:N-1,real(X_Udis_lct));hold off
legend('Reconstruction','True DLCT');
title('LCT域实部');
subplot(3,3,2); 
stem(0:N-1,imag(ReX_Udis_lct),'.');hold on
stem(0:N-1,imag(X_Udis_lct));hold off
legend('Reconstruction','True DLCT');
title('LCT域虚部');
subplot(3,3,3); 
stem(0:N-1,abs(ReX_Udis_lct),'.');hold on
stem(0:N-1,abs(X_Udis_lct));hold off
legend('Reconstruction','True DLCT');
title('LCT域幅度');

subplot(3,3,4);  
stem(0:HSN-1,real(ReHSX_Udis_lct),'.');hold on
stem(0:HSN-1,real(HSX_Udis_lct));hold off
legend('Reconstruction of HS','True DLCT of HS');
title('LCT域实部');
subplot(3,3,5); 
stem(0:HSN-1,imag(ReHSX_Udis_lct),'.');hold on
stem(0:HSN-1,imag(HSX_Udis_lct));hold off
title('LCT域虚部');
subplot(3,3,6); 
stem(0:HSN-1,abs(ReHSX_Udis_lct),'.');hold on
stem(0:HSN-1,abs(HSX_Udis_lct));hold off
title('LCT域幅度');

subplot(3,3,7);  
stem(0:WSN-1,real(ReWSX_Udis_lct),'.');hold on
stem(0:WSN-1,real(WSX_Udis_lct));hold off
legend('Reconstruction of WS','True DLCT of WS');
title('LCT域实部');
subplot(3,3,8); 
stem(0:WSN-1,imag(ReWSX_Udis_lct),'.');hold on
stem(0:WSN-1,imag(WSX_Udis_lct));hold off
legend('Reconstruction of WS','True DLCT of WS');
title('LCT域虚部');
subplot(3,3,9); 
stem(0:WSN-1,abs(ReWSX_Udis_lct),'.');hold on
stem(0:WSN-1,abs(WSX_Udis_lct));hold off
legend('Reconstruction of WS','True DLCT of WS');
title('LCT域幅度');
