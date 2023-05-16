clear
clc
%% 构造随机信号
N=512; % 信号长度
k=2; % 信号稀疏度
SNR_dB =20; % signal SNR in dB
SNR = 20^(SNR_dB/20);
sigma_n = 0; % 噪声的方差 set this as 0 for noiseless signals, otherwise set as 1
a_min = 20^(SNR_dB/20); % the minimum signal power 信号最小振幅
[SigN,A,u] = genSigOffGrid(N,k,a_min,sigma_n);
%% 构造信号
%  N=512; % 信号长度
% J=6;
% k=2;  % 信号稀疏度
% K=2*N; % 上采样点数
%  t=0:N-1; % 采样点
% t2=0:K-1; % 上采样点
% A=[270,176]; % 信号幅值
% %  fu=[51.2/n,204.8/n]; %信号真实频率
% %   fu=[50.05/n,22.05/n];  %信号真实频率
% w=[0.2*pi,0.8*pi];
% SigN=zeros(1,N);
% SigK=zeros(1,K);
% for jj=1:k
%     %      w(jj)=2*pi*fu(jj);
%     O_x=A(jj)*sin(w(jj)*t);
%     SigN=SigN+O_x;
%     Up_x=A(jj)*sin(w(jj)*t2);
%     SigK=SigK+Up_x;
% end
%% NUDFT
for i=1:N
    f(i)=Nonuniform_sampling_point(i-1,N);
end
True_freq = ndftld(SigN, N, f);
%% Min_Max_NUFFT
tic;
%% 2001年
% % step 1
% UpSample_freq=fft(SigK);
% % step 2
% Interpolate_Set=zeros(1,J);
% W=zeros(n,J);
% b=zeros(n,1);
% X=zeros(J,1);
% freq=zeros(1,n);
% for Inter=1:n
%   [Interpolate_residual,Interpolate_Set]=close_to_freq(f(Inter),K,n,J);% find the points for interpolating
%   for poi=1:J
%      for len=1:n
%         W(len,poi)=exp(1i*2*pi*Interpolate_Set(poi)*(len-1)/K);
%         b(len)=exp(1i*2*pi*(len-1)*f(Inter)/n);
%      end
%       X(poi)=UpSample_freq(Interpolate_Set(poi));
%   end
%   P=inv(W.'*W)*W.'*b;
%   freq(Inter)=P'*X;
% end
%% 2003年
% step1:对上采样点加权并FFT
L=1;
K=2*N; % 上采样点数
J=6; % 用于插值的点的个数
beta=0.19;
galma=2*pi/K;
alfa=[conj(-0.46),1,-0.46];
eita=(N-1)/2;
Weights=alfa(1)*exp(-1i*galma*beta*((0:N-1)-eita))+alfa(2)+alfa(3)*exp(1i*galma*beta*((0:N-1)-eita));
Weights_SigN=Weights.*SigN;
UpSample_freq=fft(Weights_SigN,K);

%********定义法*********
% % step2 precomputed
% freq2=zeros(1,N);
% gama=zeros(J,J);
% S=zeros(N,N);
% C=zeros(N,J);
% b=zeros(N,1);
% for p1=1:N
%     S(p1,p1)=Weights(p1);   
%     for p2=1:J
%        C(p1,p2)=exp(1i*galma*p2*(p1-eita))/sqrt(N); 
%     end
% end
% % step2: 插值
% for point=1:N
%     % 定义用于插值point的第一个点
%     if mod(J,2)~=0
%         [Int_Offset,~]=close_to_freq(f(point),K,N,J);
%     else
%         Int_Offset=floor(f(point)*K/N)-J/2;  
%     end
%     for j2=1:J
%        gama(j2,j2)=exp(-1i*( 2*pi*f(point)/N -galma*(Int_Offset+j2) )*eita);  
%     end
%     for p3=1:N
%         b(p3)=exp(1i*( 2*pi*f(point)/N -galma*Int_Offset )*(p3-eita) ) /sqrt(N);          
%     end
%     u2=gama.'*inv(C.'*(S*S.')*C)*C.'*S*b;
%     for j3=1:J
%       freq2(point)= UpSample_freq(mod(Int_Offset+j3,K)+1) * conj(u2(j3)) +freq2(point);    
%     end
% end

%*********改进的Tr法（降低运算量）**************
% step2 precomputed
gama=zeros(J,J);
T=zeros(J,J);
r=zeros(J,1);
freq=zeros(1,N);
for cc=1:J
    for j1=1:J
T(cc,j1)=alfa(1)*( conj(alfa(1))*Periodic_sinc(j1-cc,K,N) + conj(alfa(2))*Periodic_sinc(j1-cc-beta,K,N) + conj(alfa(3))*Periodic_sinc(j1-cc-2*beta,K,N) ) ...
        +alfa(2)*( conj(alfa(1))*Periodic_sinc(j1-cc+beta,K,N) + conj(alfa(2))*Periodic_sinc(j1-cc,K,N) + conj(alfa(3))*Periodic_sinc(j1-cc-beta,K,N) ) ...
        +alfa(3)*( conj(alfa(1))*Periodic_sinc(j1-cc+2*beta,K,N) + conj(alfa(2))*Periodic_sinc(j1-cc+beta,K,N) + conj(alfa(3))*Periodic_sinc(j1-cc,K,N) );
    end
end
% step3: 插值
for point=1:N
    % 定义用于插值point的第一个点
    if mod(J,2)~=0
        [Int_Offset,~]=close_to_freq(f(point),K,N,J);
    else
        Int_Offset=floor(f(point)*K/N)-J/2;  
    end
    for j2=1:J
       r(j2)=alfa(1)* Periodic_sinc( f(point)*K/N-Int_Offset-j2-beta ,K,N) ...
           + alfa(2)* Periodic_sinc( f(point)*K/N-Int_Offset-j2,K,N ) ...
           + alfa(3)* Periodic_sinc( f(point)*K/N-Int_Offset-j2+beta,K,N );
       gama(j2,j2)=exp(-1i*( 2*pi*f(point)/N -galma*(Int_Offset+j2) )*eita);  
    end
    u=gama.'*inv(T)*r;  %插值的系数向量
    for j3=1:J
      freq(point)= UpSample_freq(mod(Int_Offset+j3,K)+1) * conj(u(j3)) +freq(point);
    end
end   
toc;
figure;
% subplot(311);plot(abs(UpSample_freq));xlabel('Up sampling point in the frequency domain');ylabel('amplitude');
subplot(311);plot(f,abs(True_freq));xlabel('Ununiform sampling point in the frequency domain');ylabel('amplitude');
% subplot(312);plot(f,abs(freq2));xlabel('Ununiform Estimate point in the frequency domain by interpolator defination');ylabel('amplitude');
subplot(313);plot(f,abs(freq));xlabel('Ununiform Estimate point in the frequency domain by Tr precomputed');ylabel('amplitude');

%     %% 相对误差
%     Min_Max_err=0;
%     for ii=1:N
%      
%         Min_Max_err = Min_Max_err+abs(abs(True_freq(ii))-abs(freq(ii)))/abs(True_freq(ii));
%         
%     end
%     Min_Max_err=Min_Max_err/N;
%     disp(['L1 Error is ',  num2str(Min_Max_err)]);
