% clc
clear
n=256;      % 信号长度
k=2;        % 信号稀疏度
%% 构造时域上的线性调频信号
fs=n;          % 采样频率
ts=1/fs;       % 采样间隔
t1=(0:n-1)*ts; % 时间向量  共1s
tao=(0:n-1)*ts;
fd=[20,40];   % 初始频率
mu=[-10,30] ; % 调频率
x1=exp(1j*pi*mu(1)*(t1.^2)).*exp(1j*2*pi*fd(1)*t1);
x2=exp(1j*pi*mu(2)*(t1.^2)).*exp(1j*2*pi*fd(2)*t1);
s=x1+x2;      % 线性调频信号（时域）
%% 瞬时自相关函数
%定义
% t1_all=repmat(t1.', 1, n);
% tao_all=repmat(tao, n, 1);
% Rs=s(t1_all+tao_all).*conj(s(t1_all+tao_all));
Rs=zeros(n,n);
for i = 1:k
    s1=exp(1j*4*pi*fd(i)*tao);
    s1_all=repmat(s1, n, 1);
    s2=exp(1j*4*pi*mu(i)*t1.'*tao);
    rs=s1_all.*s2;
    Rs=Rs+rs;
end
%% DAT
S1=fft(Rs,n,2); % 按行FFT
% S1=fft(Rs,n,1); % 按行FFT
figure;imagesc(abs(S1));xlabel('Lag');ylabel('amplitude');title('a');

%% AR
alfa=13.2;
AR_S1=S1(mod(round(tan(alfa*(0:n-1))),n)+1,:);
figure;imagesc(abs(AR_S1));xlabel('Lag');ylabel('amplitude');title('a');

%% 相位补偿
com=AR_S1.*exp(1j*pi*(round(tan(alfa*(0:n-1)))));
DIRFAT=fft(com);
figure;imagesc(abs(DIRFAT));xlabel('Lag');ylabel('amplitude');title('a');

%%

%% 展示线性调频信号的时域和频域
% derta_f=fs/n;
% f=(0:n-1)*derta_f;       % 采样的频率点
% figure;
% subplot(211);plot(f,abs(x));xlabel('Sampling point in the time domain');ylabel('amplitude');title('a');
% subplot(212);plot(f,abs(X1));xlabel('Sampling point in the frequency domain');ylabel('amplitude');title('b');
% SFFT估计频域
% [out_Large,~,~,~,~]=sfft(x,k);
%% 展示sfft估计的线性调频信号的频域
%out_Large=fftshift( out_Large ); X2=out_Large/fs;
% figure;
% subplot(211);plot(abs(fft(x)));xlabel('Sampling point in the frequency domain');ylabel('amplitude');
% subplot(212);plot(abs(out_Large));  xlabel('Estimate point in the frequency domain');ylabel('amplitude');
%% 计算SFFT的误差
% [ERROR1,ERROR2]= run_error(x, n, k,LARGE_FREQ, out_Large );
% disp(ERROR1);
% disp(ERROR2);
%% NUDFT
% f=zeros(n,1);
% for i=1:n
%     f(i)=Nonuniform_sampling_point(i-1,n);
% end
% True_freq= ndftld(x, n, f).';
%  figure;
%  subplot(2,1,1);plot(abs(x));title('(a)');ylabel('amplitude');xlabel('Sampling point in time domion');
%  subplot(2,1,2);plot(f,abs(True_freq));title('(b)');ylabel('amplitude');xlabel('Sampling point in frequency domion');
% figure;plot(f,abs(ndftld(x, n, f)));title('Original Signal in frequency domion');
%% NUSFFT
% [out_Large2,BB_loc,BB_est,b_loc,b_est]=nusfft(x,k,f);
%% 展示NUsfft估计的线性调频信号的频域
% figure;
% subplot(311);plot(abs(x));ylabel('amplitude');xlabel('Sampling point in time domion');title('(a) Signal in time domain');
% subplot(312);plot(f,abs(True_freq)); xlabel('Ununiform sampling point in the frequency domain');ylabel('amplitude');title('(b) Direct method');
% subplot(313);plot(f,abs(out_Large2));xlabel('Ununiform Estimate point in the frequency domain');ylabel('amplitude');title('(c) NUSFT method');
%% 相对误差
% NUsfft_err=0;
% for i=1:n
%     if out_Large2(i) ~= 0
%        NUsfft_err = NUsfft_err+(abs(True_freq(i))-abs(out_Large2(i)))/abs(True_freq(i));
%     end
% end
% NUsfft_err=NUsfft_err/k;
% disp(['L1 Error is ',  num2str(NUsfft_err)]);