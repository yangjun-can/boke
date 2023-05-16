%%
% 相位编码信号通过非线性相位调制，其相位调制函数是离散的有限状态，属于离散
% 编码脉冲压缩信号。这种信号的突出优点是采用脉冲压缩技术后，雷达的峰值发射
% 功率得到显著降低，从而实现低截获的目的。但当回波信号与匹配滤波器有多普勒
% 失谐时，滤波器起不了脉冲压缩的作用，所以有时称之为多普勒灵敏信号，常用于
% 目标多普勒变化范围较窄的场合。
%%
%
% clc
% clear
% close all
% T = 1;%采样时间
% f = 1e5;%采样率
% t = 0:1/f:(T-1/f); % 采样点
% n = length(t); % 采样点数
% % 二进制相移键控，0到1,1到0，绝对相位控制
% binary_code = '1111100110101';
% len_code = length(binary_code);
% % 一个码元对应的脉冲长度
% f_bpsk = 1000; % 频率1000Hz，周期1ms
% T_bpsk = 1/f_bpsk; %单个码元持续时间
% A_bpsk = 10; % 幅值
% t_mayuan = 0:1/f:(T_bpsk-1/f); % 一个码元所持续的时间内的采样点时刻
% y_mayuan = A_bpsk*sin(2*pi*f_bpsk*t_mayuan); % 一个码元的波形，用于表示0
% t_bpsk = 0:1/f:(T_bpsk*len_code-1/f);
% y_bpsk = [];
% for i=1:len_code
%     tmp = str2num(binary_code(i));
%     if tmp == 0
%         y_bpsk = [y_bpsk,y_mayuan];
%     end
%     if tmp == 1
%         y_bpsk = [y_bpsk,-y_mayuan];
%     end
% end
% figure(2)
% plot(t_bpsk,y_bpsk)
% grid on;
% xlabel('t/s');
% title('‘1111100110101’的BPSK编码信号');
% axis([0 1e-2 -2*A_bpsk 2*A_bpsk]);
%
% figure
% subplot(4,1,1);      %第一个图
% t = 0:0.01:8;
% d = [0 0;0.5 0;1 1;1.5 1;2 0;2.5 0;3 0;3.5 0;4 1;4.5 1;5 0;5.5 0;6 1;6.5 1;7 1;7.5 1;8 1];
% signal = pulstran(t - 0.25, d, 'rectpuls', 0.5);   %矩形脉冲
% plot(t, signal);
% axis([0 8 -0.5 1.5]);
% title('输入的二进制序列');
% grid();
% 
% subplot(4,1,2);
% d1 = [0 1;0.5 1;1 -1;1.5 -1;2 1;2.5 1;3 1;3.5 1;4 -1;4.5 -1;5 1;5.5 1;6 -1;6.5 -1;7 -1;7.5 -1;8 -1];
% I = pulstran(t-0.25, d1, 'rectpuls', 0.5);
% plot(t, I);
% axis([0 8 -2 2]);
% title('映射过来的I路信号');
% grid();
% 
% subplot(4,1,3);
% d3 = [0 0;0.5 0;1 0;1.5 0;2 0;2.5 0;3 0;3.5 0;4 0;4.5 0;5 0;5.5 0;6 0;6.5 0;7 0;7.5 0;8 0];
% Q = pulstran(t-0.25, d3, 'rectpuls', 0.5);
% plot(t, Q);
% axis([0 8 -2 2]);
% title('映射过来的Q路信号');
% grid();
% 
% s1 = I.*cos(2*pi*5*t);
% s2 = Q.*sin(2*pi*5*t);
% output = s1 - s2;
% plot(t, output);
% axis([0 8 -2 2]);
% title('调制之后的信号');
% grid();
% subplot(4,1,4);plot(fft(output));
%%  产生一个13位巴克码编码的二相码示例  %%
% clear
% clc
% % code = [1 1 1 -1];%13位巴克码
% code = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];%13位巴克码
% tao = 0.5e-6;  %chip时宽
% fc = 20e6;  %载频(信号频率)
% % fs = 200e6;  %采样率
% fs = 20e6;  %采样率
% t_tao = 0:1/fs:tao-1/fs;  %chip时宽内采样时间点序列
% n = length(code);  %码长
% % phase = 0;  %没必要在这定义相位
% t = 0:1/fs:13*tao-1/fs;
% s = zeros(1,length(t));
% for ii = 1:n
%     if code(ii) == 1
%         phase = 0;
%     else
%         phase = pi;
%     end
%     s(1,(ii-1)*length(t_tao)+1:ii*length(t_tao)) = cos(2*pi*fc*t_tao+phase);
% end
% figure
% plot(t,s,'b-o')
% % plot(t,s)
% xlabel('t(单位：秒)');
% title('二相码（13位巴克码）');
% %  自相关函数  %%
% [a,b] = xcorr(code);
% % d = abs(a);
% figure
% plot(b,a,'r-o')
% axis tight


%% clear;
N=512; % 信号长度
t=0:N-1; % 采样点
k=3; % 信号稀疏度
%构造信号
x1=exp(1i*27.05*pi*2*t/N);
x2=exp(1i*140.05*pi*2*t/N);
x3=exp(1i*371.05*pi*2*t/N);

S_bpsk=x1+x2+x3;
%% BPSK （二相编码信号）
% g=[0 1 0 0 1 1 1 0];  %8位巴克码barker
g=[1 1 1 0];  %4位巴克码barker

fs=N;
n=N/length(g);
v=[];
M=length(g);
N=length(t);
for i=1:length(g)
    if g(i)==1
        w=ones(1,n);
    elseif g(i)==0
        w=-ones(1,n);
    end
    v=[v w];
end
% figure; plot(v);
S_bpsk1=S_bpsk.*v;
%% 构造噪声
% N_noise=(randn(1,length(t))-1i*randn(1,length(t))).*exp(1i*2*pi*fo*t);  %
% % ratio=-3;  %
% ratio=10;
% Ex1=mean(abs(S_bpsk1).^2);
% Ex2=mean(abs(N_noise).^2);
% N_A=sqrt(Ex1/(Ex2*10^(ratio/10)));
% N_noise=N_A*(randn(1,length(t))-1i*randn(1,length(t))).*(exp(1i*2*pi*fo*t)+exp(1i*2*pi*fo1*t)+exp(1i*2*pi*fo2*t));  %
%%
% S_bpsk=S_bpsk1+N_noise;
% S_bpsk=S_bpsk1;
unf=t+0.05;
freq=ndftld(S_bpsk, N, unf.');
% freq=chebfun.nufft(S_bpsk.',unf.'/N,2);
[NUSFFT_freq,~,~,~,~]=nusfft(S_bpsk,k,unf.');

% freq1=chebfun.nufft(S_bpsk1.',unf.'/N,2);
freq1=ndftld(S_bpsk1, N, unf.');
[NUSFFT_freq1,~,~,~,~]=nusfft(S_bpsk1,k,unf.');

% 画图
figure;
subplot(231)
plot(t,S_bpsk) ;xlabel('Sampling points in time domain'); ylabel('Amplitude'); title('(a)')
subplot(232)
plot(unf,abs(freq));xlabel('Sampling points in frequency domain'); ylabel('Amplitude'); title('(b)')
subplot(233)
plot(unf,abs(NUSFFT_freq));xlabel('Sampling points in frequency domain'); ylabel('Amplitude'); title('(c)')
% Braker调制后的BPSK
subplot(234)
plot(t,S_bpsk1) ;xlabel('Sampling points in time domain'); ylabel('Amplitude'); title('(a)')
subplot(235)
plot(unf,abs(freq1)) ;xlabel('Sampling points in time domain'); ylabel('Amplitude'); title('(d)')
subplot(236)
plot(unf,abs(NUSFFT_freq1));xlabel('Sampling points in frequency domain'); ylabel('Amplitude'); title('(f)')