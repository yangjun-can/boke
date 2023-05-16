clear;
clc;
close all;
n=512;    % 信号长度
k=4;      % 信号稀疏度
%% 构造时域上的线性调频信号
fs=n;     % 采样频率
ts=1/fs;  % 采样间隔
t1=(0:n-1)*ts; % 时间向量  共1s
% t1=(-n/2+1:n/2)*ts; % 时间向量  共1s
tao=(-n/2+1:n/2)*ts;
fre2=-n/2:n/2-fs/n;
fd=[100.01, 149.99, 200, 200];   % 初始频率
mu=[ -46,  55,  108,  -125];   % 调频率
s=zeros(1,n);
for ii=1:k
    x=exp(1j*pi*mu(ii)*(t1.^2)).*exp(1j*2*pi*fd(ii)*t1);
    s=x+s;      % 线性调频信号（时域）
end
figure;
subplot(211);plot(real(s));xlabel('Sampling point in the time domain');ylabel('amplitude');title('a');
subplot(212);plot(abs(fft(s)));xlabel('Sampling point in the frequency domain');ylabel('amplitude');title('b');

for ch=1:k
    alpha=atan(2*ts*mu(ch));
    %% 瞬时自相关函数
    % Rs = zeros(n,n);
    % for i = 1:k
    %     %     s1=exp(1j*2*pi*fd(i)*tao);
    %     s1=exp(1j*4*pi*fd(i)*tao);
    %     s1_all=repmat(s1, n, 1);
    %     %     s2=exp(1j*1*pi*mu(i)*t1.'*tao);
    %     s2=exp(1j*4*pi*mu(i)*t1.'*tao);
    %     rs=s1_all.*s2;
    %     Rs=Rs+rs;
    % end
    % %% DAF
    % S1=fftshift(fft(fftshift(Rs))); % 按列FFT
    % figure();
    % imagesc(tao,fre2,abs(S1));xlabel('Lag ');ylabel('Frequency (Hz)');set(gca,'ydir','normal');
    % %% AR(RDAF)
    % for i=0:n-1
    %     AR_S1(:,i+1)=S1(mod((0:n-1)+round(tan(alpha)*i),n)+1,i+1);
    % end
    % figure;imagesc(tao,fre2,abs(AR_S1));xlabel('Lag');ylabel('Frequency (Hz)');set(gca,'ydir','normal');title('a');

    %% RDAF(定义)
    AR_S1 = zeros(n,n);
    for kt=-n/2+fs/n:n/2
        for ntao=-n/2+fs/n:n/2
            for i=1:k
                RDAF=exp(1j*2*pi*(2*fd(i)+mu(i))*ntao*ts)*exp(-1j*pi*round(tan(alpha)*ntao))*exp(-1j*pi*kt)*sin(2*pi*(mu(i)*ntao*ts-(kt+round(tan(alpha)*ntao))/2))/sin(2*pi*ts*(mu(i)*ntao*ts-(kt+round(tan(alpha)*ntao))/2));
                AR_S1(kt+n/2,ntao+n/2)=AR_S1(kt+n/2,ntao+n/2)+RDAF;
            end
        end
    end
    figure;imagesc(tao,0:n-1,abs(AR_S1));
    xlabel('Lag(s)','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);
    title(['$\alpha  = $', num2str(alpha)],'interpreter','latex','FontSize',20);
    set(gca,'ydir','normal');

    %% 相位补偿
    AR=max(AR_S1(n/2-1:n/2+1,:));
    com=AR.*exp(1j*pi*(round(tan(alpha)*(-n/2+fs/n:n/2))));
    %% 相干积分
    %     DI=fft(fftshift(com)); % DI=fft(com);
    %     figure;plot(0:n-1,abs(DI));xlabel('Lag Frequency(Hz)');ylabel('Amplitude');set(gca,'ydir','normal');title('a');
    f=zeros(n,1);
    for index=1:n
        f(index)=Nonuniform_sampling_point(index-1,n);
    end
    True_freq = ndftld(com, n, f).';
    figure; plot(f,abs(True_freq));
    xlabel('Lag Frequency(Hz)','FontSize',16);ylabel('Amplitude','FontSize',16);
    set(gca,'ydir','normal');
    title(['$\alpha  = $', num2str(alpha)],'interpreter','latex','FontSize',20);

    [DI,BB_loc,BB_est,b_loc,b_est]=nusfft(com,1,f);
    figure; plot(f,abs(DI));
    xlabel('Lag Frequency(Hz)','FontSize',16);ylabel('Amplitude','FontSize',16);
    set(gca,'ydir','normal');
    title(['$\alpha  = $', num2str(alpha)],'interpreter','latex','FontSize',20);
end