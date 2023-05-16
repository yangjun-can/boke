
N=512; % 信号长度
t=0:N-1; % 采样点
k=6; % 信号稀疏度
%构造信号
x1=2.5*sin(51.05*pi*2*t/N);
x2=3.0*sin(52.45*pi*2*t/N);
x3=1.9*sin(215*pi*2*t/N);
Sig = x1 + x2 + x3; % 信号（时域）
figure;plot(t,abs(Sig));ylabel('amplitude');xlabel('Sampling point in time domion');
%DFT
freq1=fft(Sig)/N;
ii=find(t>=N/2);
freq1(ii)=NaN;
figure;
plot(t,abs(freq1));ylabel('amplitude');xlabel('Sampling point in frequency domion');
% axis off; title('(a)');

%非均匀采样点
% for i=0:N-1
%     if 0<=i && i<=N/2
%         fs(i+1)=0.5*i;
%     else
%         if i>N/2 && i<3*N/4
%             fs(i+1)=i-N/4;
%         else
%              fs(i+1)=2*i+1-N;
%         end
%     end
% end

% for i=0:N-1
%     if 0<=i && i<=N/2
%         fs(i+1)=0.125*i;
%     else
%         if i>N/2 && i<15*N/16
%             fs(i+1)=i-7*N/16;
%         else
%             fs(i+1)=8*i+7-7*N;
%         end
%     end
% end

%NSFT
unf=(0:N-1)+0.05;
[out_Large2,~,~,~,~]=nusfft(Sig,k,unf.');
out_Large2=out_Large2/N;
iun=find(unf>=N/2);
out_Large2(iun)=NaN;
figure;
plot(unf,abs(out_Large2));title('(b)');ylabel('amplitude');xlabel('Sampling point in frequency domion');
h1=axes('position',[0.32 0.4 0.1 0.5]);%局部图的位置
axis(h1);
plot(h1, unf(51:55),abs(out_Large2(51:55)));

figure;
subplot(121);plot(abs(outlarge));xlabel('Sampling point in the frequency domain');ylabel('amplitude');title('(a) SFT method');
subplot(122);plot(abs(out_Large2));xlabel('Sampling point in the frequency domain');ylabel('amplitude');title('(b) NUSFT method');
% 
% NDFT
% fs=zeros(N,1);
% for i=0:N-1
%     if 0<=i && i<=N/2
%         fs(i+1)=0.25*i;
%     else
%         if i>N/2 && i<7*N/8
%             fs(i+1)=i-3*N/8;
%         else
%             fs(i+1)=4*i+3-3*N;
%         end
%     end
% end
% 
% freq2=ndftld(Sig, N, fs)/N;
% i=find(fs>=N/2);
% freq2(i)=NaN;
%  figure;
% plot(fs,abs(freq2));title('(b)');ylabel('amplitude');xlabel('Sampling point in frequency domion');
% axis off
% 局部放大
% h1=axes('position',[0.32 0.4 0.1 0.5]);%局部图的位置
% axis(h1);
% plot(h1, fs(203:214),abs(freq2(203:214)));
