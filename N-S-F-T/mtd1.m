clc;clear;
N=2^8;
fs=1e3;
derta_f=fs/N;
f=(0:N-1)*derta_f;%���ȹ�һ��
R0=9e3;%��ʱ���״���Ŀ��ľ���
lamda=0.04;%�״ﲨ��
c0=3e8;
prf=1100;%�����ظ�Ƶ��
v=140;%�״��ػ�ƽ̨�ٶ�
D=4;%���߿׾�����
Lsar=lamda*R0/D;%�ϳɿ׾�����
Tsar=Lsar/v;%�ϳɿ׾�ʱ��
ts=1/fs;
t1=(-N/2+1:N/2)*ts;
t=t1.*(0<t1&t1<Tsar);

TargetNumber=3;%��R0���뵥Ԫ�ڵ��˶�Ŀ�����

v_r(1:TargetNumber)=[1 2.5 5];%Ŀ�꾶���ٶ�m/s
v_c(1:TargetNumber)=[5.5 -4 6];%Ŀ�������ٶ�
a_r(1:TargetNumber)=[3.97 5.2 2.2];%Ŀ�꾶����ٶ�
snr(1:TargetNumber)=[4 -4 6];
dalta(1:TargetNumber)=[2 0.8 2];%Ŀ��ɢ��ϵ��

fd=2*v_r/lamda;%����������Ƶ��
K=2*(-(v-v_c).^2+R0*a_r)/(lamda*R0);%�����յ�Ƶ��

s=zeros(1,N);

for i=1:TargetNumber
    s2=dalta(i)*(exp(1i*2*pi*fd(i)*(Tsar+t)).*exp(1i*pi*K(i)*(Tsar+t).^2));
    s2=awgn(s2,snr(i),'measured');
    s=s+s2;%�ز��ź�
end
k=3;
a=0.9*K./fd;
b=9.3./fd;
c=3.3;
d=(1+b*c)./a;
h=dalta*((exp(-1i*a/b*ones(3,1)*t).^2-1i*a/(2*b)*(2*R0/c0).^2*ones(3,N)+1i*a/b*ones(3,1)*t).*(2*R0/c0).*(ones(3,1)*conj(s)));%��������ƥ���˲���
[H,us] = LCT(s,a(2),b(2),c,d(2), N, ts);
[X_lct1,us1] = LCT(s,a(2),b(2),c,d(2), N, ts);
Y=X_lct1.*conj(H);
y=ILCT(Y,a(2),b(2),c,d(2), N, us1);
figure(1);plot(f,s);xlabel('the sampling point in time domain');ylabel('amplitude');
fnu=Nonuniform_sampling_point(0:N-1,N).';
[out_Large,BB_loc,BB_est,b_loc,b_est]=nusfft(y,3,fnu);
 figure(2);
 plot(fnu,fftshift(abs(out_Large)));xlabel('the sampling point in linear canonical domain');ylabel('amplitude');%axis([0,1000,0,350]);

% Z=LCT(y,a(2),b(2),c,d(2), N, ts);
% figure(2);
% plot(f,fftshift(abs(Z)));xlabel('the sampling point in linear canonical domain');ylabel('amplitude');%axis([0,1000,0,350]);

% Z1=LCT_sfft(y, a(2),b(2),c,d(2), N, ts,k);
%figure(3);
% plot(f,fftshift(abs(Z1)));xlabel('the sampling point in sparse linear canonical domain');ylabel('amplitude');

% a_best1=acot(K(1))*2/pi;       % �ԡ�Ϊ��λ
% a_best2=acot(K(2))*2/pi;       % �ԡ�Ϊ��λ
% a_best3=acot(K(3))*2/pi;       % �ԡ�Ϊ��λ
% p=30*a_best2;
% [X_frft,us2] = dfrft_sfft(y,p, N, ts, k);
% figure(4);
% plot(f,fftshift(abs(X_frft)));xlabel('the sampling point in sparse fractional fourier domain');ylabel('amplitude');axis([0,1000,0,120]);

% X=sfft(y,3);
% figure(6);plot(f,fftshift(abs(X)));
% xlabel('the sampling point in sparse fourier domain');ylabel('amplitude');axis([0,1000,0,120]);
