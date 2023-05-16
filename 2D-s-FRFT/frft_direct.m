function [out] = frft_direct(matrix,p1,p2,ts1,ts2)
%
% computes 2-D FRFT of given matrix with given angles
%
% IN : matrix: matrix to be transformed
%      p1 p2 : order of FrFT in x and y direction
%      ts1 ts2 : 时域采样间隔

[N1,N2]=size(matrix);
temp=zeros(N1,N2);
out=zeros(N1,N2);
for i = 1:N1
    temp(i,:) = Dfrft(matrix(i,:),p1,ts1);
end
for i = 1:N2
    out(:,i) = Dfrft(temp(:,i),p2,ts2);
end
end

%% 采样型DFRFT(Soo_Chang Pei, Jian-Jiun Ding Closed-Form Discrete Fractional and Affine Fourier Transforms)
function F = Dfrft(f, p, ts)
% f 信号
% p 分数阶数
% ts 连续信号的时域采样间隔
% F DFRFT
f=f(:);
N=size(f,1);
M=N;
F=zeros(M,1);
p = mod(p,4);                   %按周期性将p定义在[0,3]
alpha=p*pi/2;                   % 旋转角度（旋转角以2pi为周期）
S1=sign(sin(alpha));
us=2*pi*S1*sin(alpha)/(ts*M);   % 分数域采样间隔
% shft = rem((0:N-1)+fix(N/2),N)+1; %此项等同于fftshift(1:N)，起到翻转坐标轴的作用
sN = sqrt(N);                   % 看原文中对离散傅里叶变换的定义，有这个乘积项
%% 特殊情况直接处理
if (p==0), F = f; return; end %自身
if (p==2), F=f(mod(-(0:N-1),N)+1); return; end%f(-x)
if (p==1), F = fft(f)/sN; return; end%f的傅里叶变换
if (p==3), F = ifft(f)*sN; return; end%f的逆傅里叶变换
%% 采样型定义
K=zeros(M,N);
for m=0:M-1
    for n=0:N-1
%         K(m+1,n+1)=exp(-1i*csc(alpha)*m*n*ts*us)*exp(1i/2*cot(alpha)*(ts*n)^2);
        K(m+1,n+1)=exp(-1i*S1*2*pi*m*n/M)*exp(1i/2*cot(alpha)*(ts*n)^2);
    end
%     chip_u=sqrt((1-1i*cot(alpha))/2/pi)*ts*exp(1i/2*cot(alpha)*(us*m)^2);
     chip_u=sqrt(S1*(sin(alpha)-1i*cos(alpha))/M)*exp(1i/2*cot(alpha)*(us*m)^2);
    F(m+1)=chip_u*K(m+1,:)*f;
end
end