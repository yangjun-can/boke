function [out] = DFRFT_2D_fft(matrix,p1,p2,ts1,ts2)
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
     temp(i,:) = Dfrft(matrix(i,:),p1,ts1); %Pei
end
for i = 1:N2
     out(:,i) = Dfrft(temp(:,i),p2,ts2); %pei
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
% K=zeros(M,N);
% for m=0:M-1
%     for n=0:N-1
% %         K(m+1,n+1)=exp(-1i*csc(alpha)*m*n*ts*us)*exp(1i/2*cot(alpha)*(ts*n)^2);
%         K(m+1,n+1)=exp(-1i*S1*2*pi*m*n/M)*exp(1i/2*cot(alpha)*(ts*n)^2);
%     end
% %     chip_u=sqrt((1-1i*cot(alpha))/2/pi)*ts*exp(1i/2*cot(alpha)*(us*m)^2);
%      chip_u=sqrt(S1*(sin(alpha)-1i*cos(alpha))/M)*exp(1i/2*cot(alpha)*(us*m)^2);
%     F(m+1)=chip_u*K(m+1,:)*f;
% end

%% 分解为fft
t=ts*(0:N-1);                 % 时域采样点
chirp_t = exp(1i/2*cot(alpha)*t.^2); %
x_t=f.'.*chirp_t;                    % 第一次调制
if sin(alpha)>0
    temp=fft(x_t);
else
    temp=ifft(x_t)*N;
end
u=us*(0:M-1);
chirp_u = exp(1i/2*cot(alpha)*u.^2);
k1=sqrt(S1*(sin(alpha)-1i*cos(alpha))/M);
F=k1*chirp_u.*temp;     % 第二次调制

end

%% % 基于特征值分解的FRFT
function [y ,E]= Disfrft(f,a,p)
% Computes discrete fractional Fourier transform
% of order a of vector f
% p (optional) is order of approximation, default N/2

N = length(f);  % 信号长度
even = ~rem(N,2); % rem取余数   even偶数为1，奇数为0
shft = rem((0:N-1) + fix(N/2),N)+1; % fix向零方向取整
f = f(:);
if (nargin == 2), p = N/2; end
p = min(max(2,p),N-1);
E = dFRFT(N,p);
y(shft,1) = E*(exp(-1i*pi/2*a*([0:N-2 N-1+even])).' .*(E'*f(shft)));
end
function E = dFRFT(N,p)
%
% function E = dFRFT(N,a,p) returns the NxN eigenvectors of the
% Fourier transform matrix
% The optional argument p is the order of approximation

global E_saved p_saved

if (length(E_saved) ~= N | p_saved ~= p),
    E = make_E(N,p);
    E_saved = E; p_saved = p;
else
    E = E_saved;
end;
end
function E = make_E(N,p)

% Returns sorted eigenvectors and eigenvalues of corresponding vectors

% Construct matrix H, use approx order ord

d2 = [1 -2 1]; d_p = 1; s = 0; st = zeros(1,N);
for k = 1:p/2,
    d_p = conv(d2,d_p);
    st([N-k+1:N,1:k+1]) = d_p; st(1) = 0;
    temp = [1:k;1:k]; temp = temp(:)'./[1:2*k];
    s = s + (-1)^(k-1)*prod(temp)*2*st;
end;

% H = circulant + diagonal

col = (0:N-1)'; row = (N:-1:1);
idx = col(:,ones(N,1)) + row(ones(N,1),:);
st = [s(N:-1:2).';s(:)];
H = st(idx) + diag(real(fft(s)));

% Construct transformation matrix V

r = floor(N/2);
even = ~rem(N,2);
V1 = (eye(N-1) + flipud(eye(N-1))) / sqrt(2);
V1(N-r:end,N-r:end) = -V1(N-r:end,N-r:end);
if (even), V1(r,r) = 1; end
V = eye(N); V(2:N,2:N) = V1;

% Compute eigenvectors

VHV = V*H*V';
E = zeros(N);
Ev = VHV(1:r+1,1:r+1);           Od = VHV(r+2:N,r+2:N);
[ve,ee] = eig(Ev);               [vo,eo] = eig(Od);

%
% malab eig returns sorted eigenvalues
% if different routine gives unsorted eigvals, then sort first
%
% [d,inde] = sort(diag(ee));      [d,indo] = sort(diag(eo));
% ve = ve(:,inde');               vo = vo(:,indo');
%

E(1:r+1,1:r+1) = fliplr(ve);     E(r+2:N,r+2:N) = fliplr(vo);
E = V*E;

% shuffle eigenvectors

ind = [1:r+1;r+2:2*r+2]; ind = ind(:);
if (even), ind([N,N+2]) = []; else ind(N+1) = []; end
E = E(:,ind');
end

%%  基于采样的FRFT(Haldun M. Ozaktas, Orhan, Digital Computation of the Fractional Fourier Transform)

function Faf = frft(f, a)
% 将分数傅里叶变换分解为一个调频卷积和两个调频乘积
%分数阶傅里叶变换函数
%输入参数f为原始信号，a为阶数
%输出结果为原始信号的a阶傅里叶变换
f = f(:);
N = length(f);%总采样点数
shft = rem((0:N-1)+fix(N/2),N)+1; % 实际信号坐标[-N/2+1，N/2]
sN = sqrt(N);%看原文中对离散傅里叶变换的定义，有这个乘积项
a = mod(a,4);%按周期性将a定义在[0,4]

%特殊情况直接处理
if (a==0), Faf = f; return; end%自身
if (a==2), Faf = flipud(f); return; end%f(-x)  flipud 矩阵上下翻转
if (a==1), Faf(shft,1) = fft(f(shft))/sN; return; end%f的傅里叶变换
if (a==3), Faf(shft,1) = ifft(f(shft))*sN; return; end%f的逆傅里叶变换

%利用叠加性将阶数变换到0.5 < a < 1.5
if (a>2.0), a = a-2; f = flipud(f); end%a=2是反转
if (a>1.5), a = a-1; f(shft,1) = fft(f(shft))/sN; end%a=1是傅里叶变换
if (a<0.5), a = a+1; f(shft,1) = ifft(f(shft))*sN; end%a=-1是逆傅里叶变换

%开始正式的变换
alpha = a*pi/2;
tana2 = tan(alpha/2);
sina = sin(alpha);

f = [zeros(N-1,1) ; interp(f) ; zeros(N-1,1)];%使用香农插值，拓展为4N
% 以下操作对应原论文中公式（29）
% 线性调频预调制
chrp = exp(-1i*pi/N*tana2/4*(-2*N+2:2*N-2)'.^2);
f = chrp.*f;
% 线性调频卷积
c = pi/N/sina/4;
Faf = fconv(exp(1i*c*(-(4*N-4):4*N-4)'.^2),f);
Faf = Faf(4*N-3:8*N-7)*sqrt(c/pi);
% 线性调频后调制
Faf = chrp.*Faf;
% 乘以最前面的A_Phi项
Faf = exp(-1i*(1-a)*pi/4)*Faf(N:2:end-N+1);
end

function xint=interp(x)%香农插值
% sinc interpolation
N = length(x);
y = zeros(2*N-1,1);
y(1:2:2*N-1) = x;
xint = fconv(y(1:2*N-1), sinc([-(2*N-3):(2*N-3)]'/2));%计算卷积
xint = xint(2*N-2:end-2*N+3);
end

function z = fconv(x,y)%利用fft快速计算卷积
N = length([x(:);y(:)])-1;%计算最大点数
P = 2^nextpow2(N);%补零
z = ifft( fft(x,P) .* fft(y,P));%频域相乘，时域卷积
z = z(1:N);%去零
end

