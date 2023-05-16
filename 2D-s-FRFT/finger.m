clear;
ts=0.01;
fs=1/ts;
F=imread('image\fingerprint.bmp'); % 256*256
%F=rgb2gray(F);
%F=imcrop(F,[1,1,255,255]); % 251*251 剪切图像
F=double(F);
% F(210:220,150:200) = 255;% 删除部分数据
% F(25:30,120:170) = 255;
% F(155:230,238:242) = 255;
% F(75:145,86:90) = 255;

%  F=imnoise( F , 'salt & pepper' ); % 添加盐和胡椒噪声
%  F=imnoise( F , 'gaussian' ); % 添加高斯噪声

[n1,n2]=size(F);
figure;colormap("gray");imagesc((0:n1-1),(0:n2-1),real(F));  %二维图，颜色深浅表述幅度大小
set(gca,'xtick',[],'xticklabel',[]);%隐藏x轴的刻度与数字
set(gca,'ytick',[],'yticklabel',[]);%隐藏y轴的刻度与数字

%% 2个1D DFRFT
p1=1;
p2=1;
X_frft = DFRFT_2D_fft(F,p1,p2,ts,ts);X_frft=fftshift(X_frft);
figure;colormap("colorcube");mesh((0:n1-1)*fs/n1,(0:n2-1)*fs/n2,abs(X_frft));zlabel('Amplitude');axis([0 fs 0 fs]);%title("Decomposition method");
% figure;colormap("colorcube");imagesc((0:n1-1)*fs/n1,(0:n2-1)*fs/n2,abs(X_frft));  %二维图，颜色深浅表述幅度大小
%%
%[X_sfrft,us1,us2] = DFRFT_2D_sft(F,p1,p2,n1,n2,ts,ts,1);
X_sfrft=zeros(n1,n2);
[Omega,hA] =MARS_SFT_vote(F, ones(n1,1)*ones(n2,1).', n1, n2, 1, 5e+3, 1*1e5, 5, 2); %epsolon 无噪声5e+3，高斯噪声5e+1，盐和胡椒噪声9e+1
[Omega, ind] = sortrows(Omega); % 估计的位置（升序）
hA = hA(ind);   % 估计的频率值
X=zeros(n1,n2);
u = Omega(:,1);%（行位置）
v = Omega(:,2);%（列位置）
X_sfrft(sub2ind([n1,n2],u+1,v+1)) = hA/sqrt(n1*n2); X_sfrft=fftshift(X_sfrft);
figure;colormap("colorcube");mesh((0:n1-1)*fs/n1,(0:n2-1)*fs/n2,abs(X_sfrft));zlabel('Amplitude');axis([0 fs 0 fs]);
% figure;colormap("colorcube");imagesc((0:n1-1)*fs/n1,(0:n2-1)*fs/n2,abs(X_sfrft));  %二维图，颜色深浅表述幅度大小

%title(["p1=", num2str(p1)," p2=",  num2str(p2) ]);
% clear X_frft
