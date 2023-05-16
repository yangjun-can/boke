% 使用2D SFRFT实现图像融合
% clc
clear
p1=0.88;
p2=0.88;
mul=0.002;
ts1=1e-5;
ts2=1e-5;
sigma_n=1;

A = imread('image\p3.png');  % 读取图像
A=double(A);
% figure;imshow(A/max(max(A)));
B = imread('image\p4.png'); % 读取图像
B=double(B);
% figure;imshow(B/max(max(B)));
[N1,N2]=size(A);
ind1=A>=B;
% num1=sum(sum(ind1~=0));
% AB=num1/N1/N2*A+(1-num1/N1/N2)*B;
AB=zeros(N1,N2);
AB(ind1)=A(ind1);
AB(~ind1)=B(~ind1);
% figure;imshow(AB/max(max(AB)));

% [A_frft,~,~] = DFRFT_2D_fft2(A,p1,p2,N1,N2,ts1,ts2);
% figure;mesh(0:N1-1,0:N2-1,abs(A_frft));
[A_frft,us1,us2] = DFRFT_2D_sft(A,p1,p2,N1,N2,ts1,ts2,sigma_n);
E_A=A-A_frft;
% [B_frft,us1,us2] = DFRFT_2D_fft2(B,p1,p2,N1,N2,ts1,ts2);
[B_frft,us1,us2] = DFRFT_2D_sft(B,p1,p2,N1,N2,ts1,ts2,sigma_n);
% figure;mesh(0:N1-1,0:N2-1,abs(B_frft));
E_B=B-B_frft;

ind2=E_A>=E_B;
% num2=sum(sum(ind2~=0));
% E_AB=num2/N1/N2*E_A+(1-num2/N1/N2)*E_B;
E_AB=zeros(N1,N2);
E_AB(ind2)=E_A(ind2);
E_AB(~ind2)=E_B(~ind2);
% figure;imshow(abs(E_AB/max(max(E_AB))));
Fusion=AB+mul*E_AB;
figure;imshow(abs(Fusion/max(max(Fusion))));
outval = Hero_GGBOND(Fusion)  %图像平均梯度
SFRUENCY=sfrquency(Fusion)    %图像空间频率

%% 图像平均梯度 —— 衡量图像清晰度、表达细节的能力
function outval = Hero_GGBOND(img) 
% OUTVAL = AVG_GRADIENT(IMG) 
if nargin == 1 
    img = double(img); 
    % Get the size of img 
    [r,c,b] = size(img); 
    dx = 1; 
    dy = 1; 
    for k = 1 : b 
        band = img(:,:,k); 
        [dzdx,dzdy] = gradient(band,dx,dy); 
        s = sqrt((dzdx .^ 2 + dzdy .^2) / 2); 
        g(k) = sum(sum(s)) / ((r - 1) * (c - 1)); 
    end 
    outval = mean(g); 
else 
    error('Wrong number of input!'); 
end
end

%% 图像空间频率 —— 衡量图像质量、细节丰富度
function SFRUENCY=sfrquency(img)
[m,n]=size(img);
rf=0;
cf=0;
for i=1:m
    for j=2:n
        rf=rf+(img(i,j)-img(i,j-1))^2;
    end
end
RF=(rf/(m*n))^(1/2);
for i=2:m
    for j=1:n
        cf=cf+(img(i,j)-img(i-1,j))^2;
    end
end
CF=(cf/(m*n))^(1/2);
SF=(RF^2+CF^2)^(1/2);
SFRUENCY=SF;
end