 %主函数main

%close all;

tic

clear;

thin1=fingeryuchuli('image\54_4.tif');
thin2=fingeryuchuli('image\54_5.tif');
figure;
txy1=point(thin1);
txy2=point(thin2);
% [w1,txy1]=guanghua(thin1,txy1);
% [w2,txy2]=guanghua(thin2,txy2);
H=fspecial(gaussian);
w1=imfliter(thin1,H);
w2=imfliter(thin2,H);
thin1=w1;
thin2=w2;
txy1=cut(thin1,txy1);
txy2=cut(thin2,txy2);
[pxy31,error2]=last1(thin1,8,txy1,60)

[pxy32,error2]=last1(thin2,8,txy2,60)

error=1;

num=20;

cxy1=pxy31;

cxy2=pxy32;

d1=distance(cxy1(1,1),cxy1(1,2),num,thin1);

d2=distance(cxy2(1,1),cxy2(1,2),num,thin2);

f=(sum(abs((d1./d2)-1)));

if ff<=1
error=0;

else

error=1;

end

c11=find_point(cxy1(1,1),cxy1(1,2),txy1,1);

c12=find_point(cxy1(1,1),cxy1(1,2),txy1,2);

c21=find_point(cxy2(1,1),cxy2(1,2),txy2,1);

c22=find_point(cxy2(1,1),cxy2(1,2),txy2,2);

cxy1(2,:)=c11;

cxy1(3,:)=c12(2,:);

cxy2(2,:)=c21;

cxy2(3,:)=c22(2,:);

x11=cxy1(1,1); y11=cxy1(1,2);

x12=cxy1(2,1); y12=cxy1(2,2);

x13=cxy1(3,1); y13=cxy1(3,2);

x21=cxy2(1,1); y21=cxy2(1,2);

x22=cxy2(2,1); y22=cxy2(2,2);

x23=cxy2(3,1); y23=cxy2(3,2);

dd1(1)=juli(x11,y11,x12,y12);

dd1(2)=juli(x12,y12,x13,y13);

dd1(3)=juli(x13,y13,x11,y11);

dd2(1)=juli(x21,y21,x22,y22);

dd2(2)=juli(x22,y22,x23,y23);

dd2(3)=juli(x23,y23,x21,y21);

ff=(sum(abs((dd1./dd2)-1)))

if ff<=1

error=0;

else

error=1;

end

cxy1(2:41,:)=find_point(pxy31(1,1),pxy31(1,2),txy1,40);

cxy2(2:41,:)=find_point(pxy32(1,1),pxy32(1,2),txy2,40);

f11=length(find(cxy1(:,3)==2));

f12=length(find(cxy1(:,3)==6));

f21=length(find(cxy2(:,3)==2));

f22=length(find(cxy2(:,3)==6));

fff=abs(f11-f21)/(f11+f12);

toc

%% 预处理
function img = fingeryuchuli(path)
M=0;var=0;
I=double(imread(path));
[m,n,p]=size(I);
for x=1:m
for y=1:n

M=M+I(x,y);

end

end

M1=M/(m*n);

for x=1:m

for y=1:n

var=var+(I(x,y)-M1).^2;

end

end

var1=var/(m*n);

for x=1:m

for y=1:n

if I(x,y)>=M1

I(x,y)=150+sqrt(2000*(I(x,y)-M1)/var1);

else

I(x,y)=150-sqrt(2000*(M1-I(x,y))/var1);

end

end

end

% figure; imshow(I(:,:,3)./max(max(I(:,:,3))));title('归一化')
figure; imshow(I./max(max(I)));title('归一化')

%************************************************************************

M =3; %3*3

H = m/M; L= n/M;

aveg1=zeros(H,L);

var1=zeros(H,L); %计算每一块的平均值

for x=1:H

for y=1:L

aveg=0;var=0;

for i=1:M

for j=1:M

aveg=I(i+(x-1)*M,j+(y-1)*M)+aveg;

end

end

aveg1(x,y)=aveg/(M*M); %计算每一块的方差

for i=1:M

for j=1:M

var=(I(i+(x-1)*M,j+(y-1)*M)-aveg1(x,y)).^2+var;

end

end

var1(x,y)=var/(M*M);

end

end

Gmean=0;Vmean=0;

for x=1:H

for y=1:L

Gmean=Gmean+aveg1(x,y);

Vmean=Vmean+var1(x,y);

end

end

Gmean1=Gmean/(H*L); %所有块的平均值

Vmean1=Vmean/(H*L); %所有块的方差

gtemp=0;gtotle=0;vtotle=0;vtemp=0;

for x=1:H

for y=1:L

if Gmean1>aveg1(x,y)

gtemp=gtemp+1;

gtotle=gtotle+aveg1(x,y);

end

if Vmean1<var1(x,y)

vtemp=vtemp+1;

vtotle=vtotle+var1(x,y);

end

end

end

G1=gtotle/gtemp;V1=vtotle/vtemp;

gtemp1=0;gtotle1=0;vtotle1=0;vtemp1=0;

for x=1:H

for y=1:L

if G1<aveg1(x,y)

gtemp1=gtemp1-1;

gtotle1=gtotle1+aveg1(x,y);

end

if 0<var1(x,y)<V1

vtemp1=vtemp1+1;

vtotle1=vtotle1+var1(x,y);

end

end

end

G2=gtotle1/gtemp1;V2=vtotle1/vtemp1;

e=zeros(H,L);

for x=1:H

for y=1:L

if aveg1(x,y)>G2 && var1(x,y)<V2

e(x,y)=1;

end

if aveg1(x,y)< G1-100 && var1(x,y)< V2

e(x,y)=1;

end

end

end

for x=2:H-1

for y=2:L-1

if e(x,y)==1

if e(x-1,y) + e(x-1,y+1) +e(x,y+1) + e(x+1,y+1) + e(x+1,y) + e(x+1,y-1) + e(x,y-1) + e(x-1,y-1) <=4

e(x,y)=0;

end

end

end

end

Icc = ones(m,n);

for x=1:H

for y=1:L

if e(x,y)==1

for i=1:M

for j=1:M

I(i+(x-1)*M,j+(y-1)*M)=G1;

Icc(i+(x-1)*M,j+(y-1)*M)=0;

end

end

end

end

end

% figure, imshow(I(:,:,3)./max(max(I(:,:,3))));title('分割');
figure, imshow(I./max(max(I)));title('分割');
%************************************************************************

temp=(1/9)*[1 1 1;1 1 1;1 1 1]; % 模板系数、均值滤波

Im=double(I);

In=zeros(m,n);

for a=2:m-1;

for b=2:n-1;

In(a,b)=Im(a-1,b-1)*temp(1,1)+Im(a-1,b)*temp(1,2)+Im(a-1,b+1)*temp(1,3)+Im(a,b-1)*temp(2,1)+Im(a,b)*temp(2,2)+Im(a,b+1)*temp(2,3)+Im(a+1,b-1)*temp(3,1)+Im(a+1,b)*temp(3,2)+Im(a+1,b+1)*temp(3,3);

end

end

I=In;

Im=zeros(m,n);

for x=5:m-5;

for y=5:n-5;

sum1=I(x,y-4)+I(x,y-2)+I(x,y+2)+I(x,y+4);

sum2=I(x-2,y+4)+I(x-1,y+2)+I(x+1,y-2)+I(x+2,y-4);

sum3=I(x-2,y+2)+I(x-4,y+4)+I(x+2,y-2)+I(x+4,y-4);

sum4=I(x-2,y+1)+I(x-4,y+2)+I(x+2,y-1)+I(x+4,y-2);

sum5=I(x-2,y)+I(x-4,y)+I(x+2,y)+I(x+4,y);

sum6=I(x-4,y-2)+I(x-2,y-1)+I(x+2,y+1)+I(x+4,y+2);

sum7=I(x-4,y-4)+I(x-2,y-2)+I(x+2,y+2)+I(x+4,y+4);

sum8=I(x-2,y-4)+I(x-1,y-2)+I(x+1,y+2)+I(x+2,y+4);

sumi=[sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8];

summax=max(sumi);

summin=min(sumi);

summ=sum(sumi);

b=summ/8;

if (summax+summin+ 4*I(x,y))> (3*summ/8)

sumf = summin;

else

sumf =summax;

end

if sumf > b

Im(x,y)=128;

else

Im(x,y)=255;

end

end

end

for i=1:m

for j =1:n

Icc(i,j)=Icc(i,j)*Im(i,j);

end

end

for i=1:m

for j =1:n

if (Icc(i,j)==128)

Icc(i,j)=0;

else

Icc(i,j)=1;

end

end

end

figure,imshow(double(Icc));title('二值化');

%************************************************************************

u=Icc;

[m,n]=size(u) %去空洞和毛刺

for x=2:m-1

for y=2:n-1

if u(x,y)==0

if u(x,y-1)+u(x-1,y)+u(x,y+1)+u(x+1,y)>=3

u(x,y)=1;

end

else u(x,y)=u(x,y);

end

end

end

figure,imshow(u) ;title('去毛刺')

for a=2:m-1

for b=2:n-1

if u(a,b)==1

if abs(u(a,b+1)-u(a-1,b+1))+abs(u(a-1,b+1)-u(a-1,b))+abs(u(a-1,b)-u(a-1,b-1))+abs(u(a-1,b-1)-u(a,b-1))+abs(u(a,b-1)-u(a+1,b-1))+abs(u(a+1,b-1)-u(a+1,b))+abs(u(a+1,b)-u(a+1,b+1))+abs(u(a+1,b+1)-u(a,b+1))~=1 %去空洞

if(u(a,b+1)+u(a-1,b+1)+u(a-1,b))*(u(a,b-1)+u(a+1,b-1)+u(a+1,b))+(u(a-1,b)+u(a-1,b-1)+u(a,b-1))*(u(a+1,b)+u(a+1,b+1)+u(a,b+1))==0 %去毛刺

u(a,b)=0;

end

end

end

end

end

figure,imshow(u) ;title('去空洞')
%*************************************************************************
v=~u;

se=strel('square',3);

fo=imopen(v,se);

v=imclose(fo,se); %对图像开操作和闭操作

img=bwmorph(v,'thin',Inf); %对图像进行细化

figure,imshow(img)
title('细化图')
end

function j = P (img, x, y, i)
switch (i)
case {1, 9}
j = img(x+1, y);
case 2
j = img(x + 1, y-1);
case 3

j = img(x, y - 1);

case 4

j = img(x - 1, y - 1);

case 5

j = img(x - 1, y);

case 6

j = img(x - 1, y + 1);

case 7

j = img(x, y + 1);

case 8

j = img(x + 1, y + 1);
end
end

%point函数
function txy=point(thin)
count = 1;
txy(count, :) = [0,0,0];
siz=min(size(thin,1),size(thin,2));
for x=40:siz - 40
for y=40:siz - 40
if (thin(y, x) )
CN = 0;
for i = 1:8
CN = CN + abs (P(thin, y, x, i) - P(thin, y, x, i + 1));
end
if (CN == 2)
txy(count, :) = [x, y,2];
count = count + 1;
end
if (CN == 6)
txy(count, :) = [x, y,6];
count = count + 1;
end
end
end
end
for i=1:count - 1
x(i) =txy(i, 1);
y(i)= txy(i, 2);
end
imshow(double(thin));
hold on;
plot(x,y,'.');
end
%cut函数
function txy=cut(thin,txy)
s(8,8)=0;
delta(8,8)=0;
n=size(txy,1);
txy=txy(find(txy(:,1)),:);
plot(txy(:,1),txy(:,2),'ro');
end


%% 指纹图像特征点代码
%single_point函数
function [pxy2,error]=single_point(txy,r)

error=0;

x=txy(:,1);

y=txy(:,2);

n=length(x);

d(1:n,1:n)=0;

for j=1:n

for i=1:n

if (i~=j)

d(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);

else

d(i,j)=2*r;

end

end

end

[a,b]=min(d);

c=find(a>r);

pxy2=txy(c,:);

pxy2=pxy2(find(pxy2(:,3)==2),:);

t=size(pxy2,1);

if t==0

error=1

else

plot(x,y,'b.');

hold on

plot(pxy2(:,1),pxy2(:,2),'r.');

end
end
%walk函数
function [error,a,b]=walk(thin,x0,y0,num)

error=0;

thin(y0,x0)=0;

t1=0;

for n=1:num

if error==1

break;

else

x=x0;

y=y0;

for x=x0-1:x0+1

if error==1

break;

else

for y=y0-1:y0+1

t1=sum(sum(thin(y0-1:y0+1,x0-1:x0+1)));

if (t1==0||t1>=2)

error=1;

a=x0;

b=y0;

break;

else

if (thin(y,x)==1&&(x-x0)^2+(y-y0)^2~=0)

if (t1>=2 )

error=1;

break ;

else

thin(y,x)=0;

x0=x;

y0=y;

a=x0;

b=y0;

plot(x0,y0,'r.')

end

end

end

end

end

end

end

end
end
%last1函数
function [pxy3,error2]=last1(thin,r,txy,num)

error=0;

[pxy2,error]=single_point(txy,r);

n=size(pxy2,1);

l=1;

error2=0;
end
%% 特征点匹配代码
%distance函数程序
function d=distance(x0,y0,num,thin)

num2=fix(num/5);

for i=1:num2

[error,a,b]=walk(thin,x0,y0,5*i);

if error~=1

d(i)=sqrt((a-x0)^2+(b-y0)^2);

else

break;

end

end
end
%A函数
function pxy=A(x0,y0,txy,num)

x=txy(:,1);

y=txy(:,2);

n=length(x);

l(1,n)=0;

lnn=1;

pxy(num,:)=[0,0,0];

for i=1:n

l(i)=sqrt((x(i)-x0)^2+(y(i)-y0)^2);

end

ll=sort(l);

for i=1:num

xiao=ll(i+lnn);

nn=find(l==xiao);

lnn=length(nn);

pxy(i,:)=[x(nn(1)),y(nn(1)),txy(nn(1),3)];

end

plot(x0,y0,'bo');

x0;

y0;

hold on

plot(pxy(:,1),pxy(:,2),'ro');
end
