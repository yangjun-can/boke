clc;
clear;
P=20;N=512;
er_m=0.5;%|er_m|<=0.5
K=1024;%K>N
Km=512;%Km=0:K-1
bi=5;
t=0:511;
K=length(t);
w1=0.2*pi; 
w2=0.3*pi;
x1=150*sin(w1*t);
x2=276*sin(w2*t);
SigK = x1 + x2;
t2=0:255;
n=length(t2);
x12=150*sin(w1*t2);
x22=276*sin(w2*t2);
SigN = x12 + x22;

for i=0:P
    y=(bi*(pi*(2*n+1-N)/K)^i).*SigN;
end
for n=0:N-1
    Bi=y*exp(2*pi*n*Km/K);
end

% for i=0:P
%     for n=0:N-1
%         Xm=exp(pi*er_m*(n-1)/K)*er_m^i*(bi*(pi*(2*n+1-N)/K)^i)*SigN*exp(2*pi*n*Km/K);
%     end
% end