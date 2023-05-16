function  freq = Min_Max (SigN,N,f)
% step1:对上采样点加权并FFT
L=1;
K=10*N; % 上采样点数
J=6; % 用于插值的点的个数
beta=0.19;
galma=2*pi/K;
alfa=[conj(-0.46),1,-0.46];
eita=(N-1)/2;
Weights=alfa(1)*exp(-1i*galma*beta*((0:N-1)-eita))+alfa(2)+alfa(3)*exp(1i*galma*beta*((0:N-1)-eita));
Weights_SigN=Weights.*SigN;
UpSample_freq=fft(Weights_SigN,K);
% step2 precomputed
gama=zeros(J,J);
T=zeros(J,J);
r=zeros(J,1);
freq=zeros(1,N);
for cc=1:J
    for j1=1:J
T(cc,j1)=alfa(1)*( conj(alfa(1))*Periodic_sinc(j1-cc,K,N) + conj(alfa(2))*Periodic_sinc(j1-cc-beta,K,N) + conj(alfa(3))*Periodic_sinc(j1-cc-2*beta,K,N) ) ...
        +alfa(2)*( conj(alfa(1))*Periodic_sinc(j1-cc+beta,K,N) + conj(alfa(2))*Periodic_sinc(j1-cc,K,N) + conj(alfa(3))*Periodic_sinc(j1-cc-beta,K,N) ) ...
        +alfa(3)*( conj(alfa(1))*Periodic_sinc(j1-cc+2*beta,K,N) + conj(alfa(2))*Periodic_sinc(j1-cc+beta,K,N) + conj(alfa(3))*Periodic_sinc(j1-cc,K,N) );
    end
end
% step2: 插值
for point=1:N
    % 定义用于插值point的第一个点
    if mod(J,2)~=0
        [Int_Offset,~]=close_to_freq(f(point),K,N,J);
    else
        Int_Offset=floor(f(point)*K/N)-J/2;  
    end
    for j2=1:J
       r(j2)=alfa(1)* Periodic_sinc( f(point)*K/N-Int_Offset-j2-beta ,K,N) ...
           + alfa(2)* Periodic_sinc( f(point)*K/N-Int_Offset-j2,K,N ) ...
           + alfa(3)* Periodic_sinc( f(point)*K/N-Int_Offset-j2+beta,K,N );
       gama(j2,j2)=exp(-1i*( 2*pi*f(point)/N -galma*(Int_Offset+j2) )*eita);
    end
    u=gama.'*inv(T)*r;  %插值的系数向量
    for j3=1:J
      freq(point)= UpSample_freq(mod(Int_Offset+j3,K)+1) * conj(u(j3)) +freq(point);      
    end
end   
end