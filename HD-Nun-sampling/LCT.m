% Pei DLCT  一维离散线性正则变换
function [X_lct,us] = LCT(x,a,b,c,d, N, ts)
%% 参数
% x : signal
% a,b,c,d : parameters
% N :  signal size
% ts : Sampling interval
% X_lct: DLCT of x
% us :  Sampling interval in LCTD
%%
M=N;  % signal size in LCTD
us=2*pi*abs(b)/(ts*M);
% n=-N/2+1:N/2;
% m=-M/2+1:M/2;  
n=0:N-1;
m=0:M-1;
% A=[a,b;c,d];
% if det(A)==1.0000
    %% 特殊情况
    % Special cases 1: DFT
    if a==0 && b==1 && c==-1 && d==0
        X_lct=fft(x);
        return;
    end
    % Special cases 2
    if a==0 && b==-1 && c==1 && d==0
        X_temp=fft(x);
        X_lct = X_temp(M:-1:1);
        return;
    end
     % Special cases 3: b=0
    if b==0      
        u=us*m;
        remainder=mod(d,1); % 取模  
        if remainder==0  % d is an intager
            X_lct = sqrt(d)*exp(1i/2*c*d*u.^2)*x(d*m); 
        else
            y=sqrt(1/N)*fft(x);
            us=N*abs(a)*ts/M;
            u=us*m;
            z=sqrt(1/M)*y*exp(1i*c/(2*a)*u.^2);
            if a<0
                X_lct =fft(z);
            else
                X_lct =ifft(z);
            end
        end
        return;
    end
    %%  DLCT
    % 第一次chirp调制
    t=ts*n;
    chirp_t = exp(1i*a/(2*b)*t.^2);
    x_t=x.*chirp_t;
    % FFT
    if b>0
        X_lct=fft(x_t);
    else
        X_lct=ifft(x_t)*N;
    end
    % 第二次chirp调制
    u=us*m;
    chirp_u = exp(1i*d/(2*b)*u.^2);
    X_lct=sqrt(1/M)*chirp_u.*X_lct;
% k1=1/sqrt(1i*2*pi*b);
% X_lct=k1*chirp_u.*X_lct;
end
