function [hat_x]=csreconstruction(x,N,K)
M=64;     %  测量数
Phi=randn(M,N);                                   %  测量矩阵(高斯分布白噪声)
% s=Phi*x_noise.';                                        %  获得线性测量 
s=Phi*x.'; 
m=2*K;                                            %  算法迭代次数(m>=K)
Psi=fft(eye(N,N))/sqrt(N);                        
T=Phi*Psi';                                      
hat_y=zeros(1,N);                                 %  待重构的谱域(变换域)向量                     
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=s;                                                       

for times=1:m                                 
    for col=1:N                               
        product(col)=abs(T(:,col)'*r_n);           
    end
    [val,pos]=max(product);                      
    Aug_t=[Aug_t,T(:,pos)];                      
    T(:,pos)=zeros(M,1);                         
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;          
    r_n=s-Aug_t*aug_y;                            %  残差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
end
hat_y(pos_array)=aug_y;                           %  重构的谱域向量
hat_x=real(Psi'*hat_y.');                         %  做逆傅里叶变换重构得到时域信号
end