function [Sig, A, u, v] = genSigOnGrid(N0,N1,K,a_min,sigma_n)
% Generate 2-D signals with on-grid uniformly distributed frequencies randomly 

% Input: 
%   N0, N1: signal length of the two dimensions
%   K: sparsity level
%   a_min: the minimum sigal amplitude
%   sigma_n: noise variance 

% Output: 
%   Sig: signal matrix (N0-by-N1)
%   A: amplitude of each frequency (K-length vector)
%   u,v: frequency locations

%% 
%随机生成频率的K个格点位置
    u = randi(N0,[K,1])-1;% 随机生成K个（0，N0-1）之间均匀分布的数；（行）
    v = randi(N1,[K,1])-1;% 随机生成K个（0，N1-1）之间均匀分布的数；（列）
      %rand()函数用于生成取值在（0~1）之间均匀分布的伪随机数。
      %randn()函数用于生成标准正态分布的伪随机数。
      %randi(imax,m,n)函数用于生成均匀分布的m*n的伪随机整数矩阵，范围为（imin~imax）（开区间），若imin缺省，默认为1。
    ind = [u,v];
    uInd = unique(ind,'rows');%频率位置（每个位置只出现一次）
    %unique用法：
    %1.C = unique(A) 返回与A中相同的数据，但是不包含重复项。且向量C已按照从小到大排序。
    %2.C = unique(A,___,'rows') 和 C = unique(A,'rows',___) 将A中的每一行视为单个实体，并按排序顺序返回A中的唯一行。必须指定A，而setOrder和occurrence是可选的。'rows' 选项不支持元胞数组。
    %3.[C,ia,ic] = unique(A,___) 还可使用上述任何语法返回索引向量 ia 和 ic。ia为矩阵C中的元素在矩阵A中的位置，ic为矩阵A中的元素在矩阵C中的位置。
    while size(uInd,1)<K
        %若频率位置个数（第一列的行数）<K,继续随机生成位置。
        u = randi(N0,[K,1])-1;
        v = randi(N1,[K,1])-1;
        uInd = [uInd; [u,v]];
        uInd = unique(uInd,'rows');
    end
    rdI = randperm(length(uInd));%将频率位置随机打乱顺序
    %y = randperm(n) :把1到n这些数随机打乱得到的一个数字序列（行向量）。
    %y=randperm(n,m)，1-n中随机选择m个数，n一定大于等于m
    %size 获取数组的行数和列数
    %length 获取数组长度，即行数和列数中的较大值，相当于max(size(a))
    %numel 返回元素总数
    uInd = uInd(rdI(1:K),:);%取定前K个位置    
    oriInd = sortrows(uInd);%对位置按行升序排列
    u = oriInd(:,1);%（行位置）
    v = oriInd(:,2);%（列位置）
%%
 %生成频域信号及时域信号
    A = a_min*exp(1i*2*pi*rand(K,1));%随机生成信号在频域的振幅；1i复数，pi圆周率
    Sig=zeros(N0,N1);
    Sig(sub2ind([N0,N1],u+1,v+1)) = A; %稀疏的频域信号  

    %ind2sub把数组或者矩阵的线性索引转化为相应的下标；
    %sub2ind则正好相反，将下标转化为线性索引。ind = sub2ind(sz,row,col) 针对大小为 sz 的矩阵返回由 row 和 col 指定的行列下标的对应线性索引 ind    
    %Sig = ifft2(Sig);%对应的时域信号   %iff2 二维快速傅里叶逆变换
    %Sig = fft2(Sig);%对应的频域信号  
    p1=0.8;p2=0.8;ts1=0.005;ts2=0.005;
    alpha=p1*pi/2;     S1=sign(sin(alpha));
    beta=p2*pi/2;     S2=sign(sin(beta));
    us1=2*pi*S1*sin(alpha)/(ts1*N0);
    us2=2*pi*S2*sin(beta)/(ts2*N1);
    [Sig,~,~] = DFRFT_2D_fft2(Sig,-p1,-p2,N0,N1,us1,us2); %时域信号
    if  sigma_n==1
    noise = sigma_n*sqrt(2)/2*(randn(N0,N1)+1i*randn(N0,N1));  % 时域加复高斯噪声
    % 计算时域信噪比
%     Ps=sum(abs(Sig).^2);
%     Pn=sum(abs(noise).^2);
%     SNR_time=10*log10(Ps/Pn)
    SNR_t=snr(Sig,noise)
    Sig = Sig+noise;%含噪时域信号
    end
%% Visualization       
%     figure,mesh(0:N1-1,0:N0-1,abs(fft2(Sig)));
    %abs函数:数值的绝对值和复数的幅值
    %mesh语句画网格图片，实际上就是给出一对坐标(x,y)，来画矩阵z(x,y)的值。
%     xlabel('$v$','Interpreter','LaTex'), ylabel('$u$','Interpreter','LaTex'), zlabel('Amplitude');
%     axis tight; % 自动设置x轴和y轴的范围使图形区域正好占满整个显示空间
%     set(gca,'FontSize',20);%x,y轴标注文字都会改变大小
end