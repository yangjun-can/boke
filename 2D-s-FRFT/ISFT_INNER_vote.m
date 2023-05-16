function [Omega, A] = ISFT_INNER_vote(Sig, Win, Omega1, A1, N0, N1, epsilon, gamma, n_d, n_s)
% Input:
%   Sig: singal matrix 
%   Win: window matrix
%   Omega1, A1: location and amplitude of previously recovered frequenceis
%   N0, N1: signal length of the two dimensions
%   epsilon: the threshold of detecting significant frequencies in a slice
%   gamma: the threshold for 1-sparsity detection
%   n_d, n_s: parameters of $n_d$-out-of-$n_s$ detection

% Output: 
%   Omega: frequency locations
%   A: complex amplitude of each frequency 
%   P: debug info 调试信息

S = Sig; 
% P = [];
Omega = [];
A = [];
oc = [];
it = 0;
L = round(lcm(N0,N1)); % 切片长度 
c0 = L/N0;
c1 = L/N1; 
%% n_d次内循环
while it < n_d   
    it = it+1;    
    %% 构造随机切片的参数
    tau0 = randi(N0,1)-1;
    tau1 = randi(N1,1)-1;
    alpha0 = randi(N0,1)-1; 
    alpha1 = randi(N1,1)-1;   
      % randi(imax,n)产生一个n*n矩阵，这个矩阵的元素都是小于等于imax的正整数。
    while gcd(alpha0,alpha1) ~= 1 || gcd(alpha0,c1) ~= 1 || gcd(alpha1,c0) ~= 1
        % gcd 求两个整数的最大公约数，返回值是正整数
        alpha0 = randi(N0,1)-1; 
        alpha1 = randi(N1,1)-1; 
    end
    %% 利用切片估计
    [ind,ha] = triSlicing2(S, Win, Omega1, A1 ,alpha0, alpha1, tau0,tau1, epsilon, gamma,L);
        % ind 储存估计频率的坐标  ha 储存估计频率的值  
%     P{it,1} = [alpha0,alpha1];
%     P{it,2} = [tau0, tau1];
%     P{it,3} = ind;
%     P{it,4} = ha;
    %% 如果此次估计的频率集非空，将估计结果与其他内循环的估计结果进行储存和统计
    if ~isempty(ha)
        if isempty(A)     % 之前内循环没有估计频率，直接将此次结果并入，投票数记为1
            Omega = [Omega;ind];
            A = [A; ha.'];
            oc = [oc; ones(length(ha),1)];
        else             % 之前内循环有估计频率，相同的估计位置幅度相加 投票数+1 ，剩余的直接将结果加入
            [C,ia,ib] = intersect(Omega,ind,'rows');   
            if ~isempty(C)
                A(ia) = A(ia)+ha(ib).';  % 重复的估计位置幅度相加
                oc(ia) = oc(ia)+1;       % 投票数+1 
                df = setdiff(1:size(ind,1), ib);  % 估计集中除去重复频率的下标 
                Omega = [Omega; ind(df,:)];
                A = [A; ha(df).'];
                oc = [oc; ones(length(df),1)];
            else
                Omega = [Omega;ind];
                A = [A; ha.'];
                oc = [oc; ones(length(ha),1)];
            end                      
        end       
    end
end

%visual_votes(V,3);
i = find(oc>=n_s);   % 投票次数超过n_s的位置，返回其坐标及平均幅度
Omega = Omega(i,:);
A = A(i,:)./oc(i);

end
%% 利用三个平行切片的IDFT估计频率
function [ind, ha] = triSlicing2(Sig, Win, I1, hA1, alpha0, alpha1, tau0, tau1, epsilon, gamma,L)
    [N0,N1] = size(Sig);
    u = [];   % 储存估计频率的横坐标
    v = [];   % 储存估计频率的纵坐标
    ha = [];  % 储存估计频率的幅值
    %% 三次切片分桶
    [in0, hs0] = SLICING(Sig, Win, I1, hA1, alpha0, alpha1, tau0, tau1, epsilon,L);    % in 大桶下标
    [in1, hs1] = SLICING(Sig, Win, I1, hA1, alpha0, alpha1, tau0+1, tau1, epsilon,L);  % hs 大桶桶值
    [in2, hs2] = SLICING(Sig, Win, I1, hA1, alpha0, alpha1, tau0, tau1+1, epsilon,L);  
    %% 大桶的位置和桶值
    in = intersect(in0, in1);   
    in = intersect(in, in2);    % 找到三次切片都是大桶的桶坐标
    
    hs = zeros(3, length(in));  % 存储大桶的三次桶值，每一行代表一次分桶    
    if length(in) ~= length(in0) || length(in) ~= length(in1) || length(in) ~= length(in2)
        [~,ia0,~] = intersect(in0,in); % 三次切片都是大桶的桶在大桶集in0中的位置
        [~,ia1,~] = intersect(in1,in);
        [~,ia2,~] = intersect(in2,in);        
        hs(1,:) = hs0(ia0);            % 三次切片都是大桶的桶在第一次切片中的值
        hs(2,:) = hs1(ia1);
        hs(3,:) = hs2(ia2);        
    else                               % 若三次找到的大桶相同，直接存储相应桶值
        hs(1,:) = hs0;
        hs(2,:) = hs1;
        hs(3,:) = hs2;
    end
 %%     对于1稀疏的大桶，估计频率的位置和值
    ahs = abs(hs);
    for i = 1:length(in)  % 对于大桶i 如果是1稀疏的，估计大频率
        % 1-sparse detection
        if var(ahs(:,i))< gamma
            u0 = wrapTo2Pi(-angle(hs(2,i)/hs(1,i)))*N0/(2*pi); % 大频率的横坐标 % angle 相位角 wrapTo2Pi弧度角转换到[0,2pi]
            u(end+1) = mod(round(u0), N0);    % 储存大频率的横坐标（消除噪声影响，整数） end+1:在当前的基础上增加一列
            v0 = wrapTo2Pi(-angle(hs(3,i)/hs(1,i)))*N1/(2*pi); 
            v(end+1) = mod(round(v0), N1);    % 储存大频率的纵坐标
            %ha(end+1) = hs(1,i)*exp(-1i*2*pi*(u(end)*tau0/N0+v(end)*tau1/N1));
            ha(end+1) = ((hs(1,i)+ hs(2,i)*exp(1i*2*pi*u(end)/N0) + ...   % 一行写不下时,可以用省略号开启另外一行
                hs(3,i)*exp(1i*2*pi*v(end)/N1)))*exp(1i*2*pi*(u(end)*tau0/N0+v(end)*tau1/N1))/3;  % 储存大频率的值（三次切片的均值）
        end
    end
    ind = [u; v]';             
end

%% 切片（减去已估计频率的影响）并返回大桶的位置和桶值
function [in, hs] = SLICING(Sig, Win, I, hA,alpha0, alpha1, tao0,tao1, epsilon,L)
    [N0, N1]=size(Sig);
    % sampling on a line
    l = 0:L-1;
    x = mod(alpha0*l+tao0,N0)+1;
    y = mod(alpha1*l+tao1,N1)+1;
    ind = sub2ind([N0,N1], x, y);  % 返回位置（x，y）在大小[N0,N1]的矩阵中的按列索引号
    sl = Sig(ind); % 按照参数获取到的切片
    
    % add window
%     win = Win(ind);
%     sl = sl.*win; % 加窗的信号
    
    % construct the line from found frequencies
    cl = zeros(1, L);  % 储存已估计频率的系数
    if ~isempty(hA)  % 若有已估计频率，找出其时域影响
        cl=CONSTRUCTION(I,hA,alpha0,alpha1,N0,N1,tao0,tao1,L); 
    end
    dl = sl; % 时域信号
    
    % detection of significant frequencies
    fsl = ifft(dl)-cl; % 切片的IDFT减去已估计频率的影响
    ft= abs(fsl);
    
%     figure,stem(0:N0-1,ft,'LineWidth',2),grid;
%     xlabel('$\omega$','Interpreter','LaTex');
%     ylabel('Amplitude');
%     %legend('DFT grids','Ground truth','Estimated');
%     set(gca,'FontSize',20);
%     axis tight;
    
    in = find(ft>=epsilon); % the first detection 大桶
    hs = fsl(in); % 大桶的桶值
    in  = in-1;  % 大桶的坐标
end

%% 已估计点的影响
function f = CONSTRUCTION(I,f_hat,alpha0,alpha1,N0,N1,tao0,tao1,L)
    f_hat = f_hat.* exp(-1i*2*pi*(I(:,1)/N0*tao0 + I(:,2)/N1*tao1));  % 已估计频率按照参数alpha0,alpha1,tao0,tao1的投影值
    k= mod(round(I*[alpha0*L/N0 alpha1*L/N1]'),L)+1;   % 已估计频率按照参数alpha0,alpha1,tao0,tao1的分桶
    fhat1=accumarray(k,f_hat,[L,1]); % 对分到同一桶中的投影值累加，得到总的已估计点的投影向量
    f=fhat1.';
end
