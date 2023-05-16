function [Omega,A] =MARS_SFT_corr(Sig, Win, N0, N1, T, epsilon, gamma, gama)                     
% Robutst MARS_SFT algorithm; 定位误差矫正法

% Input:
%   Sig: singal matrix 
%   Win: window matrix
%   N0, N1: signal length of the two dimensions
%   epsilon: the threshold of detecting significant frequencies in a slice
%   gamma: the threshold for 1-sparsity detection
%   gama： 定位误差（由相位误差的KDE估计得到）

% Output: 
%   Omega: frequency locations
%   A: complex amplitude of each frequency 

Omega = [];
A = [];
i = 0;  
while i<T  % 迭代T次
    i = i+1;
    [ind,ha] = SFT_INNER_corr(Sig,Win,Omega,A,N0,N1,epsilon, gamma,gama);
    % ind 内循环估计的位置
    % ha 内循环估计的幅度
    if ~isempty(ha)    % 这次内循环有被估计的大频率
        if isempty(A)  % 之前没有被估计的大频率
            Omega = [Omega;ind]; 
            A = [A; ha];
        else
            [C,ia,ib] = intersect(Omega,ind,'rows');            
           % intersect(A,B,___,'rows') 和 C = intersect(A,B,'rows',___) 将A和B的每一行都视为单个实体，并返回A和B的共有行，但是不包含重复项。
           % 必须指定 A 和 B，setOrder 是可选的。C = A(ia,:) 且 C = B(ib,:)。
           % C 重复估计的位置  ia重复位置在A中的坐标， ib重复位置在B中的坐标
            if ~isempty(C)  % 如果有重复定位，该位置的幅度叠加
                A(ia) = A(ia)+ha(ib);
                df = setdiff(1:size(ind,1), ib);
                %c = setdiff(A, B) 返回在A中有，而B中没有的值，结果向量将以升序排序返回。
                Omega = [Omega; ind(df,:)]; % 将剩余不重复的位置加入已估计集
                A = [A; ha(df)];
            else       % 如果没有重复定位，将所有位置加入已估计集
                Omega = [Omega;ind];
                A = [A; ha];
            end            
        end    
    end   
end
end

function [Omega, A] = SFT_INNER_corr(Sig, Win, Omega1, A1, N0, N1, epsilon, gamma,gama)
% 随机切片估计含噪的频率  误差矫正法
% Input:
%   Sig: singal matrix
%   Win: window matrix
%   Omega1, A1: location and amplitude of previously recovered frequenceis
%   N0, N1: signal length of the two dimensions
%   epsilon: the threshold of detecting significant frequencies in a slice
%   gamma: the threshold for 1-sparsity detection

% Output:
%   Omega: frequency locations
%   A: complex amplitude of each frequency
%   P: debug info 调试信息

Omega = [];
A = [];
L = round(lcm(N0,N1)); % 切片长度
c0 = L/N0;
c1 = L/N1;
%% 构造随机切片的参数
tau0 = randi(N0,1)-1; % randi(imax,n)产生一个n*n矩阵，这个矩阵的元素都是小于等于imax的正整数。
tau1 = randi(N1,1)-1;
alpha0 = randi(N0,1)-1;
alpha1 = randi(N1,1)-1;
while gcd(alpha0,alpha1) ~= 1 || gcd(alpha0,c1) ~= 1 || gcd(alpha1,c0) ~= 1
    % gcd 求两个整数的最大公约数，返回值是正整数
    alpha0 = randi(N0,1)-1;
    alpha1 = randi(N1,1)-1;
end
disp([' tau0 = ',num2str(tau0)]);
disp([' tau1 = ',num2str(tau1)]);
disp([' alpha0 = ',num2str(alpha0)]);
disp([' alpha1 = ',num2str(alpha1)]);
%% 利用切片估计
[ind,ha] = triSlicing2(Sig, Win, Omega1, A1 ,alpha0, alpha1, tau0,tau1, epsilon, gamma,L,gama);
% ind 储存估计频率的坐标  ha 储存估计频率的值  res 切片的平均能量
% 如果此次估计的频率集非空，将估计结果进行储存和统计
if ~isempty(ha)
    Omega = [Omega;ind];
    A = [A; ha.'];
end
end


%% 利用三个平行切片的DFT估计频率
function [ind, ha] = triSlicing2(Sig, Win, I1, hA1, alpha0, alpha1, tau0, tau1, epsilon, gamma,L,gama)
[N0,N1] = size(Sig);
u = [];   % 储存估计频率的横坐标
v = [];   % 储存估计频率的纵坐标
ha = [];  % 储存估计频率的幅值
Icand=[]; % 储存定位矫正时可能的频率坐标
%% 三次切片分桶
[in0, hs0] = SLICING(Sig, Win, I1, hA1, alpha0, alpha1, tau0, tau1, epsilon,L);    % in 大桶下标
[in1, hs1] = SLICING(Sig, Win, I1, hA1, alpha0, alpha1, tau0+1, tau1, epsilon,L);  % hs 大桶桶值
[in2, hs2] = SLICING(Sig, Win, I1, hA1, alpha0, alpha1, tau0, tau1+1, epsilon,L);  % res 切片总能量（不含已估计频率）

% bucket0=zeros(1,L);
% bucket0(in0) = hs0;
% bucket1=zeros(1,L);
% bucket1(in1) = hs1;
% bucket2=zeros(1,L);
% bucket2(in2) = hs2;
% figure;
% subplot(131); plot((0:L-1)/(2^-4),abs(bucket0)); xlabel('Buckets ${{\hat Q}_{{\sigma _1},{\sigma _2},{\tau _1},{\tau _2}}}$', 'Interpreter','latex');ylabel('amplitude');
% set(gca,'FontSize',15)
% subplot(132); plot((0:L-1)/(2^-4),abs(bucket1)); xlabel('Buckets ${{\hat Q}_{{\sigma _1},{\sigma _2},{\tau _1+1},{\tau _2}}}$', 'Interpreter','latex');ylabel('amplitude');
% set(gca,'FontSize',15)
% subplot(133); plot((0:L-1)/(2^-4),abs(bucket2)); xlabel('Buckets ${{\hat Q}_{{\sigma _1},{\sigma _2},{\tau _1},{\tau _2+1}}}$', 'Interpreter','latex');ylabel('amplitude');
% set(gca,'FontSize',15)

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
%%  对于1稀疏的大桶，估计频率的位置和值
ahs = abs(hs);
times_all=0;% 统计内循环次数
for i = 1:length(in)  % 对于大桶i 如果是1稀疏的，估计大频率
    % 1-sparse detection
    if var(ahs(:,i))< gamma
        u0 = round(wrapTo2Pi(angle(hs(2,i)/hs(1,i)))*N0/(2*pi)); % 大频率的横坐标 % angle 相位角 wrapTo2Pi弧度角转换到[0,2pi]
        v0 = round(wrapTo2Pi(angle(hs(3,i)/hs(1,i)))*N1/(2*pi));
        % 定位矫正
        u_corr=[]; % 可能的行位置集
        v_corr=[]; % 可能的列位置集
        for u_prob=(u0-gama):(u0+gama)
            for v_prob=(v0-gama):(v0+gama)
                u1=mod(u_prob, N0);
                v1=mod(v_prob, N1);
                if mod(L*alpha0/N0*u1+L*alpha1/N1*v1,L)==in(i)
                    u_corr(end+1)=u1;
                    v_corr(end+1)=v1;
                end
            end
        end
        if isempty(u_corr)
            continue;  % 没有找到可能的位置集：1.大桶检测错误，不存在大频率 2.
        elseif length(u_corr)==1     % 若只有一个可能的位置，记为恢复的位置
            u(end+1) = u_corr(1);    % 储存大频率的横坐标（消除噪声影响，整数） end+1:在当前的基础上增加一列
            v(end+1) = v_corr(1);    % 储存大频率的纵坐标
            ha(end+1) = ((hs(1,i)+ hs(2,i)*exp(-1i*2*pi*u(end)/N0) + ...   % 一行写不下时,可以用省略号开启另外一行
                hs(3,i)*exp(-1i*2*pi*v(end)/N1)))*N0*N1/L*exp(-1i*2*pi*(u(end)*tau0/N0+v(end)*tau1/N1))/3;  % 储存大频率的值（三次切片的均值）
        else      %  若有多个可能的位置，重新分桶并删去分在小桶中的位置,直到仅剩下一个位置
            Icand = [u_corr;v_corr].'; % 储存待定位置集Icand
            Icand = unique(Icand,'rows'); % 所有的可能位置集
            times=0;
            while size(Icand,1)~=1
                %% 找到合适的参数，使得Icand中的位置分到不同的桶
                times=times+1;
                if times>=10
                     break; 
                end
                in_cand=zeros(size(Icand,1),1); % 存储位置集Icand的分桶结果
                while length(unique(in_cand))~=length(in_cand)  % 各位置要分在不同的桶中
                    alpha10 = randi(N0,1)-1;
                    alpha11 = randi(N1,1)-1;
                    while gcd(alpha10,alpha11) ~= 1 || gcd(alpha10,L/N1) ~= 1 || gcd(alpha11,L/N0) ~= 1
                        % gcd 求两个整数的最大公约数，返回值是正整数
                        alpha10 = randi(N0,1)-1;
                        alpha11 = randi(N1,1)-1;
                    end
                    const=[L*alpha10/N0;L*alpha11/N1];
                    in_cand=mod(Icand*const,L); % 各位置分桶结果
                end
                %% 按照选取的参数进行分桶，仅保留仍分在大桶中的位置
                [in3, ~] = SLICING(Sig, Win, I1, hA1, alpha10, alpha11, tau0, tau1, epsilon,L);  % in3大桶位置 hs3大桶桶值
                [~, ia3_cand,~]= intersect(in_cand, in3); % 大桶
                Icand=Icand(ia3_cand,:);
            end
            times_all=times_all+times;
            %disp([' 该点的内循环次数 = ',num2str(times)]);
            if size(Icand,1)==1
            u(end+1) = Icand(1,1);    % 储存大频率的横坐标（消除噪声影响，整数） end+1:在当前的基础上增加一列
            v(end+1) = Icand(1,2);    % 储存大频率的纵坐标
            ha(end+1) = ((hs(1,i)+ hs(2,i)*exp(-1i*2*pi*u(end)/N0) + ...   % 一行写不下时,可以用省略号开启另外一行
                hs(3,i)*exp(-1i*2*pi*v(end)/N1)))*N0*N1/L*exp(-1i*2*pi*(u(end)*tau0/N0+v(end)*tau1/N1))/3;  % 储存大频率的值（三次切片的均值）
            end
        end
      % ha(end+1) = hs(1,i)*exp(-1i*2*pi*(u(end)*tau0/N0+v(end)*tau1/N1));
    end
end
disp([' 总内循环次数 = ',num2str(times_all)]);
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
%     % add window
%     win = Win(ind);
%     sl = sl.*win; % 加窗的信号

% construct the line from found frequencies
cl = zeros(1, L);  % 储存已估计频率的分桶值
if ~isempty(hA)  % 若有已估计频率，找出其时域影响
    cl=CONSTRUCTION(I,hA,alpha0,alpha1,N0,N1,tao0,tao1,L);
end

% detection of significant frequencies
fsl = fft(sl)-cl; % 切片的DFT
% figure;plot((0:L-1)/(2^-4),abs(fsl)); 
% xlabel('Buckets ${{\hat Q}_{{\sigma _1},{\sigma _2},{\tau _1},{\tau _2}}}$', 'FontSize',15,'Interpreter','latex');
% ylabel('amplitude','FontSize',15);



ft= abs(fsl).^2;
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

%% 已估计频率的影响
function f = CONSTRUCTION(I,f_hat,alpha0,alpha1,N0,N1,tao0,tao1,L)
f_hat = L/N0/N1*f_hat.* exp(1i*2*pi*(I(:,1)/N0*tao0 + I(:,2)/N1*tao1));  % 已估计频率按照参数alpha0,alpha1,tao0,tao1的投影值
k= mod(round(I*[alpha0*L/N0 alpha1*L/N1]'),L)+1;   % 已估计频率按照参数alpha0,alpha1,tao0,tao1的分桶
fhat1=accumarray(k,f_hat,[L,1]); % 对分到同一桶中的投影值累加，得到总的已估计频率的投影向量
f=fhat1.';
end
