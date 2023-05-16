function [out_Large]=sft(x,k)
% ERROR1 所有值平均的误差
% ERROR2 大值点的误差
%% 参数初始化
n=length(x);           % 信号长度
%k =2;                 % 信号稀疏度
repetitions = 1;	   % 外循环运行次数
Bcst_loc=8;            % loc分桶常量
Bcst_est=8;            % est分桶常量
Comb_cst=16;           % SFFT2.0混叠滤波器的分桶常量
loc_loops =3;		   % loc循环次数
est_loops =8;          % est循环次数
threshold_loops =2;	   % 大值的被标记次数阈值
Comb_loops = 1;        % SFFT2.0中混叠滤波的个数
simulate = 0;		   % 不运行完整个程序，仿真SFFT的运算时间和错误
snr=100000;            % 信噪比
std_noise = 0;         % 噪声放大值，std_noise = sqrt(k/(2*snr));
FFTW_OPT = false;      % 使用FFTW库
tolerance_loc = 1e-6;  % loc噪声（泄露）容限
tolerance_est = 1e-6;  % est噪声（泄露）容限

BB_loc = fix(Bcst_loc*sqrt(n*k/(log2(n))));
B_loc = floor_to_pow2(BB_loc);                 % loc所用B
BB_est = fix(Bcst_est*sqrt(n*k/(log2(n))));
B_est = floor_to_pow2(BB_est);                 % est所用B

B_thresh = 2*k;                                % 大桶的个数

lobefrac_loc = 0.5 / BB_loc;				   % loc通带区域[-lobefrac_loc,lobefrac_loc]
lobefrac_est = 0.5 / BB_est;                   % est通带区域[-lobefrac_est,lobefrac_est]

b_loc = fix(1.2*1.1*(n/BB_loc));
b_est = fix(1.4*1.1*(n/BB_est));

W_Comb = floor_to_pow2(Comb_cst*n/B_loc);      % 混叠滤波器所用B
WITH_COMB = false;                             % 运行SFTT2.0
ALGORITHM1 = true;                             % loc_loops=0;仅使用混叠滤波器

% word=sprintf('\n\nRUNNING EXPERIMENT: n=%d, k=%d.\n', n, k);
% disp(word);                                  % 输出n,k
% word=sprintf('\n\nSimulation:\n');
% disp(word);
% word=sprintf('*******************************************************\n');
% disp(word);

%% SFFT估计频率
%% 大桶数不超过桶数，标记次数不超过循环次数
if B_thresh > B_loc ||threshold_loops >loc_loops 
    return;
end
%% 构造平坦窗滤波器
[filtert,w_loc ]=make_dolphchebyshev_t(lobefrac_loc, tolerance_loc);                 % 构造 filtert 切比雪夫滤波器，w_loc：窗长
[filter_time filter_sizet filter_freq]=make_multiple_t(filtert, w_loc, n, b_loc);  % 构造平坦窗滤波器，filter_time：滤波器时域 filter_sizet：滤波器窗长 filter_freq：滤波器频域
% figure;plot(abs(filter_time))
% figure;plot(abs(filter_freq))
% aaaaa=fft(filtert);
% figure;plot(abs(aaaaa(1:8000)));title('dolphchebyshev filter time domain - window function');
% figure;plot(abs(filter_time));title('dolphchebyshev filter time domain - flat window');
% figure;plot(abs(filter_freq(1:8000)));title('dolphchebyshev filter freq domain - flat window');
[filtert_est w_est]=make_dolphchebyshev_t(lobefrac_est, tolerance_est);
[filter_est_time filter_est_sizet filter_est_freq]=make_multiple_t(filtert_est, w_est, n, b_est);
% word=sprintf(' Window size: Location Filter : %d; Estimation Filter : %d;\n', w_loc, w_est);
% disp(word);

filter_noise = 0;                       % 平坦窗滤波器的泄露
filter_noise_est = 0;
for i=0:9
    filter_noise = max(filter_noise,max(abs(filter_freq(n/2+i+1)),abs(filter_freq(n/2-i+1))));
    filter_noise_est = max(filter_noise_est,max(abs(filter_est_freq(n/2+i+1)),abs(filter_est_freq(n/2-i+1))));
end
% word=sprintf('Noise in filter: Location Filter : %d; Estimation Filter %d\n', filter_noise, filter_noise_est);
% disp(word);
% word=sprintf('****************************************************************\n\n');
% disp(word);
% word=sprintf('sFFT Results\n**************************************************************');
% disp(word);
%% 外循环
loops=loc_loops+est_loops;          % 总的内循环次数
for ii=1:repetitions
    % 外循环，每次得到heavy值的位置和大小
    tic;
    [I,out_I,out]= outer_loop(x, n,filter_time,filter_sizet, filter_freq, filter_est_time,...
        filter_est_sizet,filter_est_freq, B_est, B_thresh, B_loc, W_Comb, ...
        Comb_loops, threshold_loops, loc_loops, loops , WITH_COMB, ALGORITHM1);
    toc
    % I     大值的位置（升序）
    % out_I 对应的估值
    % out   hits位置集对应的估值
end
%%
num_candidates= length(out_I);      % 所估计的大值个数
x_f_Large = zeros(1,n);                 
out_Large = zeros(1,n);
candidates=zeros(2,num_candidates);
counter=0;
%% 用估计结果构造信号的频域
for i=1:num_candidates                    % 将算法估计的结果保存在二维矩阵，第二行是位置，第一行是相应的幅值
    counter=counter+1;
    candidates(1,counter)=abs(out_I(i));  % 大值的幅度序列
    candidates(2,counter)=I(i);           % 大值的位置序列
end
temp=sortrows((candidates).').';          % 大值按幅值升序排序
for i=1:k                                 % 取K个最大的幅值，组成信号频域out_Large
    key=temp(2,num_candidates - k+ i);
    out_Large(key+1) =out(key+1)*n;       
end

end


%% 输入的数据是2的n次幂
function [out] = floor_to_pow2 (in)
    out = 1;
    while out <= in
        out = out*2;
    end
    out=out/2;
end


%% 构造切比雪夫滤波器
function [ filter ,w] = make_dolphchebyshev_t(lobefrac,tolerance)
%% 参数
% lobefrac 通带截频
% tolerance 泄露容限
%% 窗长
w =fix((1 / pi) * (1/lobefrac) * acosh(1/tolerance)); %  窗长  fix():向0取整
if ~mod(w,2)==1                                       % ~：非
    w=w-1;						                      % 保证w是奇数，滤波器是对称的
end
%% 滤波器
filter = zeros(1,w);
t0 = cosh(acosh(1/tolerance) / (w-1));
for i=0:(w-1)
    filter(i+1) = Cheb(w-1, t0*cos(pi*i/w))*tolerance; % 切比雪夫滤波器的时域
end
%figure;plot(abs(filter));title('win');
temp=fft(filter, w);			                     % 切比雪夫滤波器的频域
%figure;plot(real(temp));title('FREwin');
filter=fftshift(temp);                               % 将零频点移到频谱的中间（将fft处理之后的pi-2pi部分搬移至-pi-0，从而使零频分量居于频谱的中心位置。）
%figure;plot(real(temp));title('FREwin_half');

for i=1:w
    filter(i) = real(filter(i));                	 % 返回的是实部
end
end

%% 利用切比雪夫窗函数生成的平坦窗函数
function [out_time,out_sizet,out_freq]=make_multiple_t(filtert, w, n, b)
%% 参数
% filtert 切比雪夫滤波器（频域）
% w 滤波器的窗长
% n 信号长度
% b 滤波器个数
% out_time  平坦窗滤波器（时域）
% out_sizet 平坦窗滤波器窗长
% out_freq  平坦窗滤波器（频域）
%%  滤波器的窗长和数量均不超过n
if b >= n || w >= n
    return;
end
%% 
g=zeros(1,n);
h=zeros(1,n);
w1=fix(w/2);
g(1:w-w1)=filtert((w1+1):w);
g(n-w1+1:n)=filtert(1:w1);
% figure;plot(abs(g));title('win');
g=fft(g,n);
s = 0;
for i=1:b
    s = g(i)+s ;									% 滤波器能量
end
offset = fix(b/2);
for i=0:(n-1)
    h(mod((i+offset),n)+1) = s;
    s = s + (g(mod((i+ b),n)+1) - g(i+1));
end

h_max=max(abs(h));
for i=1:n
    h(i) =h(i)/ h_max;
end
offsetc = 1;
step=exp(-2j*pi*w1/n);
for i=1:n
    h(i)= h(i)*offsetc;
    offsetc= offsetc *step;							% 在频域上对每个幅值进行移位
end
g=ifft(h,n);
filtert(1:w)=g(1:w);

out_time=filtert;								    % 时域信号
out_sizet = w;										% 滤波器频宽
out_freq= h;										% 频域信号
end

%% 稀疏傅立叶变换的外循环
function [I,out_I,out]=outer_loop(origx, n, filter_time,filter_sizet, filter_freq, ...
    filter_est_time,filter_est_sizet, filter_est_freq, B2,num,B,...
    W_Comb, Comb_loops, loop_threshold, location_loops, loops,  WITH_COMB, ALGORITHM1)
%% 参数
% origx 原始信号
%  n    信号长度
% filter_time 平坦窗滤波器（时域）
% filter_sizet 平坦窗滤波器窗长
% filter_freq 平坦窗滤波器（频域）
% B  loc的B 
% B2  est的B
% num 大桶的个数
% W_Comb  混叠滤波器所用B
% Comb_loops  SFFT2.0中混叠滤波的个数
% location_loops  loc循环次数
% loop_threshold 大值的被标记次数阈值
% loops   内循环总次数
% WITH_COMB 运行SFTT2.0
% ALGORITHM1  仅使用混叠滤波器，loc_loops=0;
% LARGE_FREQ
%output:
% I     大值的位置（升序）
% out_I 对应的估值
% out   估值关于位置的函数

%%
permutea =zeros(1,loops);
permuteb =zeros(1,loops);
x_samp=zeros(loops,max(B,B2));
score =zeros(1,n);
hits_ini=zeros(1,n);
hits_found_ini=0;
% a_INI=[3181497 3978157 3467391 2445571 1072031 2119695 1002225 3676519 3995001 ...
%     40375 1777721 610525 2569107 1769277 4168711];
% a_INI=[2029 10759 4377  12349  15373  4097  1123  8219  3823  13727  15029];

for i=1:loops
    %% 设定参数
    a = 0;                                  % 随机重排参数sigma的模逆
    b = 0;							        % 重排时域下标平移参数 （相当于课本上的-a*sigma）
    while(undivide(a, n)~= 1 || mod(a,2)==0 || mod(a,5)==0)
        a = mod(fix(rand()*n) ,n);			% a是n的质数，结束循环
    end
%     a
%     a=a_INI(i);
    ai = mod_inverse(a, n);			        % 模逆 a*ai mod n = 1
    permutea(i)=ai;						    % 保存每次循环的随机重排参数ai
    permuteb(i) = b;					    % 保存每次循环的随机重排参数 b
 
    perform_location =(location_loops>=i);   
    if perform_location                     % 设置定位循环窗函数
        cur_filter_time=filter_time;
        cur_filter_sizet=filter_sizet;
        cur_filter_freq=filter_freq;
        cur_B=B;
    else                                     % 设置估值循环窗函数
        cur_filter_time=filter_est_time;
        cur_filter_sizet=filter_est_sizet;
        cur_filter_freq=filter_est_freq;
        cur_B=B2;
    end
   
    %% 分桶：
    [x_samp(i,1:cur_B),J ]=inner_loop_locate(origx,n,cur_filter_time,cur_filter_sizet,cur_filter_freq, ...
        num, cur_B,	a,ai, b);              
    % x_samp：降采样后的频率值，桶值。
    % J：大桶的下标集
    % figure;plot(abs(fft(x_samp(i,1:cur_B))))；
    %%  累计大值原始位置集：
     if  perform_location   
        [hits,hits_found,score_out ]=inner_loop_filter_regular(J, n, num, cur_B, a, ai, b, loop_threshold, score,hits_ini, hits_found_ini);
        hits_ini=hits;                  % 找到的所有的大值位置集
        hits_found_ini=hits_found;      % 找到的大值位置总数
        score=score_out;                % 标记分数集 
    end
end
% word=sprintf('Number of candidates: %d\n', hits_found);
% disp(word);
%% 估值
[I,out_I,out] = estimate_values(hits, hits_found, x_samp,  loops, n, permutea, ...
    B,B2,filter_time,filter_sizet, filter_freq, filter_est_time, ...
    filter_est_sizet, filter_est_freq,  location_loops);

end

%% 切比雪夫窗函数
function [out]= Cheb ( m, x)
if abs(x) <= 1
    out=cos(m * acos(x)); % acos(x):反余弦arccos(x) 
else
    out=real(cosh(m * acosh(x))); % cosh(x):双曲余弦函数=(e^x+e^-x)/2
end
end

function [I,out_I,out] =estimate_values(hits, hits_found,x_samp,loops, n, ...
    permutea, B, B2,	filter_time,filter_sizet, filter_freq, filter_est_time, ...
    filter_est_sizet, filter_est_freq, location_loops)
%% 参数
% hits  大值频率位置集
% hits_found 大值频率个数
% x_samp 降采样后的频率值，桶值。
% loops 总循环次数
% n 信号长度
% permutea 保存每次循环的随机重排参数ai
% B loc的B 
% B2  est的B
% filter_time 平坦窗滤波器（时域）
% filter_sizet 平坦窗滤波器窗长
% filter_freq 平坦窗滤波器（频域）
% location_loops loc循环次数
%output:
% I     大值的位置（升序）
% out_I 对应的估值
% out   估值关于位置的函数
%% 估值
out=zeros(1,hits_found);
I=zeros(1,hits_found);
out_val=zeros(2,loops);                      % 记录内循环的估值结果，第一行real，第二行imag
location = fix((loops - 1) / 2);                                  
for i=0:(hits_found-1)                       % 第i个大值的位置
    r_loops=0; 
    for j=0:(loops-1)                        % 对第j次内循环估值
        %% 平坦窗滤波器参数设置
        if j < location_loops                                      
            cur_filter_time=filter_time;
            cur_filter_sizet=filter_sizet;
            cur_filter_freq=filter_freq;
            cur_B=B;
        else
            cur_filter_time=filter_est_time;
            cur_filter_sizet=filter_est_sizet;
            cur_filter_freq=filter_est_freq;
            cur_B=B2;
        end
        %% 分桶结果：桶值、哈希和偏移函数
        permuted_index= timesmod(permutea(j+1), hits(i+1), n);		% 重排后的下标
        hashed_to = fix(permuted_index / (n / cur_B));				% hash function 分到的桶下标
        dist = fix(rem(permuted_index , (n / cur_B)));				% offest function 桶内偏移量
        if dist > ((n/cur_B)/2)                                     % 若偏移量超过滤波器的通带截频，要分到下一个桶中
            hashed_to = rem((hashed_to + 1),cur_B);
            dist =dist- n/cur_B;
        end
        %% 估计频率值
        dist = rem((n - dist) , n);              % 分桶时对应的滤波器下标
        filter_value = cur_filter_freq(dist+1);  % 分桶时的滤波器值
%         x_samp(j+1,hashed_to+1) 
        out_val(1,r_loops+1) =real(x_samp(j+1,hashed_to+1) / filter_value); % 频率值实部
        out_val(2,r_loops+1) =imag(x_samp(j+1,hashed_to+1) / filter_value); % 频率值虚部
        r_loops=r_loops+1;     
    end
    %% 在所有估值中取中值
    for ii=1:2
        out_val(ii,:)=sort(out_val(ii,1:r_loops)); % 将位置i的所有估值按升序排序
    end
    realv = out_val(1,location+1);                 % 实部取中值
    imagv = out_val(2,location+1);                 % 虚部取中值
    out(hits(i+1)+1) = realv + 1j*imagv;           % 第i个大值最终的频率估值
    if hits(i+1)~=0
        I(i+1)=hits(i+1);                          % 第i个大值的位置                   
    end
end
I=sort(I);                                         % 将位置按升序排序
out_I=out(I+1);                                    % 排序后对应的频率值
I=[0,I];
out_I=[0,out_I];

end

%%   找到大桶，将下标保存在out（升序）
function [out]=find_largest_indices( num, samples,cur_B)
%% 参数
% samples   功率谱
% num   大桶个数 
% cur_B  分桶个数
% out 大值频域可能存在的坐标集
%% 
count = 0;
% out=zero(1,num);
cutoff = nth_element_immutable(samples, cur_B, cur_B-num);  %大桶的功率阈值
% 先找严格大于阈值的大桶
for i=1:cur_B
    if samples(i) > cutoff      
        count=count+1;
        out(count) = i-1;
    end
end
% 确保找到num个大桶
if count < num						% 若找到的点数小于B_thresh（即num），可以再增加一些。
    for i=1:cur_B
        if samples(i) == cutoff            
            count=count+1;
            out(count) = i-1;		% 对应图2中的阴影部分
            if count >= num 
                break;
            end
        end
    end
    out=sort(out);                  % 按照升序排列
end
end

function [x_samp ]=inner_loop_est(x,filter_est_time,filter_est_sizet, ...
   cur_B,cur_i,cur_index,x_sampt)

for i=1:filter_est_sizet
    x_sampt(cur_i(i)) = x(cur_index(i))*filter_est_time(i)+x_sampt(cur_i(i));			%滤波器在时域中体现为点乘形式，对应z的时域值
%     x_sampt(cur_i(i)) = x(cur_index(i))+x_sampt(cur_i(i));			%滤波器在时域中体现为点乘形式，对应z的时域值
end

x_samp=fft(x_sampt,cur_B);


end

function [hits, hits_found,score_out]=inner_loop_filter_regular(J, n, num, cur_B, a, ...
    ai, b,loop_threshold, score,hits_ini, hits_found_ini)
%% 参数
 % J 大桶下标集
 % n 信号长度
 % num 大桶数量
 % cur_B 分桶数
 % a ai的模逆
 % ai 重排时域下标伸缩参数（相当于课本上的sigma）
 % b 重排时域下标平移参数 （相当于课本上的-a*sigma）
 % loop_threshold 大值标记次数阈值
 % score  标记次数集
 % hits_ini %已找到的大值位置集
 % hits_found_ini 已找到的大值位置个数
%output
 % hits 找到的所有的大值位置集
 % hits_found 找到的大值位置总数
 % score_out 标记分数集 
 
%%
hits=hits_ini;
hits_found=hits_found_ini;
for i=1:num         
    % 找分到大桶J(i)的原始位置集 ，进行标记
    low  = mod(ceil((J(i)-0.5)*n/cur_B + n), n);	    % ceil函数:朝正无穷大方向取整 
    high = mod(ceil((J(i)+0.5)*n/cur_B + n), n);		% 分到桶J(i)的频率位置集[low,high]
    loc = timesmod(low, a, n);                          % 重排后位置low 对应原始信号的位置loc=a*low-a*b
    j = low;
    while  j ~= high                                    % 对分到桶J(i)的原始频率位置标记一分
        score(loc+1)=score(loc+1)+1;                   
        if  score(loc+1)==loop_threshold                % 若被标记次数达到阈值，将该位置存在hits
            hits_found =hits_found+1;
            hits(hits_found)=loc;
        end
        loc = mod((loc + a),n);
        j = mod((j + 1),n);
    end
end

score_out=score;

end


%% 定位内循环
function [x_samp,J ]=inner_loop_locate(origx, n, cur_filter_time,cur_filter_sizet, ...
    cur_filter_freq,num,cur_B,  a, ai,b)
%% 参数
% origx 原始信号
% n 信号长度
% cur_filter_time 平坦窗滤波器（时域）
% cur_filter_sizet 平坦窗滤波器窗长
% cur_filter_freq 平坦窗滤波器（频域）
% num 大桶个数
% cur_B 分桶数
% a ai的模逆
% ai 重排时域下标伸缩参数（相当于课本上的sigma）
% b 重排时域下标平移参数 （相当于课本上的-a*sigma）
% x_samp ：桶值
% J ：num个大桶坐标集
%%  要求信号长度是分桶数的整数倍，否则输出Warning
if mod(n,cur_B)       
    word=sprintf('Warning: n is not divisible by cur_B, which algorithm expects.\n');
    disp(word);
end
%% 分桶操作
x_sampt =zeros(1,n);  
index=b;
for i=0:(cur_filter_sizet-1)
    x_sampt(mod(i,cur_B)+1) = origx(index+1) * cur_filter_time(i+1) + x_sampt(mod(i,cur_B)+1); 
    index =rem((index+ai),n);	
    % origx(index+1)：随机重排，重排坐标index = b+ai*i,故 ai=sigma,b=-a*sigma
    % origx(index+1) * cur_filter_time(i+1)：滤波（时域相乘，频域卷积）
    % x_sampt(mod(i,cur_B)+1) =....： 频率降采样（对应的时域值）
end
x_samp=fft(x_sampt,cur_B);                    % 降采样频域上的采样值

%% 寻找大桶位置集
samples = zeros(1,cur_B);                     % 降采样后的功率谱
for i=1:cur_B
    samples(i) = (abs(x_samp(i)))^2;							
end
% figure;plot(samples);
J=find_largest_indices( num, samples ,cur_B);  % num个大桶的坐标
end



%% 求模逆 a*out mod n = 1
function [out]=mod_inverse(a, n)
i = n;
out = 0;
d = 1;
while a>0
    t = fix(i/a);
    x = a;
    a = mod(i,x);
    i = x;
    x = d;
    d = out - t*x;
    out = x;
end
out=rem(out,n);         % rem(x,y):求整除x/y的余数
if (out<0)
    out = rem((out+n),n);  
end
end


%% 输入离散信号的第n个元素
function [out,coordinate]=nth_element(a,b)
len=length(a);
for i=1:len
    swap=0;
    for j=1:len-1
        if a(j)>a(j+1)
            temp1=a(j);
            a(j)=a(j+1);
            a(j+1)=temp1;
            temp2=b(j);
            b(j)=b(j+1);
            b(j+1)=temp2;
            swap=1;
        end
    end
    if ~swap
        break
    end
end
out=a;
coordinate=b;

end


%% 获取大桶的功率阈值
function [out]= nth_element_immutable(input,  cur_B, num1)
% input 降采样后的功率谱
% cur_B  分桶数
% num1   非大值桶的数目
temp=input;
temp=sort(temp); % sort(X)：对X的元素进行升序排序；
out=temp(num1);  % 大桶的最小功率值
end



function [out]=undivide(a,b)

if mod(a,b)==0
    out=b;
else
    temp=undivide(b, mod(a,b));
    out=temp;
end

end

%% 重排前位置out,重排后位置x
function [out]=timesmod(x, a, n) 
out=fix(rem(x*a,n));      % rem:取余函数；fix:向0靠近取整;
end
