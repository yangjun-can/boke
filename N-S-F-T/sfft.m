function [out_Large,BB_loc,BB_est,b_loc,b_est]=sfft(x,k)
%% 参数初始化
n=length(x);           % 信号长度
%k =2;                 % 信号稀疏度
repetitions = 1;	   % 外循环运行次数
Bcst_loc=2;            % loc分桶常量
Bcst_est=2;            % est分桶常量
loc_loops =3;		   % loc循环次数
est_loops =8;          % est循环次数
threshold_loops =fix(loc_loops/2);	   % 大值的被标记次数阈值
% simulate = 0;		      % 不运行完整个程序，仿真SFFT的运算时间和错误
% snr=100000;             % 信噪比
% std_noise = 0;          % 噪声放大值，std_noise = sqrt(k/(2*snr));
% FFTW_OPT = false;       % 使用FFTW库
tolerance_loc = 1e-6;     % loc噪声（泄露）容限
tolerance_est = 1e-6;     % est噪声（泄露）容限
loops=loc_loops+est_loops;% 总的内循环次数
BB_loc = fix(Bcst_loc*sqrt(n*k/(log2(n))));
B_loc = floor_to_pow2(BB_loc);                 % loc所用B
BB_est = fix(Bcst_est*sqrt(n*k/(log2(n))));
B_est = floor_to_pow2(BB_est);                 % est所用B
B_thresh = 2*k;                                % 大桶的个数

lobefrac_loc = 0.5 / BB_loc;				   % loc通带区域[-lobefrac_loc,lobefrac_loc]
lobefrac_est = 0.5 / BB_est;                   % est通带区域[-lobefrac_est,lobefrac_est]

b_loc = fix(1.2*1.1*(n/BB_loc));
b_est = fix(1.4*1.1*(n/BB_est));

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
% aaaaa=fft(filtert);
% figure;plot(abs(aaaaa(1:8000)));title('dolphchebyshev filter time domain - window function');
% figure;plot(abs(filter_time));title(' flat window dolphchebyshev filter time domain ');
% figure;plot(abs(filter_freq(1:8000)));title('flat window dolphchebyshev filter freq domain ');
[filtert_est w_est]=make_dolphchebyshev_t(lobefrac_est, tolerance_est);
[filter_est_time filter_est_sizet filter_est_freq]=make_multiple_t(filtert_est, w_est, n, b_est);
% word=sprintf(' Window size: Location Filter : %d; Estimation Filter : %d;\n', w_loc, w_est);
% disp(word);
 %% 平坦窗滤波器的泄露
filter_noise = 0;                      
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

for ii=1:repetitions
    % 外循环，每次得到heavy值的位置和大小
    [I,out_I,out]= outer_loop(x, n,filter_time,filter_sizet, filter_freq, filter_est_time,...
        filter_est_sizet,filter_est_freq, B_est, B_thresh, B_loc, ...
        threshold_loops, loc_loops, loops );
    % I     大值的位置（升序）
    % out_I 对应的估值
    % out   估值关于位置的函数
end
%%
num_candidates= length(out_I);      % 所估计的大值个数
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