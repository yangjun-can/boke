function [out_Large,BB_loc,BB_est,b_loc,b_est]=nusfft(x,k,f)
%% 参数初始化
% profile on;
x=single(x);
f=single(f);
n=length(x);           % 信号长度
%k                     % 信号稀疏度
%repetitions = 1;	   % 外循环运行次数
Bcst_loc=2;            % loc分桶常量
Bcst_est=2;            % est分桶常量
loc_loops =16;		   % loc循环次数
est_loops =6;          % est循环次数
tolerance_loc = 1e-4;     % loc噪声（泄露）容限 test中n=64,取1e-4； test中n>64,取1e-6；
tolerance_est = 1e-4;     % est噪声（泄露）容限
threshold_loops =fix(loc_loops/2);	   % 大值的被标记次数阈值

loops=loc_loops+est_loops;% 总的内循环次数
BB_loc = fix(Bcst_loc*sqrt(n*k/(log2(n))));
B_loc = floor_to_pow2(BB_loc);                 % loc所用B
BB_est = fix(Bcst_est*sqrt(n*k/(log2(n))));
B_est = floor_to_pow2(BB_est);                 % est所用B
B_thresh = 2*k;                                % 大桶的个数
lobefrac_loc = 0.5 / BB_loc;				   % loc通带区域[-lobefrac_loc,lobefrac_loc]
lobefrac_est = 0.5 / BB_est;                   % est通带区域[-lobefrac_est,lobefrac_est]
b_loc = fix(1.4*1.1*(n/BB_loc));
b_est = fix(1.2*1.1*(n/BB_est));
% 预计算重排矩阵
% permuteMateix=zeros(n,n,loops);
% sequence1=gpuArray(single(exp(-2*1i*pi*(0:n-1).*f(1:n)/n))); 
%  sequence1=gpuArray(exp(-2i*pi*(0:n-1).*f(1:n)/n));
 sequence1=exp(-2i*pi*(0:n-1).*f(1:n)/n); 
pe_major=n*ifft(sequence1);
inv_pe_major=inv(pe_major);
% 预计算下采样矩阵
% sequence2=gpuArray(single(exp(-2i*pi*(0:B_loc-1).*(f((1:B_loc)*n/B_loc)-(1:B_loc)'*(n/B_loc-1))/B_loc)));  % 预计算下采样矩阵
%  sequence2=gpuArray(exp(-2i*pi*(0:B_loc-1).*(f((1:B_loc)*n/B_loc)-(1:B_loc)'*(n/B_loc-1))/B_loc));
 sequence2=exp(-2i*pi*(0:B_loc-1).*(f((1:B_loc)*n/B_loc)-(1:B_loc)'*(n/B_loc-1))/B_loc);  % 预计算下采样矩阵
sa_major1=B_loc*ifft(sequence2); 
sa_major=repmat(sa_major1,[1,n/B_loc]);
sa_major=sa_major/pe_major;

% word=sprintf('\n\nRUNNING EXPERIMENT: n=%d, k=%d.\n', n, k);
% disp(word);                                  % 输出n,k
% word=sprintf('\n\nSimulation:\n');
% disp(word);
% word=sprintf('*******************************************************\n');
% disp(word);

%% NUSFFT估计频率
%% 大桶数不超过桶数，标记次数不超过循环次数
if B_thresh > B_loc ||threshold_loops >loc_loops
    return;
end
%% 构造平坦窗滤波器
[filtert,w_loc ]=make_dolphchebyshev_t(lobefrac_loc, tolerance_loc);                 % 构造 filtert 切比雪夫滤波器，w_loc：窗长
[filter_time, filter_sizet, filter_freq]=make_multiple_t(filtert, w_loc, n, b_loc);    % 构造平坦窗滤波器，filter_time：滤波器时域 filter_sizet：滤波器窗长 filter_freq：滤波器频域
filter_time=single(filter_time);
filter_freq=single(filter_freq);
% aaaaa=fft(filtert);
% figure;plot(abs(aaaaa));title('dolphchebyshev filter time domain - window function');
% figure;
% subplot(2,2,1);plot(abs(filter_time));title('(a)');ylabel('amplitude');
% subplot(2,2,3);plot(abs(filter_freq));title('(c)');ylabel('amplitude');
% subplot(2,2,2);semilogy(abs(filter_time));title('(b)');ylabel('amplitude');
% subplot(2,2,4);semilogy(abs(filter_freq));title('(d)');ylabel('amplitude');
[filtert_est, w_est]=make_dolphchebyshev_t(lobefrac_est, tolerance_est);
[filter_est_time ,filter_est_sizet ,filter_est_freq]=make_multiple_t(filtert_est, w_est, n, b_est);
filter_est_time=single(filter_est_time);
filter_est_freq=single(filter_est_freq);
% word=sprintf(' Window size: Location Filter : %d; Estimation Filter : %d;\n', w_loc, w_est);
% disp(word);

%% 平坦窗滤波器的泄露
% filter_noise = 0;
% filter_noise_est = 0;
% for i=0:9
%     filter_noise = max(filter_noise,max(abs(filter_freq(n/2+i+1)),abs(filter_freq(n/2-i+1))));
%     filter_noise_est = max(filter_noise_est,max(abs(filter_est_freq(n/2+i+1)),abs(filter_est_freq(n/2-i+1))));
% end
% word=sprintf('Noise in filter: Location Filter : %d; Estimation Filter %d\n', filter_noise, filter_noise_est);
% disp(word);
% word=sprintf('****************************************************************\n\n');
% disp(word);
% word=sprintf('sFFT Results\n**************************************************************');
% disp(word);
%% 外循环
% for ii=1:repetitions
    % 外循环，每次得到heavy值的位置和大小
        [I,out_I,out]= nusfft_outer_loop(x, n,filter_time,filter_sizet, filter_freq, filter_est_time,...
        filter_est_sizet,filter_est_freq, B_est, B_thresh, B_loc,threshold_loops, ...
        loc_loops,loops,pe_major,inv_pe_major,sa_major,f);
    
    % I     大值的位置（升序）
    % out_I 对应的估值
    % out   估值关于位置的函数
% end
%%
num_candidates= length(out_I);      % 所估计的大值个数
out_Large = zeros(1,n);
candidates=zeros(2,num_candidates);
% counter=0;
%% 用估计结果构造信号的频域
% for i=1:num_candidates                   
%     counter=counter+1;
%     candidates(1,counter)=abs(out_I(i));  % 大值的幅度序列
%     candidates(2,counter)=I(i);           % 大值的位置序列
% end 
% 将算法估计的结果保存在二维矩阵，第二行是位置，第一行是相应的幅值
candidates(1,:)=abs(out_I);  % 大值的幅度序列
candidates(2,:)=I;           % 大值的位置序列
temp=sortrows((candidates).').';          % 大值按幅值升序排序
% for i=1:k                               % 取K个最大的幅值，组成信号频域out_Large
%     key=temp(2,num_candidates - k+ i);
key=temp(2,num_candidates - k+ (1:k));
out_Large(key+1) =out(key+1)*n;
% end
% profile viewer;
end