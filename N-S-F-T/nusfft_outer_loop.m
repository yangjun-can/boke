%% 外循环
function [I,out_I,out]=nusfft_outer_loop(origx, n, filter_time,filter_sizet,  ...
    filter_freq,filter_est_time,filter_est_sizet, filter_est_freq, B2,num,B,...
    loop_threshold, location_loops, loops,pe_major,inv_pe_major,sa_major,f)
% function [I,out_I,out]=nusfft_outer_loop(origx, n, filter_time,filter_sizet,  ...
%     filter_freq,filter_est_time,filter_est_sizet, filter_est_freq, B2,num,B,...
%     loop_threshold, location_loops, loops,sa_major,permutea,permuteb,permuteMateix,f)
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
%output:
% I     大值的位置（升序）
% out_I 对应的估值
% out   估值关于位置的函数

%%
tic;
x_samp=zeros(loops,max(B,B2),'single');          % 保存每次循环的桶值
permutea=zeros(1,loops);
permuteb=zeros(1,loops);
score =zeros(1,n,'single');
hits_ini=zeros(1,n,'single');
hits_found_ini=0;
%% 内循环
for i=1:loops
    %     a=permutea(i);
    %     b=permuteb(i);
    %     permute=permuteMateix(:,:,i);
    %% 设定参数
    a = 0;                                  % 随机重排参数sigma的模逆
    b = 0;							        % 重排时域下标平移参数 （相当于课本上的-a*sigma）
    while(undivide(a, n)~= 1 || mod(a,2)==0 || mod(a,5)==0)
        a = mod(fix(rand()*n) ,n);			% a是n的质数，结束循环
    end
    ai = mod_inverse(a, n);     	        % 模逆 a*ai mod n = 1
    permutea(i)=ai;						    % 保存每次循环的随机重排参数ai
    permuteb(i) = b;					    % 保存每次循环的随机重排参数 b
    
    perform_location =(location_loops>=i);
    if perform_location                     % 设置定位循环窗函数
        cur_filter_time=filter_time;
        cur_filter_sizet=filter_sizet;
        %         cur_filter_freq=filter_freq;
        cur_B=B;
    else                                     % 设置估值循环窗函数
        cur_filter_time=filter_est_time;
        cur_filter_sizet=filter_est_sizet;
        %         cur_filter_freq=filter_est_freq;
        cur_B=B2;
    end
    flat_filter_time=zeros(n,1,'single');
    flat_filter_time(1:cur_filter_sizet)=cur_filter_time;
    %% 分桶,并找到大桶：
    [x_samp(i,1:cur_B),J]=nusfft_inner_loop_locate(origx,n,flat_filter_time, ...
        num, cur_B,ai,b,pe_major,inv_pe_major,sa_major,f);
    %     [x_samp(i,1:cur_B),J]=nusfft_inner_loop_locate(origx,n,flat_filter_time, ...
    %         num, cur_B,sa_major,permute,f);
    
    % x_samp：降采样后的频率值，桶值。
    % J：大桶的下标集
    %%  累计大值原始位置集：
    if  perform_location
        [hits,hits_found,score_out ]=inner_loop_filter_regular(J, n, num, cur_B, a, loop_threshold, score,hits_ini, hits_found_ini);
        hits_ini=hits;                  % 找到的所有的大值位置集
        hits_found_ini=hits_found;      % 找到的大值位置总数
        score=score_out;                % 标记分数集
    end
end
% word=sprintf('Number of candidates: %d\n', hits_found);
% disp(word);
%% 估值
[I,out_I,out] = estimate_values(hits, hits_found, x_samp,  loops, n, permutea, ...
    B, B2, filter_freq, filter_est_freq,  location_loops);
%% 计时
% disp('NUSFFT的运行时间为：');
toc
end