%% 大值点的估计操作
function [I,out_I,out] =estimate_values(hits, hits_found,x_samp,loops, n, ...
    permutea, B, B2,filter_freq, filter_est_freq, location_loops)
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
location = fix((loops - 1) / 2); %中值的位置
% out=zeros(1,hits_found);
% I=zeros(1,hits_found);
% out_val=zeros(2,loops);                      % 记录内循环的估值结果，第一行real，第二行imag
% for i=0:(hits_found-1)                       % 第i个大值的位置
%     r_loops=0;
%     for j=0:(loops-1)                        % 对第j次内循环估值
%         %% 平坦窗滤波器参数设置
%         if j < location_loops
%             %             cur_filter_time=filter_time;
%             %             cur_filter_sizet=filter_sizet;
%             cur_filter_freq=filter_freq;
%             cur_B=B;
%         else
%             %             cur_filter_time=filter_est_time;
%             %             cur_filter_sizet=filter_est_sizet;
%             cur_filter_freq=filter_est_freq;
%             cur_B=B2;
%         end
%         %% 分桶结果：桶值、哈希和偏移函数
%         %      permuted_index= timesmod(permutea(j+1), hits(i+1), n);		% 重排后的下标
%         permuted_index= mod(permutea(j+1)* hits(i+1), n);
%         hashed_to = fix(permuted_index / (n / cur_B));				% hash function 分到的桶下标
%         dist = fix(rem(permuted_index , (n / cur_B)));				% offest function 桶内偏移量
%         if dist > ((n/cur_B)/2)                                     % 若偏移量超过滤波器的通带截频，要分到下一个桶中
%             hashed_to = rem((hashed_to + 1),cur_B);
%             dist =dist- n/cur_B;
%         end
%         %% 估计频率值
%         dist = rem((n - dist) , n);              % 分桶时对应的滤波器下标
%         filter_value = cur_filter_freq(dist+1);  % 分桶时的滤波器值
%         %         x_samp(j+1,hashed_to+1)
%         out_val(1,r_loops+1) =real(x_samp(j+1,hashed_to+1) / filter_value); % 频率值实部
%         out_val(2,r_loops+1) =imag(x_samp(j+1,hashed_to+1) / filter_value); % 频率值虚部
%         r_loops=r_loops+1;
%     end
%     %% 在所有估值中取中值
%     out_val=sort(out_val(:,1:r_loops),2); % 将位置i的所有估值按升序排序
%     realv = out_val(1,location+1);                 % 实部取中值
%     imagv = out_val(2,location+1);                 % 虚部取中值
%     out(hits(i+1)+1) = realv + 1i*imagv;           % 第i个大值最终的频率估值
%     if hits(i+1)~=0
%         I(i+1)=hits(i+1);                          % 第i个大值的位置
%     end
% end
% I=sort(I);                                         % 将位置按升序排序
% out_I=out(I+1);                                    % 排序后对应的频率值
% I=[0,I];
% out_I=[0,out_I];

%% 平坦窗滤波器参数设置
cur_B=[repmat(B,hits_found,location_loops),repmat(B2,hits_found,loops-location_loops)];
%% 分桶结果：桶值、哈希和偏移函数
permuted_index= mod(hits(1:hits_found).'*permutea , n);           % 重排后的下标
hashed_to = fix(permuted_index .* cur_B /n );				% hash function 分到的桶下标
dist = fix(rem(permuted_index , (n ./ cur_B)));				% offest function 桶内偏移量
% 若偏移量超过滤波器的通带截频，要分到下一个桶中
pand = dist > (n./cur_B)/2 ;
hashed_to(pand) = rem((hashed_to(pand) + 1),cur_B(pand));
dist(pand) =dist(pand)- n./cur_B(pand);
%% 估计频率值
dist = rem((n - dist) , n);              % 分桶时对应的滤波器下标
dist1=dist(:,1:location_loops);
dist2=dist(:,location_loops+1:loops);
filter_value_1=[filter_freq(dist1+1).';filter_est_freq(dist2+1).'].';
% out_valy(1,j+1) =real(x_samp(j+1,hashed_to(i+1,j+1)+1) / filter_value(i+1,j+1)); % 频率值实部
out_val=zeros(2,hits_found,loops);
filter_value=reshape(filter_value_1,1,[]);
hashed_to1=reshape(hashed_to,1,[]);
out_val1=real(x_samp(:,hashed_to1+1) ./ filter_value); % 频率值实部
out_val2=imag(x_samp(:,hashed_to1+1) ./ filter_value); % 频率值虚部
for jj=1:loops
    out_val(1,:,jj)=out_val1(jj,(jj-1)*hits_found+1:jj*hits_found);
    out_val(2,:,jj)=out_val2(jj,(jj-1)*hits_found+1:jj*hits_found);
end
%% 在所有估值中取中值
out_val=sort(out_val,3); % 将位置i的所有估值按升序排序
realv = out_val(1,:,location+1);                 % 实部取中值
imagv = out_val(2,:,location+1);                 % 虚部取中值
out(hits((1:hits_found))+1) = realv + 1i*imagv;  % 第i个大值最终的频率估值
ton= hits~=0;  
I(ton)=hits(ton);                                % 第i个大值的位置
I=sort(I);                                       % 将位置按升序排序
out_I=out(I+1);                                  % 排序后对应的频率值
I=[0,I];
out_I=[0,out_I];
end



