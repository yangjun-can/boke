%% 累计大值位置集（每次循环所得位置的并集）
function [hits, hits_found,score_out]=inner_loop_filter_regular(J, n, num,...
        cur_B, a,loop_threshold, score,hits_ini, hits_found_ini)
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
%     % 找分到大桶J(i)的原始位置集 ，进行标记
%     low  = mod(ceil((J(1:num)-0.5)*n/cur_B + n), n);	    % ceil函数:朝正无穷大方向取整 
%     high = mod(ceil((J(1:num)+0.5)*n/cur_B + n), n);		% 分到桶J(i)的频率位置集[low,high]
%     loc = timesmod(low, a, n);                          % 重排后位置low 对应原始信号的位置loc=a*low-a*b
%     len = mod(high-low,n);
% %     for i=1:num
% %         score(sub2ind(size(score),[mod(loc(i)+(0:len)*a,n)+1]))=score(sub2ind(size(score),[mod(loc(i)+(0:len)*a,n)+1]))+1;      
% %     end        
%   for i=1:num  
%     j = low(i);
%     while  j ~= high(i)                                    % 对分到桶J(i)的原始频率位置标记一分
%         score(loc(i)+1)=score(loc(i)+1)+1;                   
%         if  score(loc(i)+1)==loop_threshold                % 若被标记次数达到阈值，将该位置存在hits
%             hits_found =hits_found+1;
%             hits(hits_found)=loc(i);
%         end
%         loc(i) = mod((loc(i)+ a),n);
%         j = mod((j + 1),n);
%     end
%   end

score_out=score;

end


