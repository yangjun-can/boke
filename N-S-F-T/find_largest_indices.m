%%   找到大桶，将下标保存在out（升序）
function [out]=find_largest_indices( num, samples,cur_B)
%% 参数
% samples   功率谱
% num   大桶个数 
% cur_B  分桶个数
% out 大值频域可能存在的坐标集
%% 
% count = 0;
out=zeros(1,num);
cutoff = nth_element_immutable(samples, cur_B-num);  %大桶的功率阈值
%% 先找严格大于阈值的大桶
% for i=1:cur_B
%     if samples(i) > cutoff      
%         count=count+1;
%         out(count) = i-1;
%     end
% end
lage=samples>cutoff;
count=size(samples(lage));
out(1:count) = find(lage)-1 ;

%% 确保找到num个大桶
if count < num						% 若找到的点数小于B_thresh（即num），可以再增加一些。
%     for i=1:cur_B
%         if samples(i) == cutoff            
%             count=count+1;
%             out(count) = i-1;		% 对应图2中的阴影部分
%             if count >= num 
%                 break;
%             end
%         end
%     end
 out(count+1:num)= find(samples==cutoff,num-count) -1 ;
end
out=sort(out);                  % 按照升序排列
end