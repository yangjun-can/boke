%%   �ҵ���Ͱ�����±걣����out������
function [out]=find_largest_indices( num, samples,cur_B)
%% ����
% samples   ������
% num   ��Ͱ���� 
% cur_B  ��Ͱ����
% out ��ֵƵ����ܴ��ڵ����꼯
%% 
% count = 0;
out=zeros(1,num);
cutoff = nth_element_immutable(samples, cur_B-num);  %��Ͱ�Ĺ�����ֵ
%% �����ϸ������ֵ�Ĵ�Ͱ
% for i=1:cur_B
%     if samples(i) > cutoff      
%         count=count+1;
%         out(count) = i-1;
%     end
% end
lage=samples>cutoff;
count=size(samples(lage));
out(1:count) = find(lage)-1 ;

%% ȷ���ҵ�num����Ͱ
if count < num						% ���ҵ��ĵ���С��B_thresh����num��������������һЩ��
%     for i=1:cur_B
%         if samples(i) == cutoff            
%             count=count+1;
%             out(count) = i-1;		% ��Ӧͼ2�е���Ӱ����
%             if count >= num 
%                 break;
%             end
%         end
%     end
 out(count+1:num)= find(samples==cutoff,num-count) -1 ;
end
out=sort(out);                  % ������������
end