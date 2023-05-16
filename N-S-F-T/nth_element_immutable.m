%% 获取大桶的功率阈值
function [out]= nth_element_immutable(input, num1)
% input 降采样后的功率谱
% cur_B  分桶数
% num1   非大值桶的数目
% temp=input;
temp=sort(input); % sort(X)：对X的元素进行升序排序；
out=temp(num1);  % 大桶的最小功率值
end