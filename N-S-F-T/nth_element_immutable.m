%% ��ȡ��Ͱ�Ĺ�����ֵ
function [out]= nth_element_immutable(input, num1)
% input ��������Ĺ�����
% cur_B  ��Ͱ��
% num1   �Ǵ�ֵͰ����Ŀ
% temp=input;
temp=sort(input); % sort(X)����X��Ԫ�ؽ�����������
out=temp(num1);  % ��Ͱ����С����ֵ
end