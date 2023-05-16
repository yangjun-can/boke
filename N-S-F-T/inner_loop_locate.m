%% ��λ��ѭ��
function [x_samp,J ]=inner_loop_locate(origx, n, cur_filter_time,cur_filter_sizet, ...
    num,cur_B, ai,b)
%% ����
% origx ԭʼ�ź�
% n �źų���
% cur_filter_time ƽ̹���˲�����ʱ��
% cur_filter_sizet ƽ̹���˲�������
% cur_filter_freq ƽ̹���˲�����Ƶ��
% num ��Ͱ����
% cur_B ��Ͱ��
% a ai��ģ��
% ai ����ʱ���±������������൱�ڿα��ϵ�sigma��
% b ����ʱ���±�ƽ�Ʋ��� ���൱�ڿα��ϵ�-a*sigma��
% x_samp ��Ͱֵ
% J ��num����Ͱ���꼯
%%  Ҫ���źų����Ƿ�Ͱ�������������������Warning
if mod(n,cur_B)       
   fprintf('Warning: n is not divisible by cur_B, which algorithm expects.\n'); 
end
%% ��Ͱ����
x_sampt =zeros(1,n,'single');  
index=b;
for i=0:(cur_filter_sizet-1)   
    x_sampt(mod(i,cur_B)+1) = origx(index+1) * cur_filter_time(i+1) + x_sampt(mod(i,cur_B)+1); 
    index =rem((index+ai),n);	
    % origx(index+1)��������ţ���������index = b+ai*i,�� ai=sigma,b=-a*sigma
    % origx(index+1) * cur_filter_time(i+1)���˲���ʱ����ˣ�Ƶ������
    % x_sampt(mod(i,cur_B)+1) =....�� Ƶ�ʽ���������Ӧ��ʱ��ֵ��
end
x_samp=fft(x_sampt,cur_B);                    % ������Ƶ���ϵĲ���ֵ

%% Ѱ�Ҵ�Ͱλ�ü�
samples = zeros(1,cur_B);                     % ��������Ĺ�����
for i=1:cur_B
    samples(i) = (abs(x_samp(i)))^2;							
end
% figure;plot(samples);
J=find_largest_indices( num, samples ,cur_B);  % num����Ͱ������
end