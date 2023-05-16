%% ��λ��ѭ�� (��Ͱ)
function [x_samp,J]=nusfft_inner_loop_locate(origx, n, flat_filter_time, ...
    num,cur_B,ai,b,pe_major,inv_pe_major, sa_major,f)

% function [x_samp,J]=nusfft_inner_loop_locate(origx, n, flat_filter_time, ...
%     num,cur_B,sa_major,permute,f)

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
% x2=zeros(n,1);
% index=b;
%% �������
x1=inv_pe_major*origx.';
% for i=0:(n-1)
%     x2(i+1)=x1(index+1);
%     index =rem((index+ai),n);
% end
% x2=x1(rem(b+(0:n-1)*ai,n)+1);
% x3=permute*origx.';
x3=pe_major*x1(rem(b+(0:n-1)*ai,n)+1);
%  figure;
%  subplot(2,1,1);plot(abs(x3));title('(a)');ylabel('amplitude');xlabel('Sampling point in time domion');
%  subplot(2,1,2);plot(f,abs(ndftld(x3.', n, f)));title('(b)');ylabel('amplitude');xlabel('Sampling point in frequency domion');
 %% ��ƽ̹���˲�
x4=x3.* flat_filter_time;
%% Ƶ���²���
% x5=sqrt(n/cur_B)*sa_major*x4.';
x5=n/cur_B*sa_major*x4;
% x5=x5.';
%% NUDFT
f2=Nonuniform_sampling_point(0:cur_B-1,cur_B);
f2=f2(:);
x_samp =chebfun.nufft(x5,f2/cur_B,2);
% x_samp = ndftld(x5,cur_B,f2).';
% x_samp = fft(x5);
%  figure;
%  subplot(2,1,1);plot(abs(x5));title('(a)');ylabel('amplitude');xlabel('Sampling point in time domion');
%  subplot(2,1,2);stem(f2,abs(x_samp));title('(b)');ylabel('amplitude');xlabel('Sampling point in frequency domion');
%% Ѱ�Ҵ�Ͱλ�ü�
% samples = zeros(1,cur_B);                     % ��������Ĺ�����
% for i=1:cur_B
%     samples(i) = (abs(x_samp(i)))^2;
% end
samples = (abs(x_samp)).^2;
J=find_largest_indices( num, samples ,cur_B);  % num����Ͱ������
end