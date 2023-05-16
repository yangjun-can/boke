%%  ���
function [ERROR1,ERROR2]=run_error(x, n, k,LARGE_FREQ, out_Large )
%% �����б�
% x �ź�
% n �źų���
% k �ź�ϡ���
% LARGE_FREQ ��ʵ��ֵƵ��λ��
% x_f��fft(x)
% out_Large���㷨���õ��ź�Ƶ�� 
% ERROR1 ����ֵƽ�������
% ERROR2 ��ֵ������
%% ʵ��ϡ��Ƶ��
x_f_Large = zeros(1,n);            
x_f=fft(x);
for i=1:k
    x_f_Large(LARGE_FREQ(i)+1)=x_f(LARGE_FREQ(i)+1);  % ʵ�ʵĴ�ֵλ�ú�ֵ
end
 %% �ҵ��Ĵ�ֵ����
% large_found = 0;
% FOUND=0;
% for i=1:k
%     if out(LARGE_FREQ(i)+1)==0
%         FOUND= FOUND +0;
%     else
%         FOUND= FOUND +1;
%     end
%     large_found = large_found +(out_Large((LARGE_FREQ(i)+1)) ~= 0);
% end
%% ������ƽ������
ERROR =0;
for i=1:n
    ERROR = ERROR +abs(out_Large(i)- x_f_Large(i));  % �������ܷ���
end
ERROR1=ERROR/k;   % ����ƽ������
%% ƽ��������Ȳ�
ERROR2 =0;
for i=1:k 
    ERROR2 = ERROR2 +abs(abs(out_Large(LARGE_FREQ(i)+1))- abs(x_f_Large(LARGE_FREQ(i)+1))); %������ȵ�����
end
ERROR2=ERROR2/k;
% word=sprintf('���ơ���ʵ��ֵ����\n ʵ�ʴ�ֵ����=%f ���ƴ�ֵ����=%f ����������ʧ)\n',abs(x_f_Large(LARGE_FREQ(i)+1)),abs(out_Large(LARGE_FREQ(i)+1)));
% disp(word);
% word=sprintf('***********************************************************************\n\n');

% word=sprintf('ERROR:\nK=%d; MISSED_Sample=%d; ��ֵ����ʧ���� ERROR= %f (%f ����������ʧ)\n',k, k-large_found, ERROR2, ERROR1);
% disp(word);
% word=sprintf('***********************************************************************\n\n');
% disp(word);
end


