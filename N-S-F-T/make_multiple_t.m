%% �����б�ѩ�򴰺������ɵ�ƽ̹������
function [out_time,out_sizet,out_freq]=make_multiple_t(filtert, w, n, b)
%% ����
% filtert �б�ѩ���˲�����Ƶ��
% w �˲����Ĵ���
% n �źų���
% b �˲�������
% out_time  ƽ̹���˲�����ʱ��
% out_sizet ƽ̹���˲�������
% out_freq  ƽ̹���˲�����Ƶ��
%%  �˲����Ĵ�����������������n
if b >= n || w >= n
    return;
end
%% 
g=zeros(1,n);
h=zeros(1,n);
w1=fix(w/2);
g(1:w-w1)=filtert((w1+1):w);
g(n-w1+1:n)=filtert(1:w1);
% figure;plot(abs(g));title('win');
g=fft(g,n);

% s = 0;
% for i=1:b
%     s = g(i)+s ;									% �˲�������
% end
s=sum(g(1:b));
offset = fix(b/2);
for i=0:(n-1)
    h(mod((i+offset),n)+1) = s;
    s = s + (g(mod((i+ b),n)+1) - g(i+1));
end

h_max=max(abs(h));
% for i=1:n
%     h(i) =h(i)/ h_max;
% end
h = h/ h_max;
step=exp(-2j*pi*w1/n);

% offsetc = 1;
% for i=1:n
%     h(i)= h(i)*offsetc;
%     offsetc= offsetc *step;							% ��Ƶ���϶�ÿ����ֵ������λ
% end
h=h.*step.^(0:n-1);
g=ifft(h,n);
filtert(1:w)=g(1:w);

out_time = filtert;								    % ʱ���ź�
out_sizet = w;										% �˲���Ƶ��
out_freq = h;										% Ƶ���ź�
end