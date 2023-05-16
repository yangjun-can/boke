%% �����б�ѩ���˲���
function [ filter ,w] = make_dolphchebyshev_t(lobefrac,tolerance)
%% ����
% lobefrac ͨ����Ƶ
% tolerance й¶����
%% ����
w =fix((1 / pi) * (1/lobefrac) * acosh(1/tolerance)); %  ����  fix():��0ȡ��
if ~mod(w,2)==1                                       % ~����
    w=w-1;						                      % ��֤w���������˲����ǶԳƵ�
end
%% �˲���
% filter = zeros(1,w);
t0 = cosh(acosh(1/tolerance) / (w-1));

% for i=0:(w-1)
%     filter(i+1) = Cheb(w-1, t0*cos(pi*i/w))*tolerance; % �б�ѩ���˲�����ʱ��
% end
 filter = Cheb(w-1, t0*cos(pi*(0:w-1)/w))*tolerance;

%figure;plot(abs(filter));title('win');
temp=fft(filter, w);			                     % �б�ѩ���˲�����Ƶ��
%figure;plot(real(temp));title('FREwin');
filter=fftshift(temp);                               % ����Ƶ���Ƶ�Ƶ�׵��м䣨��fft����֮���pi-2pi���ְ�����-pi-0���Ӷ�ʹ��Ƶ��������Ƶ�׵�����λ�á���
%figure;plot(real(temp));title('FREwin_half');

% for i=1:w
%     filter(i) = real(filter(i));                	 % ���ص���ʵ��
% end
filter = real(filter);
end