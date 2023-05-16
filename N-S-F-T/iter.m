%% 第一次迭代
function [I,out_I,out]=iter(origx, n, filter_time,filter_sizet,  ...
    filter_freq, num,B,pe_major,inv_pe_major,sa_major,f)
% function [I,out_I,out]=nusfft_outer_loop(origx, n, filter_time,filter_sizet,  ...
%     filter_freq,filter_est_time,filter_est_sizet, filter_est_freq, B2,num,B,...
%     loop_threshold, location_loops, loops,sa_major,permutea,permuteb,permuteMateix,f)
%% 参数
% origx 原始信号
%  n    信号长度
% filter_time 平坦窗滤波器（时域）
% filter_sizet 平坦窗滤波器窗长
% filter_freq 平坦窗滤波器（频域）
% B  loc的B
% B2  est的B
% num 大桶的个数
% W_Comb  混叠滤波器所用B
% Comb_loops  SFFT2.0中混叠滤波的个数
% location_loops  loc循环次数
% loop_threshold 大值的被标记次数阈值
% loops   内循环总次数
%output:
% I     大值的位置（升序）
% out_I 对应的估值
% out   估值关于位置的函数
%%

%% 设定参数
loops=2; % 两次分桶，调制参数分别取0,1
tao=[0,1];
b = 0;							        % 重排时域下标平移参数 （相当于课本上的-a*sigma）
a = 0;                                  % 随机重排参数sigma的模逆
while(undivide(a, n)~= 1 || mod(a,2)==0 || mod(a,5)==0)
    a = mod(fix(rand()*n) ,n);			% a是n的质数，结束循环
end
ai = mod_inverse(a, n);			        % 模逆 a*ai mod n = 1
flat_filter_time=zeros(n,1,'single');
flat_filter_time(1:filter_sizet)=filter_time;
x_samp=zeros(loops,B,'single');         % 保存2次分桶的桶值
J=zeros(loops,num,'single');            % 保存2次分桶的大桶坐标
out=zeros(1,n);                       % 保存大值频率位置和值的函数

%% 两次分桶
tic;
for i=1:loops
    taoi=tao(i);
    x=origx(mod((0:n-1)-ai*taoi,n)+1); % 信号调制
    %% 分桶,并找到大桶：
    [x_samp(i,1:B),J(i,1:num)]=nusfft_inner_loop_locate(x,n,flat_filter_time, ...
        num, B,ai,b,pe_major,inv_pe_major,sa_major,f);
    % x_samp：降采样后的频率值，桶值。
    % J：大桶的下标集
end
%% 估计大值 
in = intersect(J(1,:), J(2,:));    % 找到2次分桶都是大桶的桶坐标
hs=x_samp(:,in);  % 大桶的两次桶值
% ahs = abs(hs);
% I=[];
% out_I=[];
% for i = 1:length(in)  % 对于大桶i 如果是1稀疏的，估计大频率
%     % 1-sparse detection
% %     if var(ahs(:,i))< 1e-5
%         u0 = wrapTo2Pi(angle(hs(1,i)/hs(2,i)))*n/(2*pi); % 大频率的坐标 % angle 相位角 wrapTo2Pi弧度角转换到[0,2pi]
%         I(end+1) = mod(round(u0*a), n);    % 储存大频率的坐标（消除噪声影响，整数） end+1:在当前的基础上增加一列
% %         off =mod(ai*(I(end)-b)-in(i)*n/B,n);  % 偏移量
% %         dist = rem((n - off) , n);              % 分桶时对应的滤波器下标
% %         out_I(end+1)=hs(1,i)/filter_freq(dist+1);  % 频率值约为桶值
%         out_I(end+1)=hs(1,i);
%         out(I+1)= out_I;
% %     end
% end
u0 = wrapTo2Pi(angle(hs(1,:)./hs(2,:)))*n/(2*pi); % 大频率的坐标 % angle 相位角 wrapTo2Pi弧度角转换到[0,2pi]
I = mod(round(u0*a), n) ;   % 储存大频率的坐标（消除噪声影响，整数） end+1:在当前的基础上增加一列
out_I =hs(1,:);
out(I+1)= out_I;
I=[0,I];
out_I=[0,out_I];
%% 计时
disp(['NUSFT的运行时间为：',num2str(toc),'秒']);
end