clear all;close all;
% 地震波频谱分析
filename_mine='sd12_signal_mine_3.txt'; % 数据文件
xx_mine=load(filename_mine);  % 提取数据
%% 原始数据预处理
% xxx_mine=xx_mine(1:8192);
xxx_mine=xx_mine(1650:3697);
%第一步，去均值
xxx_aver_mine(:,1)=xxx_mine(:,1)-mean(xxx_mine(:,1));
%第二步，去趋势
detrend_xxx_aver_mine(:,1)= detrend(xxx_aver_mine(:,1));
%第三步，emd
signal_mine(:,1)=detrend_xxx_aver_mine(:,1);
z_mine=signal_mine(:,1);   
%% 
Fs=100;  % 采样频率
N=length(z_mine);  % 信号长度
t=0:1/Fs:(N-1)/Fs; % 采样时间序列
dt=1/Fs; % 采样间隔
Xt_mine=z_mine.';
f=zeros(N,1); % 非均匀频率点
for i=1:N
    f(i)=Nonuniform_sampling_point(i-1,N);
end
 [NUSFFT_freq,BB_loc,BB_est,b_loc,b_est]=nusfft(Xt_mine,2,f); 
 Xf_mine=fft(Xt_mine);
% [ha, pos] =tight_subplot(2,1,[.1 .1],[.1 .08],[.1 .05]); % 控制图像的边距
% axes(ha(1));
figure;subplot(3,1,1);
plot((0:N-1)/Fs,Xt_mine);xlabel('Time/s');ylabel('Amplitude');
grid on; title('(a) Mine wave in the time domain');
% axes(ha(2));
subplot(3,1,2); plot((0:N-1)/(N*dt),abs(Xf_mine)*2/N);xlabel('Frequency/Hz');ylabel('Amplitude');
 grid on; xlim([0 20]);title('(b) Direct method');title('(b) Direct method');
subplot(3,1,3); plot(f/(N*dt),abs(NUSFFT_freq)*2/N);xlabel('Frequency/Hz');ylabel('Amplitude');
 grid on; xlim([0 20]);title('(c) NUSFT method');



%% 控制图像边距
function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1; 
    gap = [gap gap];
end
if numel(marg_w)==1; 
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1; 
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh; 

% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
end