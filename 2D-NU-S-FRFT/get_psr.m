function [PSR] = get_psr(af)
%% 提取峰值旁瓣比, af为频谱
% 频谱幅值
af = fftshift(abs(af));
% 找到所有的peaks
% pks = findpeaks((af)); % 一维信号的峰值
PeaksMap = imregionalmax(af);% 0-1矩阵，峰值的位置都被标识为真，其他位置均为假
PeakLocation = PeaksMap == 1; % 峰值的位置
PeaksValue  = af(PeakLocation); % 所有峰值
% 对pks进行去重排序(升序)
pks = unique(PeaksValue); 
% 返回psr
PSR = pks(end)/pks(end-1);
PSR = mag2db(PSR);%转化为dB
end