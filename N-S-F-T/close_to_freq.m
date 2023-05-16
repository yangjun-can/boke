function [Int_Offset,Interpolate_Set]=close_to_freq(point,K,n,J)
 % point 要求的频域采样点
 % K 上采样点数 
 % n 要求的采样点数
 % J 上采样点中用于插值获取point点频率的
 
 residual=abs(n/K*(0:K-1)-point);%与point的距离
 [~,I]=sort (residual);%
 Interpolate_Set=I(1:J)-1;% 距离最近的J个点
 Int_Offset=I(1)-1-(J+1)/2; % 偏移
end