 %% 频域非均匀采样点
function [kesai]=Nonuniform_sampling_point(m,n)
% m:第m个采样点[0,n-1]
% n:总采样点数
% kesai：采样点的位置
  m = mod(m,n);
%   kesai = m + 0.1 * rand(1);
kesai = m + 0.01 ;
% if m <= n/2 && m>=0
%     kesai=m/2;
% else 
%     if m>n/2 && m<3*n/4
%         kesai=m-n/4;
%     else
%         kesai= 2*m+1-n;
%     end
% end
end