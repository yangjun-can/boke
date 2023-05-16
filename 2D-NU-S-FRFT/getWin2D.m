function w = getWin2D(w1)
%%  将一维窗旋转
  N = length(w1);
  M = (N-1) / 2;
n = 2 / M * (-M:M);
[x,y] = meshgrid(n);
r = sqrt( x.^2 + y.^2 );
w = zeros(N);
w(:) = interp1(n, w1, r(:));
w(isnan(w)) = 0;
% n = length(w1);
%   if any(abs(w1-rot90(w1,2))>sqrt(eps))  %eps 计算机用于区分两个数的差的最小常数
%     msg = '1-D window must be symmetric to use Huang''s method.';
%     eid = sprintf('Images:%s:oneDWindowMustBeSymmetric',mfilename);
%     error(eid,msg);
%   end
%   if length(w1)<2
%     msg = 'Length of window must be greater than 1.';
%     eid = sprintf('Images:%s:windowLengthMustBeGreaterThanOne',mfilename);
%     error(eid,msg);
%   end
%   t = (-(n-1)/2:(n-1)/2)*(2/(n-1));
%   [t1,t2] = meshgrid(t,t);%形成格点矩阵
%   r = sqrt(t1.*t1 + t2.*t2);
%   d = find(r<t(1) | r>t(length(w1)));
%   if ~isempty(d), r(d) = zeros(size(d)); end
%   w = zeros(size(r)); w(:) = interp1(t,w1,r(:),'linear');
% %   插值函数为interp1调用格式为：  yi= interp1(x,y,xi,'method')           
% %     其中x，y为插值点，yi为在被插值点xi处的插值结果；x,y为向量， 
% %     'method'表示采用的插值方法，MATLAB提供的插值方法有几种： 
% %         'nearest'是最邻近插值， 'linear'线性插值； 'spline'三次样条插值； 'pchip'立方插值．缺省时表示线性插值
% %     注意：所有的插值方法都要求x是单调的，并且xi不能够超过x的范围
%   if ~isempty(d), w(d) = zeros(size(d)); end 
end
