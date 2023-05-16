%%  误差
function [ERROR1,ERROR2]=run_error(x, n, k,LARGE_FREQ, out_Large )
%% 参数列表
% x 信号
% n 信号长度
% k 信号稀疏度
% LARGE_FREQ 真实大值频率位置
% x_f：fft(x)
% out_Large：算法所得的信号频率 
% ERROR1 所有值平均的误差
% ERROR2 大值点的误差
%% 实际稀疏频域
x_f_Large = zeros(1,n);            
x_f=fft(x);
for i=1:k
    x_f_Large(LARGE_FREQ(i)+1)=x_f(LARGE_FREQ(i)+1);  % 实际的大值位置和值
end
 %% 找到的大值数量
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
%% 误差函数的平均幅度
ERROR =0;
for i=1:n
    ERROR = ERROR +abs(out_Large(i)- x_f_Large(i));  % 误差函数的总幅度
end
ERROR1=ERROR/k;   % 误差的平均幅度
%% 平均各点幅度差
ERROR2 =0;
for i=1:k 
    ERROR2 = ERROR2 +abs(abs(out_Large(LARGE_FREQ(i)+1))- abs(x_f_Large(LARGE_FREQ(i)+1))); %各点幅度的误差和
end
ERROR2=ERROR2/k;
% word=sprintf('估计、真实大值幅度\n 实际大值幅度=%f 估计大值幅度=%f 算上噪声损失)\n',abs(x_f_Large(LARGE_FREQ(i)+1)),abs(out_Large(LARGE_FREQ(i)+1)));
% disp(word);
% word=sprintf('***********************************************************************\n\n');

% word=sprintf('ERROR:\nK=%d; MISSED_Sample=%d; 大值点损失幅度 ERROR= %f (%f 算上噪声损失)\n',k, k-large_found, ERROR2, ERROR1);
% disp(word);
% word=sprintf('***********************************************************************\n\n');
% disp(word);
end


