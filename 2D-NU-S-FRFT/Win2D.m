n=256;
att = 40; % PSR in dB
% a_max = 70;
% SNR_dB = 40;
% sigma_n = 1;
% SNR = 10^(SNR_dB/20);
% a_max = sqrt(sigma_n^2*SNR); % 最小幅度
% % win0 = taylorwin(n); %泰勒窗
% % win0 = boxcar(n); %矩形窗
% % win0 = hanning(n); % 汉宁窗
% win0 = blackman(n); % 布莱克曼窗
% win0 = triang(n); %三角窗
% win0 = bartlett(n); %巴特利特窗
 win0 = chebwin(n,att);%切比雪夫窗
%% 乘积法
Win = win0*win0.';
Win_freq = fft2(Win);

%% 时域旋转法
w1_freq = fft(win0);
w2 = getWin2D(win0);
w2_freq = fft2(w2);

%% 频域旋转
w2_freq2 = fftshift(getWin2D(fftshift(w1_freq)));
w22 = ifft2(w2_freq2);
w22 = mapminmax(w22,0,1);
w22_freq = fft2(w22);
%% 画图
figure; 
subplot(2,1,1);plot(0:n-1,win0);
subplot(2,1,2);plot(0:n-1,abs(fftshift(w1_freq)));title("一维窗")
figure;  
subplot(2,1,1);mesh(0:n-1,0:n-1,abs(Win));
subplot(2,1,2);mesh(0:n-1,0:n-1,abs(fftshift(Win_freq)));title("乘积法")
figure;
subplot(2,1,1);mesh(0:n-1,0:n-1,abs(w2));
subplot(2,1,2);mesh(0:n-1,0:n-1,abs(fftshift(w2_freq))); title("时域旋转法")
figure; 
subplot(2,1,1);mesh(0:n-1,0:n-1,abs((w22)));
subplot(2,1,2);mesh(0:n-1,0:n-1,abs(fftshift(w22_freq)));title("频域旋转法")

% figure;mesh(0:n-1,0:n-1,abs(abs(Win)));
% figure;mesh(0:n-1,0:n-1,abs(abs(w2))); 
% figure;mesh(0:n-1,0:n-1,abs(abs(w22)));


%% PSR 
PSR_1 = get_psr(Win_freq);
disp(['直接乘积法所得2D窗的PSR为',num2str(PSR_1),'dB']);
% rho_1 = (2*norm(Win, 1)*a_max)/(sqrt(pi)*n^2*norm(Win,2)*sigma_n);
% rho_1 = mag2db(rho_1);%转化为dB % rho_1 = 20*log10(rho_1);
PSR_2 = get_psr(w2_freq);
disp(['时域旋转法所得2D窗的PSR为',  num2str(PSR_2),'dB']);
% rho_2 = 2*norm(w2, 1)*a_max/(sqrt(pi)*n^2*norm(w2,2)*sigma_n);
% rho_2 = 20*log10(rho_2);%转化为dB
PSR_3 = get_psr(w2_freq2);
disp(['频域旋转法所得2D窗的PSR为',  num2str(PSR_3),'dB']);
% rho_3 = 2*norm(w22, 1)*a_max/(sqrt(pi)*n^2*norm(w22,2)*sigma_n);
% rho_3 = 20*log10(rho_3);%转化为dB
