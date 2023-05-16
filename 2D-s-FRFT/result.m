% result of corr 2D SFRFT

%% 定位误差的阈值 gama  与SNR、定位P的关系
% % N=256
% x=[0.95 0.98 0.99 0.999 0.9999];  % P
% y=[1	2	2	2	3
% 4	5	5	7	8
% 7	9	10	12	15
% 13	16	17	22	27
% ];% 不同SNR时对应的误差阈值
% figure; plot(x,y(1,:),'.-r');
% hold on;plot(x,y(2,:),'.--b'); 
% hold on;plot(x,y(3,:),'.-g'); 
% hold on;plot(x,y(4,:),'.-.m'); 
% xlabel('Minimum probability of the correct set of positions P_1^\gamma');
% ylabel('Position error threshold');
% legend('SNR = 11.8351dB','SNR = 1.8384dB','SNR = -3.1534dB','SNR = -8.1777dB');

%% 定位误差的阈值 gama  与 N1 N2
% x1=[64	128	256	512	1024	2048];  % N1
% x2=[64	128	256	512	1024	2048];  % N2
% y=[1	2	2	1	1	1
% 1	2	2	1	1	1
% 1	1	3	3	2	2
% 0	1	2	4	5	4
% 0	1	1	3	5	7
% 0	0	1	2	4	8
% ];% 不同N2时对应的误差阈值
% figure;
% mesh(x1,x2,y);
% %surf(x1,x2,y);
% % plot(x,y(1,:),'.-r');
% % hold on;plot(x,y(2,:),'.--b'); 
% % hold on;plot(x,y(3,:),'.-g'); 
% % hold on;plot(x,y(4,:),'.-.m'); 
% % hold on;plot(x,y(5,:),'.-.k'); 
% % hold on;plot(x,y(6,:),'.-.c'); 
% ylabel('N_2');
% xlabel('N_1');
% zlabel('Position error threshold');
% colormap("hsv")
% axis([64 2048 64 2048]); 
% % legend('N_2 = 64','N_2 = 128','N_2 = 256','N_2 = 512','N_2 = 1024','N_2 = 2048');
%% 频率的SNR
%  k=[10	20	30	40	50	60	80	100];
% SNR_freq=[74	71	69	68	67	66	65	64
% 64	61	59	58	57	56	55	54
% 53	50	48	47	46	45	44	43];
% figure; plot(k,SNR_freq(1,:),'.-.r');
% hold on;plot(k,SNR_freq(2,:),'.--b'); 
% hold on;plot(k,SNR_freq(3,:),'.-.k'); 
% xlabel('Sparsity of the signals');
% ylabel('SNR of one significant frequency');
% legend('SNR = 17.7707dB','SNR = 7.7738','SNR = -3.229dB');grid


%%  平均各点的内循环次数
% k=[10	20	30	40	50	60	80	100];
% s=[0	5.2/20	9.4167/30	17.8571/40	18.875/50	33.3/60	32.2857/80	95/100
% 9.1/10	21.6/20	34.3/30	47/40	63.4167/50	78/60	117.4/80	163.6667/100
% 12.9/10	25.2941/20	46.2857/30	64.9167/40	90/50	119.9/60	197.1538/80 279/100];
% figure; plot(k,s(1,:),'.-.r');
% hold on;plot(k,s(2,:),'.--b'); 
% hold on;plot(k,s(3,:),'.-.k'); 
% xlabel('Sparsity of the signals');
% ylabel('Average number of inner loops in each position');
% legend('SNR = 17.7707dB','SNR = 7.7738','SNR = -3.229dB');

%% 算法迭代次数
% k=[10	20	30	40	50	60	80	100];
% T=[1	2	2	2	2	2	2	2
% 1	2	2	2	2	2	2	2
% 2	3	3	3	4	4	4	5];
% figure; plot(k,T(1,:),'.-.r');
% hold on;plot(k,T(2,:),'.--b'); 
% hold on;plot(k,T(3,:),'.-.k'); 
% xlabel('Sparsity of the signals');
% ylabel('Number of iterations of 2D SFRFT algorithm');
% legend('SNR = 17.7707dB','SNR = 7.7738','SNR = -3.229dB');

%% 恢复概率
% k=[10	20	30	40	50	60	80	100];
% prob=[100	99	100	100	98.5	98	97	98
% 100	100	99	100	97	96	92	89
% 100	97	99	94	97	88	54	nan ];
% figure; plot(k,prob(1,:),'.-.r');
% hold on;plot(k,prob(2,:),'.--b'); 
% hold on;plot(k,prob(3,:),'.-.k'); 
% xlabel('Sparsity of the signals');
% ylabel('Probability of correct estimation of 2D SFRFT');
% legend('SNR = 17.7707dB','SNR = 7.7738','SNR = -3.229dB');
%% 时间、误差对比
N1_all=[64,128,256,512,1024,2048,4096];
time_direct=[0.242531500000000	1.67538730000000	18.4859515000000	102.233996300000	814.551410400000	7055.80006890000	62146.0679990000];
time_Pei=[0.00287020000000000	0.00467100000000000	0.0252756000000000	0.0431075000000000	0.143060500000000	0.609236400000000	2.19535590000000];
time_sFRFT=[0.00134190000000000	0.00176820000000000	0.00542520000000000	0.0225301000000000	0.0410151000000000	0.158328400000000	0.615366200000000];

error_Pei=[4.65660583482518e-09	4.94764579495902e-09	4.18367186954048e-09	4.80212549687131e-09	4.07453344148336e-09	4.11102888445432e-09	4.65706718303474e-09];
error_sFRFT=[3.23535001419941e-13	2.22946339190397e-13	1.54059335044791e-12	1.46921729644062e-13	2.06217837045123e-12	9.67915003171528e-12	2.21040360538275e-11];

figure; plot(N1_all,time_direct,'.-.k');hold on; 
plot(N1_all,time_Pei,'.--b');
hold on; plot(N1_all,time_sFRFT,':r','LineWidth',2);
xlabel('Size of the Rows and columns');
ylabel('Execution time/s');
legend('Decomposition and direct method','Decomposition and Pei method','2D SFRFT method');
axis([64,4096,0,62400])

figure;
plot(N1_all,time_Pei,'.--b');
hold on; plot(N1_all,time_sFRFT,':r','LineWidth',2);
axis([64,4096,0,2.5])
xlabel('Size of the Rows and columns');
ylabel('Execution time/s');
legend('Decomposition and Pei method','2D SFRFT method');

figure;plot(N1_all,error_Pei,'.--b');
hold on; plot(N1_all,error_sFRFT,':r','LineWidth',2);
axis([64,4096,0,5e-9])
xlabel('Size of the Rows and columns');
ylabel('$L_2$ Errors',Interpreter='latex' );
legend('Decomposition and Pei method','2D SFRFT method');

%% 迭代次数 
% k=[100,200,300,400,500,600];
% T=[3 4 5 7 10 16    
%     2 3 4 4 5 6
%     2 3 3 4 4 5];
% figure; plot(k,T(1,:),':b','LineWidth',1.5);hold on; 
% plot(k,T(2,:),':m','LineWidth',1.5);
% hold on; plot(k,T(3,:),':g','LineWidth',1.5);
% xlabel('Signal sparsity');
% ylabel('Number of iterations');
% legend('${N_1} = 256,{N_2} = 256$','${N_1} = 512,{N_2} = 128$','${N_1} = 1024,{N_2} = 64$',Interpreter='latex');

