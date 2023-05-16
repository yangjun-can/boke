%% 拟合SNR图像
% x=-8.5:0.5:10;
% x=[x,11:1:25];
% y=[0.060226 0.10682 0.05379	0.087187 0.059532 0.073004 0.11069 0.031701 ...
%     0.070365 0.056156 0.036948 0.043834	0.034016 0.029118 0.027979 0.010564...
%     0.055375 0.044439 0.028492 0.059606 0.040557 0.022176 0.0091474	0.02219....
%     0.035766 0.010604 0.01021 0.026715 0.019633	0.0068886 0.013154 0.022161...
%     0.025988 0.030341 0.012326 0.017725	0.017766 0.019267 0.024582	0.019664...
%     0.021278	0.019998	0.020552	0.017067	0.02227	0.014579	...
%     0.025939	0.01948	0.018302	0.019781	0.023096	0.02014	0.022706];
% P = polyfit(x,y,3);
% yi = polyval(P,x);   %多项式求值
% plot(x,y);         %观测数据点
% hold on;
% plot(x,yi);
% title('Robustness vs. SNR');
% xlabel('SNR（dB）');ylabel('Average L1 error per entry');

 %% SNR和k与精度的关系
% x=[5 10	15	20	25	30	35];  % SNR
% y=[0.00038417	0.00030646	0.00020249	0.00011819	6.76e-05 3.65e-05 3.38e-05
% 0.00099663	0.00059893	0.0003291	0.00021323	0.00015141	0.00010021	9.53e-05
% 0.0055791	0.0014784	0.00076688	0.00045154	0.00034338	0.00026118	0.00017014
% 0.0066785	0.0017349	0.0010533	0.00067483	0.00046051	0.0004171	0.00033242];% 不同k时对应x的l2误差
% figure;
% plot(x,y(1,:),'.-r');
% hold on;plot(x,y(2,:),'.--b'); 
% hold on;plot(x,y(3,:),'.-g'); 
% hold on;plot(x,y(4,:),'.-.m'); 
% xlabel('SNR（dB）');ylabel('Average L2 error per entry');
% legend('k=5','k=15','k=25','k=35');

%% 各个算法的鲁棒性
% N=[2^7,2^8,2^9,2^10,2^11,2^12,2^13];
% error35=[7.36E-05	6.87E-05	7.55E-05	6.98E-05	7.15E-05	7.39E-05	7.34E-05
% 7.16E-05	6.70E-05	7.05E-05	6.78E-05	6.88E-05	6.92E-05	6.94E-05
% 0.0078162	0.0078166	0.0078089	0.0078123	0.0078128	0.0078118	0.0078112
% 9.11E-06	6.39E-06	1.09E-05	9.09E-06	4.49E-06	3.33E-06	6.76E-06]; % SNR=35时的算法误差
% figure;
% plot(N,error35(1,:),'.-m');
% hold on;plot(N,error35(2,:),'.--b'); 
% hold on;plot(N,error35(3,:),'.-g'); 
% hold on;plot(N,error35(4,:),'.-r'); 
% xlabel('Signal size');ylabel('Average L2 error per entry');
% legend('Min-Max interpolation','Low rank approximation','Fast gaussian gridding','NUSFT');
% 
% error10=[1.35E-03	1.34E-03	0.0014637	0.001443	0.0014463	0.0014318	0.0014303
% 1.24E-03	1.23E-03	0.0012736	0.001243	0.0012463	0.0012418	0.0012303
% 0.0080445	0.0079683	0.0077663	0.0079248	0.007932	0.007921	0.007821
% 3.68E-05	5.14E-05	3.77E-05	4.14E-05	5.58E-05	3.04E-05	4.94E-05]; % SNR=10时的算法误差
% figure;
% plot(N,error10(1,:),'.-m');
% hold on;plot(N,error10(2,:),'.--b'); 
% hold on;plot(N,error10(3,:),'.-g'); 
% hold on;plot(N,error10(4,:),'.-r'); 
% xlabel('Signal size');ylabel('Average L2 error per entry');
% legend('Min-Max interpolation','Low rank approximation','Fast gaussian gridding','NUSFT');
%% 各算法的效率
%  N=[2^7,2^8,2^9,2^10,2^11,2^12,2^13];
% times=[0.010338	0.012818	0.040913	0.19922	1.0388	4.3139	21.501
% 0.010928	0.01102	0.011192	0.012319	0.014273	0.022264	0.03771
% 0.007049	0.007585	0.007758	0.007797	0.008721	0.017626	0.029919];
% figure;
% plot(N,times(2,:),'.-b');
% hold on;plot(N,times(3,:),'.-r'); 
% % hold on;plot(N,times(1,:),'.-m'); 
% xlabel('Signal size');ylabel('Execution time/s');
% % legend('Direct','Min-Max interpolation','NUSFT');
% legend('Min-Max interpolation','NUSFT');
%% k,N与时间的关系
% k=1:15;
% T=[0.007049	0.007057	0.007066	0.007088	0.007115	0.007166	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
% 0.007585	0.007601	0.007611	0.007669	0.007682	0.007708	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
% 0.007758	0.007794	0.007805	0.007848	0.007813	0.007943	0.008003	0.008146	0.008171	0.008205	0.008283	0.008324	0.008355 NaN	NaN
% 0.007797	0.007822	0.00791	0.007966	0.008011	0.008053	0.008122	0.008204	0.008221	0.00825	0.008316	0.008367	0.008371	0.008426	0.008453
% 0.008721	0.00875	0.00909	0.009102	0.009292	0.009303	0.009374	0.009395	0.009454	0.0095	0.009541	0.009567	0.009641	0.009761	0.009799
% 0.017626	0.018117	0.01815	0.018498	0.018604	0.018764	0.019057	0.019133	0.01921	0.019366	0.019442	0.019513	0.019522	0.019538	0.019684];
% figure;
% plot(k,T(1,:),'.-b');
% hold on;plot(k,T(2,:),'.-g'); 
% hold on;plot(k,T(3,:),'.-r'); 
% hold on;plot(k,T(4,:),'.-k'); 
% hold on;plot(k,T(5,:),'.-c'); 
% % hold on;plot(k,T(6,:),'.-m'); 
% xlabel('Signal sparity');ylabel('Execution time/s');
% legend('$N=2^7$','$N=2^8$','$N=2^9$','$N=2^{10}$','$N=2^{11}$','Interpreter','LaTex');
%% 低SNR时的恢复概率
% SNR=[10 5 0 -2.5 -5 -7.5 -9 -10];
% Prob=[94.59	100	96.15	85.72	48.28	47.62	43.40	32.26
% 100	96.77	100	96.77	96.77	88.24	81.08	78.95
% 100	96.77	100	100	96.77	96.77	68.18	53.13];
% figure;
% plot(SNR,Prob(1,:),'.-b');
% hold on;plot(SNR,Prob(2,:),'.--k'); 
% hold on;plot(SNR,Prob(3,:),'.-r'); 
% % hold on;plot(SNR,Prob(4,:),'.-k'); 
% % hold on;plot(SNR,T(5,:),'.-c'); 
% % hold on;plot(SNR,T(6,:),'.-m'); 
% xlabel('SNR(dB)');ylabel('Probability of consistent detection(%)');
% legend('$k/N=1/64$','$k/N=1/256$','$k/N=1/512$','Interpreter','LaTex');
%% 频率值差33dB时的重排图像 
%  % 提取图像数据
lh=findall(gca,'type','line'); % 从当前图(gca)中取出曲线的handle
xc=get(lh,'xdata'); % 取出x轴数据，注意，这个x和y是以cell的数据结构保存的
yc=get(lh,'ydata');
% figure;plot(xc{2,1},yc{2,1});%xlabel('Non-uniform estimation point in the frequency');
% ylabel('amplitude');

% %加横线
% set(gca,'YTick',sort([get(gca,'YTick') 35400])); %纵坐标上标记35500
% line([0,600],[35400,35400],'linestyle','-','color','r');

 add=[4.69e-4, 6.8e-4, 4.96e-4, 4.01e-4, 3.2e-4, 2.04e-4, 0.75e-4, 0.04e-4];
 figure;
plot(xc{4,1},yc{4,1}+add,'--b');hold on;
plot(xc{3,1},yc{3,1},'-m');hold on;
plot(xc{2,1},yc{2,1},'-k');hold on;
plot(xc{1,1},yc{1,1},':r','LineWidth',2);hold off
xlabel('SNR(dB)');ylabel('Average L2 error per entry');
legend('Low rank approximation','Min-Max interpolation','Fast gaussian gridding','NUSFT');
