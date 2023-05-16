clc
clear
close
%% 切比雪夫窗的PSR 与SNR、稀疏度 K 的关系
% N=256
% x=[5 10 20 30 40 50 60];  % 切比雪夫窗PSR
% y=[68	87.5	93.5	95.5	92.5	89	83.5
% 91.5	97	99.5  100	97	95.5	93.5
% 86	94.5	98.5  97	95.5	94	92];% 不同PSR时对应的第一次迭代的一致检测概率
% figure; plot(x,y(1,:),'r-*','Linewidth', 1.5);
% hold on;plot(x,y(2,:),'b-*','Linewidth', 1.5); 
% hold on;plot(x,y(3,:),'g-*','Linewidth', 1.5); 
% % hold on;plot(x,y(4,:),'.-.m'); 
% xlabel('窗函数的PSR');
% ylabel('一致检测概率');
% legend('SNR=30dB,k=10','SNR=30dB,k=5','SNR=20dB,k=5');
% axis tight;
%% 第一次迭代的一致检测概率 vs 非均匀频率位置  SNR 最大幅值
x=[0.1  0.2 0.3 0.4 0.5 ];  % 非均匀最大距离
y=[95	93	92.5	93	94
97	93	93.5	94	91.5
85.5	87.5	85	83.5	84.5];% 不同距离时对应的第一次迭代的一致检测概率
figure; plot(x,y(1,:),'r-*','Linewidth', 1.5);
hold on;plot(x,y(2,:),'b-*','Linewidth', 1.5); 
hold on;plot(x,y(3,:),'g-*','Linewidth', 1.5); 
% hold on;plot(x,y(4,:),'.-.m'); 
xlabel('频率扰动上限');
ylabel('一致检测概率');
legend('$SNR = 30dB,a_{max} = 10$','$SNR = 30dB,a_{max} = 20$','$SNR = 20dB,a_{max} = 20$','interpreter','latex');
% axis tight;
axis([0,0.6,50,100])