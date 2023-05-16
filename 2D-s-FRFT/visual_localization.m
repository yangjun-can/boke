function [] = visual_localization(N0,N1,Omega_gt,Omega)    
% visual_localization

% Input:
%   N0, N1: signal length of the two dimensions
%   Omega_gt: groud truth locations of frequencies
%   Omega: estimated locations of frequencies

% Output: 对比图象

    V = zeros(N0,N1);
    [x0,y0] = ind2sub(size(V),find(V>=0));%V中非负元素的坐标
    % find的用法：
    %1.find(A) 返回A中非零元素的位置（线性索引）；
    %2.find(A,N)返回前N个非零元素的位置（线性索引）；
    %3.find(A**) 返回满足某条件A**的元素位置（线性索引）；
    %4.find(A,1,'last')返回A最后一个非零值的位置（线性索引）；
    %5.[a1,a2] = find(a)，找出a中非零元素所在的行和列，并分别存储在a1和a2中；
    %6.[a1,a2,v] = find(a)，找出a中非零元素所在的行和列，分别存储在a1和a2中，并将结果放在v中；
    
    %ind2sub把数组或者矩阵的线性索引转化为相应的下标；
    %sub2ind则正好相反，将下标转化为线性索引。
    
    figure,
%     scatter(x0-1,y0-1,5,'fill');%绘制DFT格点
    % scatter(x,y,sz,c，...) 绘制下(x,y)的散点图，sz表示点的大小，c表示点的颜色，'filled'表示实心点，mkr表示散点的类型（圆形（默认），菱形，+等）
%     hold on;
    scatter(Omega_gt(:,1),Omega_gt(:,2),100,'Or');%绘制真实的频率位置
    if ~isempty(Omega)  %~表示逻辑运算符“非”；sempty(A) ;判断A是否为空，如果为空，结果为1，否则为0.
        hold on;
%         scatter(Omega(:,1),Omega(:,2),50,'black','fill');%绘制估计的频率位置
        scatter(Omega(:,1),Omega(:,2),50,'b','x');%绘制估计的频率位置
    end
    hold off;
%     xlabel('$m_0$','Interpreter','LaTex');%使用latex字符：字符串中的数学公式用一对'$$'包括起来，设置参数Interpreter为Latex。
%     ylabel('$m_1$','Interpreter','LaTex');
%     legend('DFT grid','Ground truth','Estimated');%标注图标
%     set(gca,'FontSize',20);%x,y轴标注文字都会改变大小
    legend('True Location','Estimated Location');%标注图标
    axis tight; % 自动设置x轴和y轴的范围使图形区域正好占满整个显示空间
    
    %% Add a small detailed graph
%     axes('Position',[.18 .18 .3 .3]);
%     x = 70:90;
%     y = 175:195;
%     [X,Y] = meshgrid(x,y);
%     
%     box on
%     scatter(reshape(X,[1,numel(X)]),reshape(Y,[1,numel(Y)]),1,'fill'),grid;
%     hold on;scatter(Ori(:,1),Ori(:,2),100,'Or');
%     hold on;scatter(I(:,1),I(:,2),50,'black','fill');
%     axis([x(1), x(end), y(1), y(end)]);
end