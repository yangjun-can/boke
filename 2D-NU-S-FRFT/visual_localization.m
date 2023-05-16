function [] = visual_localization(N0,N1,Omega_gt,Omega)    
% visual_localization

% Input:
%   N0, N1: signal length of the two dimensions
%   Omega_gt: groud truth locations of frequencies
%   Omega: estimated locations of frequencies

% Output: �Ա�ͼ��

    V = zeros(N0,N1);
    [x0,y0] = ind2sub(size(V),find(V>=0));%V�зǸ�Ԫ�ص�����
    % find���÷���
    %1.find(A) ����A�з���Ԫ�ص�λ�ã�������������
    %2.find(A,N)����ǰN������Ԫ�ص�λ�ã�������������
    %3.find(A**) ��������ĳ����A**��Ԫ��λ�ã�������������
    %4.find(A,1,'last')����A���һ������ֵ��λ�ã�������������
    %5.[a1,a2] = find(a)���ҳ�a�з���Ԫ�����ڵ��к��У����ֱ�洢��a1��a2�У�
    %6.[a1,a2,v] = find(a)���ҳ�a�з���Ԫ�����ڵ��к��У��ֱ�洢��a1��a2�У������������v�У�
    
    %ind2sub��������߾������������ת��Ϊ��Ӧ���±ꣻ
    %sub2ind�������෴�����±�ת��Ϊ����������
    
    figure,
%     scatter(x0-1,y0-1,5,'fill');%����DFT��� % scatter(x,y,sz,c��...) ������(x,y)��ɢ��ͼ��sz��ʾ��Ĵ�С��c��ʾ�����ɫ��'filled'��ʾʵ�ĵ㣬mkr��ʾɢ������ͣ�Բ�Σ�Ĭ�ϣ������Σ�+�ȣ�
%     hold on;
    scatter(Omega_gt(:,1),Omega_gt(:,2),100,'Or');%������ʵ��Ƶ��λ��
    if ~isempty(Omega)  %~��ʾ�߼���������ǡ���sempty(A) ;�ж�A�Ƿ�Ϊ�գ����Ϊ�գ����Ϊ1������Ϊ0.
        hold on;
%         scatter(Omega(:,1),Omega(:,2),50,'black','fill');%���ƹ��Ƶ�Ƶ��λ��
        scatter(Omega(:,1),Omega(:,2),50,'b','x');%���ƹ��Ƶ�Ƶ��λ��
    end
    hold off;
    xlabel('$m_0$','Interpreter','LaTex');%ʹ��latex�ַ����ַ����е���ѧ��ʽ��һ��'$$'�������������ò���InterpreterΪLatex��
    ylabel('$m_1$','Interpreter','LaTex');
%      legend('DFT grid','Ground truth','Estimated');%��עͼ��
%    legend('��ʵƵ��λ��','����Ƶ��λ��');%��עͼ��
%     set(gca,'FontSize',20);%x,y���ע���ֶ���ı��С
%     legend('True Location','Estimated Location');%��עͼ��
%     axis tight; % �Զ�����x���y��ķ�Χʹͼ����������ռ��������ʾ�ռ�
    axis([0,2050,0,2050])
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