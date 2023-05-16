function [Sig, A, u, v] = genSigOnGrid(N0,N1,K,a_min,sigma_n)
% Generate 2-D signals with on-grid uniformly distributed frequencies randomly 

% Input: 
%   N0, N1: signal length of the two dimensions
%   K: sparsity level
%   a_min: the minimum sigal amplitude
%   sigma_n: noise variance 

% Output: 
%   Sig: signal matrix (N0-by-N1)
%   A: amplitude of each frequency (K-length vector)
%   u,v: frequency locations

%% 
%�������Ƶ�ʵ�K�����λ��
    u = randi(N0,[K,1])-1;% �������K����0��N0-1��֮����ȷֲ����������У�
    v = randi(N1,[K,1])-1;% �������K����0��N1-1��֮����ȷֲ����������У�
      %rand()������������ȡֵ�ڣ�0~1��֮����ȷֲ���α�������
      %randn()�����������ɱ�׼��̬�ֲ���α�������
      %randi(imax,m,n)�����������ɾ��ȷֲ���m*n��α����������󣬷�ΧΪ��imin~imax���������䣩����iminȱʡ��Ĭ��Ϊ1��
    ind = [u,v];
    uInd = unique(ind,'rows');%Ƶ��λ�ã�ÿ��λ��ֻ����һ�Σ�
    %unique�÷���
    %1.C = unique(A) ������A����ͬ�����ݣ����ǲ������ظ��������C�Ѱ��մ�С��������
    %2.C = unique(A,___,'rows') �� C = unique(A,'rows',___) ��A�е�ÿһ����Ϊ����ʵ�壬��������˳�򷵻�A�е�Ψһ�С�����ָ��A����setOrder��occurrence�ǿ�ѡ�ġ�'rows' ѡ�֧��Ԫ�����顣
    %3.[C,ia,ic] = unique(A,___) ����ʹ�������κ��﷨������������ ia �� ic��iaΪ����C�е�Ԫ���ھ���A�е�λ�ã�icΪ����A�е�Ԫ���ھ���C�е�λ�á�
    while size(uInd,1)<K
        %��Ƶ��λ�ø�������һ�е�������<K,�����������λ�á�
        u = randi(N0,[K,1])-1;
        v = randi(N1,[K,1])-1;
        uInd = [uInd; [u,v]];
        uInd = unique(uInd,'rows');
    end
    rdI = randperm(length(uInd));%��Ƶ��λ���������˳��
    %y = randperm(n) :��1��n��Щ��������ҵõ���һ���������У�����������
    %y=randperm(n,m)��1-n�����ѡ��m������nһ�����ڵ���m
    %size ��ȡ���������������
    %length ��ȡ���鳤�ȣ��������������еĽϴ�ֵ���൱��max(size(a))
    %numel ����Ԫ������
    uInd = uInd(rdI(1:K),:);%ȡ��ǰK��λ��    
    oriInd = sortrows(uInd);%��λ�ð�����������
    u = oriInd(:,1);%����λ�ã�
    v = oriInd(:,2);%����λ�ã�
%%
 %����Ƶ���źż�ʱ���ź�
    A = a_min*exp(1i*2*pi*rand(K,1));%��������ź���Ƶ��������1i������piԲ����
    Sig=zeros(N0,N1);
    Sig(sub2ind([N0,N1],u+1,v+1)) = A; %ϡ���Ƶ���ź�  

    %ind2sub��������߾������������ת��Ϊ��Ӧ���±ꣻ
    %sub2ind�������෴�����±�ת��Ϊ����������ind = sub2ind(sz,row,col) ��Դ�СΪ sz �ľ��󷵻��� row �� col ָ���������±�Ķ�Ӧ�������� ind    
    %Sig = ifft2(Sig);%��Ӧ��ʱ���ź�   %iff2 ��ά���ٸ���Ҷ��任
    %Sig = fft2(Sig);%��Ӧ��Ƶ���ź�  
    p1=0.8;p2=0.8;ts1=0.005;ts2=0.005;
    alpha=p1*pi/2;     S1=sign(sin(alpha));
    beta=p2*pi/2;     S2=sign(sin(beta));
    us1=2*pi*S1*sin(alpha)/(ts1*N0);
    us2=2*pi*S2*sin(beta)/(ts2*N1);
    [Sig,~,~] = DFRFT_2D_fft2(Sig,-p1,-p2,N0,N1,us1,us2); %ʱ���ź�
    if  sigma_n==1
    noise = sigma_n*sqrt(2)/2*(randn(N0,N1)+1i*randn(N0,N1));  % ʱ��Ӹ���˹����
    % ����ʱ�������
%     Ps=sum(abs(Sig).^2);
%     Pn=sum(abs(noise).^2);
%     SNR_time=10*log10(Ps/Pn)
    SNR_t=snr(Sig,noise)
    Sig = Sig+noise;%����ʱ���ź�
    end
%% Visualization       
%     figure,mesh(0:N1-1,0:N0-1,abs(fft2(Sig)));
    %abs����:��ֵ�ľ���ֵ�͸����ķ�ֵ
    %mesh��仭����ͼƬ��ʵ���Ͼ��Ǹ���һ������(x,y)����������z(x,y)��ֵ��
%     xlabel('$v$','Interpreter','LaTex'), ylabel('$u$','Interpreter','LaTex'), zlabel('Amplitude');
%     axis tight; % �Զ�����x���y��ķ�Χʹͼ����������ռ��������ʾ�ռ�
%     set(gca,'FontSize',20);%x,y���ע���ֶ���ı��С
end