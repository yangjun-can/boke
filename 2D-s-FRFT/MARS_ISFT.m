function [Omega,A] = MARS_ISFT(Sig,Win,N0,N1,T ,epsilon, gamma, n_d, n_s)
% Robutst MARS_SFT algorithm; for FPS_SFT, set Win as rect window, epsilon =1e-10
% gamma = 1e-10; n_d = n_s = 1

% Input:
%   Sig: singal matrix 
%   Win: window matrix
%   N0, N1: signal length of the two dimensions
%   epsilon: the threshold of detecting significant frequencies in a slice
%   gamma: the threshold for 1-sparsity detection
%   n_d, n_s: parameters of $n_d$-out-of-$n_s$ detection

% Output: 
%   Omega: frequency locations
%   A: complex amplitude of each frequency 
%   P: debug info

Omega = [];
A = [];
i = 0;  
while i<T  % 迭代T次
    i = i+1;
    [ind,ha] = ISFT_INNER_vote(Sig,Win,Omega,A,N0,N1 ,epsilon, gamma, n_d, n_s);
    % ind 内循环估计的位置
    % ha 内循环估计的幅度
    if ~isempty(ha)    % 这次内循环有被估计的大频率
        if isempty(A)  % 之前没有被估计的大频率
            Omega = [Omega;ind]; 
            A = [A; ha];
        else
            [C,ia,ib] = intersect(Omega,ind,'rows');            
           % intersect(A,B,___,'rows') 和 C = intersect(A,B,'rows',___) 将A和B的每一行都视为单个实体，并返回A和B的共有行，但是不包含重复项。
           % 必须指定 A 和 B，setOrder 是可选的。C = A(ia,:) 且 C = B(ib,:)。
           % C 重复估计的位置  ia重复位置在A中的坐标， ib重复位置在B中的坐标
            if ~isempty(C)  % 如果有重复定位，该位置的幅度叠加
                A(ia) = A(ia)+ha(ib);
                df = setdiff(1:size(ind,1), ib);
                %c = setdiff(A, B) 返回在A中有，而B中没有的值，结果向量将以升序排序返回。
                Omega = [Omega; ind(df,:)]; % 将剩余不重复的位置加入已估计集
                A = [A; ha(df)];
            else       % 如果没有重复定位，将所有位置加入已估计集
                Omega = [Omega;ind];
                A = [A; ha];
            end            
        end    
    end
end
end