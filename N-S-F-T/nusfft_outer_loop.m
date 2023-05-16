%% ��ѭ��
function [I,out_I,out]=nusfft_outer_loop(origx, n, filter_time,filter_sizet,  ...
    filter_freq,filter_est_time,filter_est_sizet, filter_est_freq, B2,num,B,...
    loop_threshold, location_loops, loops,pe_major,inv_pe_major,sa_major,f)
% function [I,out_I,out]=nusfft_outer_loop(origx, n, filter_time,filter_sizet,  ...
%     filter_freq,filter_est_time,filter_est_sizet, filter_est_freq, B2,num,B,...
%     loop_threshold, location_loops, loops,sa_major,permutea,permuteb,permuteMateix,f)
%% ����
% origx ԭʼ�ź�
%  n    �źų���
% filter_time ƽ̹���˲�����ʱ��
% filter_sizet ƽ̹���˲�������
% filter_freq ƽ̹���˲�����Ƶ��
% B  loc��B
% B2  est��B
% num ��Ͱ�ĸ���
% W_Comb  ����˲�������B
% Comb_loops  SFFT2.0�л���˲��ĸ���
% location_loops  locѭ������
% loop_threshold ��ֵ�ı���Ǵ�����ֵ
% loops   ��ѭ���ܴ���
%output:
% I     ��ֵ��λ�ã�����
% out_I ��Ӧ�Ĺ�ֵ
% out   ��ֵ����λ�õĺ���

%%
tic;
x_samp=zeros(loops,max(B,B2),'single');          % ����ÿ��ѭ����Ͱֵ
permutea=zeros(1,loops);
permuteb=zeros(1,loops);
score =zeros(1,n,'single');
hits_ini=zeros(1,n,'single');
hits_found_ini=0;
%% ��ѭ��
for i=1:loops
    %     a=permutea(i);
    %     b=permuteb(i);
    %     permute=permuteMateix(:,:,i);
    %% �趨����
    a = 0;                                  % ������Ų���sigma��ģ��
    b = 0;							        % ����ʱ���±�ƽ�Ʋ��� ���൱�ڿα��ϵ�-a*sigma��
    while(undivide(a, n)~= 1 || mod(a,2)==0 || mod(a,5)==0)
        a = mod(fix(rand()*n) ,n);			% a��n������������ѭ��
    end
    ai = mod_inverse(a, n);     	        % ģ�� a*ai mod n = 1
    permutea(i)=ai;						    % ����ÿ��ѭ����������Ų���ai
    permuteb(i) = b;					    % ����ÿ��ѭ����������Ų��� b
    
    perform_location =(location_loops>=i);
    if perform_location                     % ���ö�λѭ��������
        cur_filter_time=filter_time;
        cur_filter_sizet=filter_sizet;
        %         cur_filter_freq=filter_freq;
        cur_B=B;
    else                                     % ���ù�ֵѭ��������
        cur_filter_time=filter_est_time;
        cur_filter_sizet=filter_est_sizet;
        %         cur_filter_freq=filter_est_freq;
        cur_B=B2;
    end
    flat_filter_time=zeros(n,1,'single');
    flat_filter_time(1:cur_filter_sizet)=cur_filter_time;
    %% ��Ͱ,���ҵ���Ͱ��
    [x_samp(i,1:cur_B),J]=nusfft_inner_loop_locate(origx,n,flat_filter_time, ...
        num, cur_B,ai,b,pe_major,inv_pe_major,sa_major,f);
    %     [x_samp(i,1:cur_B),J]=nusfft_inner_loop_locate(origx,n,flat_filter_time, ...
    %         num, cur_B,sa_major,permute,f);
    
    % x_samp�����������Ƶ��ֵ��Ͱֵ��
    % J����Ͱ���±꼯
    %%  �ۼƴ�ֵԭʼλ�ü���
    if  perform_location
        [hits,hits_found,score_out ]=inner_loop_filter_regular(J, n, num, cur_B, a, loop_threshold, score,hits_ini, hits_found_ini);
        hits_ini=hits;                  % �ҵ������еĴ�ֵλ�ü�
        hits_found_ini=hits_found;      % �ҵ��Ĵ�ֵλ������
        score=score_out;                % ��Ƿ�����
    end
end
% word=sprintf('Number of candidates: %d\n', hits_found);
% disp(word);
%% ��ֵ
[I,out_I,out] = estimate_values(hits, hits_found, x_samp,  loops, n, permutea, ...
    B, B2, filter_freq, filter_est_freq,  location_loops);
%% ��ʱ
% disp('NUSFFT������ʱ��Ϊ��');
toc
end