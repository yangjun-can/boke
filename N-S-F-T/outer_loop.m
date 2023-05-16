%% ϡ�踵��Ҷ�任����ѭ��
function [I,out_I,out]=outer_loop(origx, n, filter_time,filter_sizet,  ...
    filter_freq,filter_est_time,filter_est_sizet, filter_est_freq, B2,num,B,...
    loop_threshold, location_loops, loops)
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
% WITH_COMB ����SFTT2.0
% ALGORITHM1  ��ʹ�û���˲�����loc_loops=0;
% LARGE_FREQ
%output:
% I     ��ֵ��λ�ã�����
% out_I ��Ӧ�Ĺ�ֵ
% out   ��ֵ����λ�õĺ���

%%
permutea =zeros(1,loops);               % ����ÿ��ѭ����������Ų��� ai
permuteb =zeros(1,loops);               % ����ÿ��ѭ����������Ų��� b
x_samp=zeros(loops,max(B,B2));          % ����ÿ��ѭ����Ͱֵ
score =zeros(1,n);
hits_ini=zeros(1,n);
hits_found_ini=0;
%% ��ѭ��
for i=1:loops
    %% �趨����
    a = 0;                                  % ������Ų���sigma��ģ��
    b = 0;							        % ����ʱ���±�ƽ�Ʋ��� ���൱�ڿα��ϵ�-a*sigma��
    while(undivide(a, n)~= 1 || mod(a,2)==0 || mod(a,5)==0)
        a = mod(fix(rand()*n) ,n);			% a��n������������ѭ��
    end
    ai = mod_inverse(a, n);			        % ģ�� a*ai mod n = 1
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
   
    %% ��Ͱ,���ҵ���Ͱ��
    [x_samp(i,1:cur_B),J ]=inner_loop_locate(origx,n,cur_filter_time,cur_filter_sizet, ...
        num, cur_B,ai, b);              
    % x_samp�����������Ƶ��ֵ��Ͱֵ��
    % J����Ͱ���±꼯
    % figure;plot(abs(fft(x_samp(i,1:cur_B))))��
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

end