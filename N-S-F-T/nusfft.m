function [out_Large,BB_loc,BB_est,b_loc,b_est]=nusfft(x,k,f)
%% ������ʼ��
% profile on;
x=single(x);
f=single(f);
n=length(x);           % �źų���
%k                     % �ź�ϡ���
%repetitions = 1;	   % ��ѭ�����д���
Bcst_loc=2;            % loc��Ͱ����
Bcst_est=2;            % est��Ͱ����
loc_loops =16;		   % locѭ������
est_loops =6;          % estѭ������
tolerance_loc = 1e-4;     % loc������й¶������ test��n=64,ȡ1e-4�� test��n>64,ȡ1e-6��
tolerance_est = 1e-4;     % est������й¶������
threshold_loops =fix(loc_loops/2);	   % ��ֵ�ı���Ǵ�����ֵ

loops=loc_loops+est_loops;% �ܵ���ѭ������
BB_loc = fix(Bcst_loc*sqrt(n*k/(log2(n))));
B_loc = floor_to_pow2(BB_loc);                 % loc����B
BB_est = fix(Bcst_est*sqrt(n*k/(log2(n))));
B_est = floor_to_pow2(BB_est);                 % est����B
B_thresh = 2*k;                                % ��Ͱ�ĸ���
lobefrac_loc = 0.5 / BB_loc;				   % locͨ������[-lobefrac_loc,lobefrac_loc]
lobefrac_est = 0.5 / BB_est;                   % estͨ������[-lobefrac_est,lobefrac_est]
b_loc = fix(1.4*1.1*(n/BB_loc));
b_est = fix(1.2*1.1*(n/BB_est));
% Ԥ�������ž���
% permuteMateix=zeros(n,n,loops);
% sequence1=gpuArray(single(exp(-2*1i*pi*(0:n-1).*f(1:n)/n))); 
%  sequence1=gpuArray(exp(-2i*pi*(0:n-1).*f(1:n)/n));
 sequence1=exp(-2i*pi*(0:n-1).*f(1:n)/n); 
pe_major=n*ifft(sequence1);
inv_pe_major=inv(pe_major);
% Ԥ�����²�������
% sequence2=gpuArray(single(exp(-2i*pi*(0:B_loc-1).*(f((1:B_loc)*n/B_loc)-(1:B_loc)'*(n/B_loc-1))/B_loc)));  % Ԥ�����²�������
%  sequence2=gpuArray(exp(-2i*pi*(0:B_loc-1).*(f((1:B_loc)*n/B_loc)-(1:B_loc)'*(n/B_loc-1))/B_loc));
 sequence2=exp(-2i*pi*(0:B_loc-1).*(f((1:B_loc)*n/B_loc)-(1:B_loc)'*(n/B_loc-1))/B_loc);  % Ԥ�����²�������
sa_major1=B_loc*ifft(sequence2); 
sa_major=repmat(sa_major1,[1,n/B_loc]);
sa_major=sa_major/pe_major;

% word=sprintf('\n\nRUNNING EXPERIMENT: n=%d, k=%d.\n', n, k);
% disp(word);                                  % ���n,k
% word=sprintf('\n\nSimulation:\n');
% disp(word);
% word=sprintf('*******************************************************\n');
% disp(word);

%% NUSFFT����Ƶ��
%% ��Ͱ��������Ͱ������Ǵ���������ѭ������
if B_thresh > B_loc ||threshold_loops >loc_loops
    return;
end
%% ����ƽ̹���˲���
[filtert,w_loc ]=make_dolphchebyshev_t(lobefrac_loc, tolerance_loc);                 % ���� filtert �б�ѩ���˲�����w_loc������
[filter_time, filter_sizet, filter_freq]=make_multiple_t(filtert, w_loc, n, b_loc);    % ����ƽ̹���˲�����filter_time���˲���ʱ�� filter_sizet���˲������� filter_freq���˲���Ƶ��
filter_time=single(filter_time);
filter_freq=single(filter_freq);
% aaaaa=fft(filtert);
% figure;plot(abs(aaaaa));title('dolphchebyshev filter time domain - window function');
% figure;
% subplot(2,2,1);plot(abs(filter_time));title('(a)');ylabel('amplitude');
% subplot(2,2,3);plot(abs(filter_freq));title('(c)');ylabel('amplitude');
% subplot(2,2,2);semilogy(abs(filter_time));title('(b)');ylabel('amplitude');
% subplot(2,2,4);semilogy(abs(filter_freq));title('(d)');ylabel('amplitude');
[filtert_est, w_est]=make_dolphchebyshev_t(lobefrac_est, tolerance_est);
[filter_est_time ,filter_est_sizet ,filter_est_freq]=make_multiple_t(filtert_est, w_est, n, b_est);
filter_est_time=single(filter_est_time);
filter_est_freq=single(filter_est_freq);
% word=sprintf(' Window size: Location Filter : %d; Estimation Filter : %d;\n', w_loc, w_est);
% disp(word);

%% ƽ̹���˲�����й¶
% filter_noise = 0;
% filter_noise_est = 0;
% for i=0:9
%     filter_noise = max(filter_noise,max(abs(filter_freq(n/2+i+1)),abs(filter_freq(n/2-i+1))));
%     filter_noise_est = max(filter_noise_est,max(abs(filter_est_freq(n/2+i+1)),abs(filter_est_freq(n/2-i+1))));
% end
% word=sprintf('Noise in filter: Location Filter : %d; Estimation Filter %d\n', filter_noise, filter_noise_est);
% disp(word);
% word=sprintf('****************************************************************\n\n');
% disp(word);
% word=sprintf('sFFT Results\n**************************************************************');
% disp(word);
%% ��ѭ��
% for ii=1:repetitions
    % ��ѭ����ÿ�εõ�heavyֵ��λ�úʹ�С
        [I,out_I,out]= nusfft_outer_loop(x, n,filter_time,filter_sizet, filter_freq, filter_est_time,...
        filter_est_sizet,filter_est_freq, B_est, B_thresh, B_loc,threshold_loops, ...
        loc_loops,loops,pe_major,inv_pe_major,sa_major,f);
    
    % I     ��ֵ��λ�ã�����
    % out_I ��Ӧ�Ĺ�ֵ
    % out   ��ֵ����λ�õĺ���
% end
%%
num_candidates= length(out_I);      % �����ƵĴ�ֵ����
out_Large = zeros(1,n);
candidates=zeros(2,num_candidates);
% counter=0;
%% �ù��ƽ�������źŵ�Ƶ��
% for i=1:num_candidates                   
%     counter=counter+1;
%     candidates(1,counter)=abs(out_I(i));  % ��ֵ�ķ�������
%     candidates(2,counter)=I(i);           % ��ֵ��λ������
% end 
% ���㷨���ƵĽ�������ڶ�ά���󣬵ڶ�����λ�ã���һ������Ӧ�ķ�ֵ
candidates(1,:)=abs(out_I);  % ��ֵ�ķ�������
candidates(2,:)=I;           % ��ֵ��λ������
temp=sortrows((candidates).').';          % ��ֵ����ֵ��������
% for i=1:k                               % ȡK�����ķ�ֵ������ź�Ƶ��out_Large
%     key=temp(2,num_candidates - k+ i);
key=temp(2,num_candidates - k+ (1:k));
out_Large(key+1) =out(key+1)*n;
% end
% profile viewer;
end