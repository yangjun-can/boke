function [out_Large]=sft(x,k)
% ERROR1 ����ֵƽ�������
% ERROR2 ��ֵ������
%% ������ʼ��
n=length(x);           % �źų���
%k =2;                 % �ź�ϡ���
repetitions = 1;	   % ��ѭ�����д���
Bcst_loc=8;            % loc��Ͱ����
Bcst_est=8;            % est��Ͱ����
Comb_cst=16;           % SFFT2.0����˲����ķ�Ͱ����
loc_loops =3;		   % locѭ������
est_loops =8;          % estѭ������
threshold_loops =2;	   % ��ֵ�ı���Ǵ�����ֵ
Comb_loops = 1;        % SFFT2.0�л���˲��ĸ���
simulate = 0;		   % ���������������򣬷���SFFT������ʱ��ʹ���
snr=100000;            % �����
std_noise = 0;         % �����Ŵ�ֵ��std_noise = sqrt(k/(2*snr));
FFTW_OPT = false;      % ʹ��FFTW��
tolerance_loc = 1e-6;  % loc������й¶������
tolerance_est = 1e-6;  % est������й¶������

BB_loc = fix(Bcst_loc*sqrt(n*k/(log2(n))));
B_loc = floor_to_pow2(BB_loc);                 % loc����B
BB_est = fix(Bcst_est*sqrt(n*k/(log2(n))));
B_est = floor_to_pow2(BB_est);                 % est����B

B_thresh = 2*k;                                % ��Ͱ�ĸ���

lobefrac_loc = 0.5 / BB_loc;				   % locͨ������[-lobefrac_loc,lobefrac_loc]
lobefrac_est = 0.5 / BB_est;                   % estͨ������[-lobefrac_est,lobefrac_est]

b_loc = fix(1.2*1.1*(n/BB_loc));
b_est = fix(1.4*1.1*(n/BB_est));

W_Comb = floor_to_pow2(Comb_cst*n/B_loc);      % ����˲�������B
WITH_COMB = false;                             % ����SFTT2.0
ALGORITHM1 = true;                             % loc_loops=0;��ʹ�û���˲���

% word=sprintf('\n\nRUNNING EXPERIMENT: n=%d, k=%d.\n', n, k);
% disp(word);                                  % ���n,k
% word=sprintf('\n\nSimulation:\n');
% disp(word);
% word=sprintf('*******************************************************\n');
% disp(word);

%% SFFT����Ƶ��
%% ��Ͱ��������Ͱ������Ǵ���������ѭ������
if B_thresh > B_loc ||threshold_loops >loc_loops 
    return;
end
%% ����ƽ̹���˲���
[filtert,w_loc ]=make_dolphchebyshev_t(lobefrac_loc, tolerance_loc);                 % ���� filtert �б�ѩ���˲�����w_loc������
[filter_time filter_sizet filter_freq]=make_multiple_t(filtert, w_loc, n, b_loc);  % ����ƽ̹���˲�����filter_time���˲���ʱ�� filter_sizet���˲������� filter_freq���˲���Ƶ��
% figure;plot(abs(filter_time))
% figure;plot(abs(filter_freq))
% aaaaa=fft(filtert);
% figure;plot(abs(aaaaa(1:8000)));title('dolphchebyshev filter time domain - window function');
% figure;plot(abs(filter_time));title('dolphchebyshev filter time domain - flat window');
% figure;plot(abs(filter_freq(1:8000)));title('dolphchebyshev filter freq domain - flat window');
[filtert_est w_est]=make_dolphchebyshev_t(lobefrac_est, tolerance_est);
[filter_est_time filter_est_sizet filter_est_freq]=make_multiple_t(filtert_est, w_est, n, b_est);
% word=sprintf(' Window size: Location Filter : %d; Estimation Filter : %d;\n', w_loc, w_est);
% disp(word);

filter_noise = 0;                       % ƽ̹���˲�����й¶
filter_noise_est = 0;
for i=0:9
    filter_noise = max(filter_noise,max(abs(filter_freq(n/2+i+1)),abs(filter_freq(n/2-i+1))));
    filter_noise_est = max(filter_noise_est,max(abs(filter_est_freq(n/2+i+1)),abs(filter_est_freq(n/2-i+1))));
end
% word=sprintf('Noise in filter: Location Filter : %d; Estimation Filter %d\n', filter_noise, filter_noise_est);
% disp(word);
% word=sprintf('****************************************************************\n\n');
% disp(word);
% word=sprintf('sFFT Results\n**************************************************************');
% disp(word);
%% ��ѭ��
loops=loc_loops+est_loops;          % �ܵ���ѭ������
for ii=1:repetitions
    % ��ѭ����ÿ�εõ�heavyֵ��λ�úʹ�С
    tic;
    [I,out_I,out]= outer_loop(x, n,filter_time,filter_sizet, filter_freq, filter_est_time,...
        filter_est_sizet,filter_est_freq, B_est, B_thresh, B_loc, W_Comb, ...
        Comb_loops, threshold_loops, loc_loops, loops , WITH_COMB, ALGORITHM1);
    toc
    % I     ��ֵ��λ�ã�����
    % out_I ��Ӧ�Ĺ�ֵ
    % out   hitsλ�ü���Ӧ�Ĺ�ֵ
end
%%
num_candidates= length(out_I);      % �����ƵĴ�ֵ����
x_f_Large = zeros(1,n);                 
out_Large = zeros(1,n);
candidates=zeros(2,num_candidates);
counter=0;
%% �ù��ƽ�������źŵ�Ƶ��
for i=1:num_candidates                    % ���㷨���ƵĽ�������ڶ�ά���󣬵ڶ�����λ�ã���һ������Ӧ�ķ�ֵ
    counter=counter+1;
    candidates(1,counter)=abs(out_I(i));  % ��ֵ�ķ�������
    candidates(2,counter)=I(i);           % ��ֵ��λ������
end
temp=sortrows((candidates).').';          % ��ֵ����ֵ��������
for i=1:k                                 % ȡK�����ķ�ֵ������ź�Ƶ��out_Large
    key=temp(2,num_candidates - k+ i);
    out_Large(key+1) =out(key+1)*n;       
end

end


%% �����������2��n����
function [out] = floor_to_pow2 (in)
    out = 1;
    while out <= in
        out = out*2;
    end
    out=out/2;
end


%% �����б�ѩ���˲���
function [ filter ,w] = make_dolphchebyshev_t(lobefrac,tolerance)
%% ����
% lobefrac ͨ����Ƶ
% tolerance й¶����
%% ����
w =fix((1 / pi) * (1/lobefrac) * acosh(1/tolerance)); %  ����  fix():��0ȡ��
if ~mod(w,2)==1                                       % ~����
    w=w-1;						                      % ��֤w���������˲����ǶԳƵ�
end
%% �˲���
filter = zeros(1,w);
t0 = cosh(acosh(1/tolerance) / (w-1));
for i=0:(w-1)
    filter(i+1) = Cheb(w-1, t0*cos(pi*i/w))*tolerance; % �б�ѩ���˲�����ʱ��
end
%figure;plot(abs(filter));title('win');
temp=fft(filter, w);			                     % �б�ѩ���˲�����Ƶ��
%figure;plot(real(temp));title('FREwin');
filter=fftshift(temp);                               % ����Ƶ���Ƶ�Ƶ�׵��м䣨��fft����֮���pi-2pi���ְ�����-pi-0���Ӷ�ʹ��Ƶ��������Ƶ�׵�����λ�á���
%figure;plot(real(temp));title('FREwin_half');

for i=1:w
    filter(i) = real(filter(i));                	 % ���ص���ʵ��
end
end

%% �����б�ѩ�򴰺������ɵ�ƽ̹������
function [out_time,out_sizet,out_freq]=make_multiple_t(filtert, w, n, b)
%% ����
% filtert �б�ѩ���˲�����Ƶ��
% w �˲����Ĵ���
% n �źų���
% b �˲�������
% out_time  ƽ̹���˲�����ʱ��
% out_sizet ƽ̹���˲�������
% out_freq  ƽ̹���˲�����Ƶ��
%%  �˲����Ĵ�����������������n
if b >= n || w >= n
    return;
end
%% 
g=zeros(1,n);
h=zeros(1,n);
w1=fix(w/2);
g(1:w-w1)=filtert((w1+1):w);
g(n-w1+1:n)=filtert(1:w1);
% figure;plot(abs(g));title('win');
g=fft(g,n);
s = 0;
for i=1:b
    s = g(i)+s ;									% �˲�������
end
offset = fix(b/2);
for i=0:(n-1)
    h(mod((i+offset),n)+1) = s;
    s = s + (g(mod((i+ b),n)+1) - g(i+1));
end

h_max=max(abs(h));
for i=1:n
    h(i) =h(i)/ h_max;
end
offsetc = 1;
step=exp(-2j*pi*w1/n);
for i=1:n
    h(i)= h(i)*offsetc;
    offsetc= offsetc *step;							% ��Ƶ���϶�ÿ����ֵ������λ
end
g=ifft(h,n);
filtert(1:w)=g(1:w);

out_time=filtert;								    % ʱ���ź�
out_sizet = w;										% �˲���Ƶ��
out_freq= h;										% Ƶ���ź�
end

%% ϡ�踵��Ҷ�任����ѭ��
function [I,out_I,out]=outer_loop(origx, n, filter_time,filter_sizet, filter_freq, ...
    filter_est_time,filter_est_sizet, filter_est_freq, B2,num,B,...
    W_Comb, Comb_loops, loop_threshold, location_loops, loops,  WITH_COMB, ALGORITHM1)
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
permutea =zeros(1,loops);
permuteb =zeros(1,loops);
x_samp=zeros(loops,max(B,B2));
score =zeros(1,n);
hits_ini=zeros(1,n);
hits_found_ini=0;
% a_INI=[3181497 3978157 3467391 2445571 1072031 2119695 1002225 3676519 3995001 ...
%     40375 1777721 610525 2569107 1769277 4168711];
% a_INI=[2029 10759 4377  12349  15373  4097  1123  8219  3823  13727  15029];

for i=1:loops
    %% �趨����
    a = 0;                                  % ������Ų���sigma��ģ��
    b = 0;							        % ����ʱ���±�ƽ�Ʋ��� ���൱�ڿα��ϵ�-a*sigma��
    while(undivide(a, n)~= 1 || mod(a,2)==0 || mod(a,5)==0)
        a = mod(fix(rand()*n) ,n);			% a��n������������ѭ��
    end
%     a
%     a=a_INI(i);
    ai = mod_inverse(a, n);			        % ģ�� a*ai mod n = 1
    permutea(i)=ai;						    % ����ÿ��ѭ����������Ų���ai
    permuteb(i) = b;					    % ����ÿ��ѭ����������Ų��� b
 
    perform_location =(location_loops>=i);   
    if perform_location                     % ���ö�λѭ��������
        cur_filter_time=filter_time;
        cur_filter_sizet=filter_sizet;
        cur_filter_freq=filter_freq;
        cur_B=B;
    else                                     % ���ù�ֵѭ��������
        cur_filter_time=filter_est_time;
        cur_filter_sizet=filter_est_sizet;
        cur_filter_freq=filter_est_freq;
        cur_B=B2;
    end
   
    %% ��Ͱ��
    [x_samp(i,1:cur_B),J ]=inner_loop_locate(origx,n,cur_filter_time,cur_filter_sizet,cur_filter_freq, ...
        num, cur_B,	a,ai, b);              
    % x_samp�����������Ƶ��ֵ��Ͱֵ��
    % J����Ͱ���±꼯
    % figure;plot(abs(fft(x_samp(i,1:cur_B))))��
    %%  �ۼƴ�ֵԭʼλ�ü���
     if  perform_location   
        [hits,hits_found,score_out ]=inner_loop_filter_regular(J, n, num, cur_B, a, ai, b, loop_threshold, score,hits_ini, hits_found_ini);
        hits_ini=hits;                  % �ҵ������еĴ�ֵλ�ü�
        hits_found_ini=hits_found;      % �ҵ��Ĵ�ֵλ������
        score=score_out;                % ��Ƿ����� 
    end
end
% word=sprintf('Number of candidates: %d\n', hits_found);
% disp(word);
%% ��ֵ
[I,out_I,out] = estimate_values(hits, hits_found, x_samp,  loops, n, permutea, ...
    B,B2,filter_time,filter_sizet, filter_freq, filter_est_time, ...
    filter_est_sizet, filter_est_freq,  location_loops);

end

%% �б�ѩ�򴰺���
function [out]= Cheb ( m, x)
if abs(x) <= 1
    out=cos(m * acos(x)); % acos(x):������arccos(x) 
else
    out=real(cosh(m * acosh(x))); % cosh(x):˫�����Һ���=(e^x+e^-x)/2
end
end

function [I,out_I,out] =estimate_values(hits, hits_found,x_samp,loops, n, ...
    permutea, B, B2,	filter_time,filter_sizet, filter_freq, filter_est_time, ...
    filter_est_sizet, filter_est_freq, location_loops)
%% ����
% hits  ��ֵƵ��λ�ü�
% hits_found ��ֵƵ�ʸ���
% x_samp ���������Ƶ��ֵ��Ͱֵ��
% loops ��ѭ������
% n �źų���
% permutea ����ÿ��ѭ����������Ų���ai
% B loc��B 
% B2  est��B
% filter_time ƽ̹���˲�����ʱ��
% filter_sizet ƽ̹���˲�������
% filter_freq ƽ̹���˲�����Ƶ��
% location_loops locѭ������
%output:
% I     ��ֵ��λ�ã�����
% out_I ��Ӧ�Ĺ�ֵ
% out   ��ֵ����λ�õĺ���
%% ��ֵ
out=zeros(1,hits_found);
I=zeros(1,hits_found);
out_val=zeros(2,loops);                      % ��¼��ѭ���Ĺ�ֵ�������һ��real���ڶ���imag
location = fix((loops - 1) / 2);                                  
for i=0:(hits_found-1)                       % ��i����ֵ��λ��
    r_loops=0; 
    for j=0:(loops-1)                        % �Ե�j����ѭ����ֵ
        %% ƽ̹���˲�����������
        if j < location_loops                                      
            cur_filter_time=filter_time;
            cur_filter_sizet=filter_sizet;
            cur_filter_freq=filter_freq;
            cur_B=B;
        else
            cur_filter_time=filter_est_time;
            cur_filter_sizet=filter_est_sizet;
            cur_filter_freq=filter_est_freq;
            cur_B=B2;
        end
        %% ��Ͱ�����Ͱֵ����ϣ��ƫ�ƺ���
        permuted_index= timesmod(permutea(j+1), hits(i+1), n);		% ���ź���±�
        hashed_to = fix(permuted_index / (n / cur_B));				% hash function �ֵ���Ͱ�±�
        dist = fix(rem(permuted_index , (n / cur_B)));				% offest function Ͱ��ƫ����
        if dist > ((n/cur_B)/2)                                     % ��ƫ���������˲�����ͨ����Ƶ��Ҫ�ֵ���һ��Ͱ��
            hashed_to = rem((hashed_to + 1),cur_B);
            dist =dist- n/cur_B;
        end
        %% ����Ƶ��ֵ
        dist = rem((n - dist) , n);              % ��Ͱʱ��Ӧ���˲����±�
        filter_value = cur_filter_freq(dist+1);  % ��Ͱʱ���˲���ֵ
%         x_samp(j+1,hashed_to+1) 
        out_val(1,r_loops+1) =real(x_samp(j+1,hashed_to+1) / filter_value); % Ƶ��ֵʵ��
        out_val(2,r_loops+1) =imag(x_samp(j+1,hashed_to+1) / filter_value); % Ƶ��ֵ�鲿
        r_loops=r_loops+1;     
    end
    %% �����й�ֵ��ȡ��ֵ
    for ii=1:2
        out_val(ii,:)=sort(out_val(ii,1:r_loops)); % ��λ��i�����й�ֵ����������
    end
    realv = out_val(1,location+1);                 % ʵ��ȡ��ֵ
    imagv = out_val(2,location+1);                 % �鲿ȡ��ֵ
    out(hits(i+1)+1) = realv + 1j*imagv;           % ��i����ֵ���յ�Ƶ�ʹ�ֵ
    if hits(i+1)~=0
        I(i+1)=hits(i+1);                          % ��i����ֵ��λ��                   
    end
end
I=sort(I);                                         % ��λ�ð���������
out_I=out(I+1);                                    % ������Ӧ��Ƶ��ֵ
I=[0,I];
out_I=[0,out_I];

end

%%   �ҵ���Ͱ�����±걣����out������
function [out]=find_largest_indices( num, samples,cur_B)
%% ����
% samples   ������
% num   ��Ͱ���� 
% cur_B  ��Ͱ����
% out ��ֵƵ����ܴ��ڵ����꼯
%% 
count = 0;
% out=zero(1,num);
cutoff = nth_element_immutable(samples, cur_B, cur_B-num);  %��Ͱ�Ĺ�����ֵ
% �����ϸ������ֵ�Ĵ�Ͱ
for i=1:cur_B
    if samples(i) > cutoff      
        count=count+1;
        out(count) = i-1;
    end
end
% ȷ���ҵ�num����Ͱ
if count < num						% ���ҵ��ĵ���С��B_thresh����num��������������һЩ��
    for i=1:cur_B
        if samples(i) == cutoff            
            count=count+1;
            out(count) = i-1;		% ��Ӧͼ2�е���Ӱ����
            if count >= num 
                break;
            end
        end
    end
    out=sort(out);                  % ������������
end
end

function [x_samp ]=inner_loop_est(x,filter_est_time,filter_est_sizet, ...
   cur_B,cur_i,cur_index,x_sampt)

for i=1:filter_est_sizet
    x_sampt(cur_i(i)) = x(cur_index(i))*filter_est_time(i)+x_sampt(cur_i(i));			%�˲�����ʱ��������Ϊ�����ʽ����Ӧz��ʱ��ֵ
%     x_sampt(cur_i(i)) = x(cur_index(i))+x_sampt(cur_i(i));			%�˲�����ʱ��������Ϊ�����ʽ����Ӧz��ʱ��ֵ
end

x_samp=fft(x_sampt,cur_B);


end

function [hits, hits_found,score_out]=inner_loop_filter_regular(J, n, num, cur_B, a, ...
    ai, b,loop_threshold, score,hits_ini, hits_found_ini)
%% ����
 % J ��Ͱ�±꼯
 % n �źų���
 % num ��Ͱ����
 % cur_B ��Ͱ��
 % a ai��ģ��
 % ai ����ʱ���±������������൱�ڿα��ϵ�sigma��
 % b ����ʱ���±�ƽ�Ʋ��� ���൱�ڿα��ϵ�-a*sigma��
 % loop_threshold ��ֵ��Ǵ�����ֵ
 % score  ��Ǵ�����
 % hits_ini %���ҵ��Ĵ�ֵλ�ü�
 % hits_found_ini ���ҵ��Ĵ�ֵλ�ø���
%output
 % hits �ҵ������еĴ�ֵλ�ü�
 % hits_found �ҵ��Ĵ�ֵλ������
 % score_out ��Ƿ����� 
 
%%
hits=hits_ini;
hits_found=hits_found_ini;
for i=1:num         
    % �ҷֵ���ͰJ(i)��ԭʼλ�ü� �����б��
    low  = mod(ceil((J(i)-0.5)*n/cur_B + n), n);	    % ceil����:�����������ȡ�� 
    high = mod(ceil((J(i)+0.5)*n/cur_B + n), n);		% �ֵ�ͰJ(i)��Ƶ��λ�ü�[low,high]
    loc = timesmod(low, a, n);                          % ���ź�λ��low ��Ӧԭʼ�źŵ�λ��loc=a*low-a*b
    j = low;
    while  j ~= high                                    % �Էֵ�ͰJ(i)��ԭʼƵ��λ�ñ��һ��
        score(loc+1)=score(loc+1)+1;                   
        if  score(loc+1)==loop_threshold                % ������Ǵ����ﵽ��ֵ������λ�ô���hits
            hits_found =hits_found+1;
            hits(hits_found)=loc;
        end
        loc = mod((loc + a),n);
        j = mod((j + 1),n);
    end
end

score_out=score;

end


%% ��λ��ѭ��
function [x_samp,J ]=inner_loop_locate(origx, n, cur_filter_time,cur_filter_sizet, ...
    cur_filter_freq,num,cur_B,  a, ai,b)
%% ����
% origx ԭʼ�ź�
% n �źų���
% cur_filter_time ƽ̹���˲�����ʱ��
% cur_filter_sizet ƽ̹���˲�������
% cur_filter_freq ƽ̹���˲�����Ƶ��
% num ��Ͱ����
% cur_B ��Ͱ��
% a ai��ģ��
% ai ����ʱ���±������������൱�ڿα��ϵ�sigma��
% b ����ʱ���±�ƽ�Ʋ��� ���൱�ڿα��ϵ�-a*sigma��
% x_samp ��Ͱֵ
% J ��num����Ͱ���꼯
%%  Ҫ���źų����Ƿ�Ͱ�������������������Warning
if mod(n,cur_B)       
    word=sprintf('Warning: n is not divisible by cur_B, which algorithm expects.\n');
    disp(word);
end
%% ��Ͱ����
x_sampt =zeros(1,n);  
index=b;
for i=0:(cur_filter_sizet-1)
    x_sampt(mod(i,cur_B)+1) = origx(index+1) * cur_filter_time(i+1) + x_sampt(mod(i,cur_B)+1); 
    index =rem((index+ai),n);	
    % origx(index+1)��������ţ���������index = b+ai*i,�� ai=sigma,b=-a*sigma
    % origx(index+1) * cur_filter_time(i+1)���˲���ʱ����ˣ�Ƶ������
    % x_sampt(mod(i,cur_B)+1) =....�� Ƶ�ʽ���������Ӧ��ʱ��ֵ��
end
x_samp=fft(x_sampt,cur_B);                    % ������Ƶ���ϵĲ���ֵ

%% Ѱ�Ҵ�Ͱλ�ü�
samples = zeros(1,cur_B);                     % ��������Ĺ�����
for i=1:cur_B
    samples(i) = (abs(x_samp(i)))^2;							
end
% figure;plot(samples);
J=find_largest_indices( num, samples ,cur_B);  % num����Ͱ������
end



%% ��ģ�� a*out mod n = 1
function [out]=mod_inverse(a, n)
i = n;
out = 0;
d = 1;
while a>0
    t = fix(i/a);
    x = a;
    a = mod(i,x);
    i = x;
    x = d;
    d = out - t*x;
    out = x;
end
out=rem(out,n);         % rem(x,y):������x/y������
if (out<0)
    out = rem((out+n),n);  
end
end


%% ������ɢ�źŵĵ�n��Ԫ��
function [out,coordinate]=nth_element(a,b)
len=length(a);
for i=1:len
    swap=0;
    for j=1:len-1
        if a(j)>a(j+1)
            temp1=a(j);
            a(j)=a(j+1);
            a(j+1)=temp1;
            temp2=b(j);
            b(j)=b(j+1);
            b(j+1)=temp2;
            swap=1;
        end
    end
    if ~swap
        break
    end
end
out=a;
coordinate=b;

end


%% ��ȡ��Ͱ�Ĺ�����ֵ
function [out]= nth_element_immutable(input,  cur_B, num1)
% input ��������Ĺ�����
% cur_B  ��Ͱ��
% num1   �Ǵ�ֵͰ����Ŀ
temp=input;
temp=sort(temp); % sort(X)����X��Ԫ�ؽ�����������
out=temp(num1);  % ��Ͱ����С����ֵ
end



function [out]=undivide(a,b)

if mod(a,b)==0
    out=b;
else
    temp=undivide(b, mod(a,b));
    out=temp;
end

end

%% ����ǰλ��out,���ź�λ��x
function [out]=timesmod(x, a, n) 
out=fix(rem(x*a,n));      % rem:ȡ�ຯ����fix:��0����ȡ��;
end
