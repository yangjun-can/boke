%% ��ֵ��Ĺ��Ʋ���
function [I,out_I,out] =estimate_values(hits, hits_found,x_samp,loops, n, ...
    permutea, B, B2,filter_freq, filter_est_freq, location_loops)
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
location = fix((loops - 1) / 2); %��ֵ��λ��
% out=zeros(1,hits_found);
% I=zeros(1,hits_found);
% out_val=zeros(2,loops);                      % ��¼��ѭ���Ĺ�ֵ�������һ��real���ڶ���imag
% for i=0:(hits_found-1)                       % ��i����ֵ��λ��
%     r_loops=0;
%     for j=0:(loops-1)                        % �Ե�j����ѭ����ֵ
%         %% ƽ̹���˲�����������
%         if j < location_loops
%             %             cur_filter_time=filter_time;
%             %             cur_filter_sizet=filter_sizet;
%             cur_filter_freq=filter_freq;
%             cur_B=B;
%         else
%             %             cur_filter_time=filter_est_time;
%             %             cur_filter_sizet=filter_est_sizet;
%             cur_filter_freq=filter_est_freq;
%             cur_B=B2;
%         end
%         %% ��Ͱ�����Ͱֵ����ϣ��ƫ�ƺ���
%         %      permuted_index= timesmod(permutea(j+1), hits(i+1), n);		% ���ź���±�
%         permuted_index= mod(permutea(j+1)* hits(i+1), n);
%         hashed_to = fix(permuted_index / (n / cur_B));				% hash function �ֵ���Ͱ�±�
%         dist = fix(rem(permuted_index , (n / cur_B)));				% offest function Ͱ��ƫ����
%         if dist > ((n/cur_B)/2)                                     % ��ƫ���������˲�����ͨ����Ƶ��Ҫ�ֵ���һ��Ͱ��
%             hashed_to = rem((hashed_to + 1),cur_B);
%             dist =dist- n/cur_B;
%         end
%         %% ����Ƶ��ֵ
%         dist = rem((n - dist) , n);              % ��Ͱʱ��Ӧ���˲����±�
%         filter_value = cur_filter_freq(dist+1);  % ��Ͱʱ���˲���ֵ
%         %         x_samp(j+1,hashed_to+1)
%         out_val(1,r_loops+1) =real(x_samp(j+1,hashed_to+1) / filter_value); % Ƶ��ֵʵ��
%         out_val(2,r_loops+1) =imag(x_samp(j+1,hashed_to+1) / filter_value); % Ƶ��ֵ�鲿
%         r_loops=r_loops+1;
%     end
%     %% �����й�ֵ��ȡ��ֵ
%     out_val=sort(out_val(:,1:r_loops),2); % ��λ��i�����й�ֵ����������
%     realv = out_val(1,location+1);                 % ʵ��ȡ��ֵ
%     imagv = out_val(2,location+1);                 % �鲿ȡ��ֵ
%     out(hits(i+1)+1) = realv + 1i*imagv;           % ��i����ֵ���յ�Ƶ�ʹ�ֵ
%     if hits(i+1)~=0
%         I(i+1)=hits(i+1);                          % ��i����ֵ��λ��
%     end
% end
% I=sort(I);                                         % ��λ�ð���������
% out_I=out(I+1);                                    % ������Ӧ��Ƶ��ֵ
% I=[0,I];
% out_I=[0,out_I];

%% ƽ̹���˲�����������
cur_B=[repmat(B,hits_found,location_loops),repmat(B2,hits_found,loops-location_loops)];
%% ��Ͱ�����Ͱֵ����ϣ��ƫ�ƺ���
permuted_index= mod(hits(1:hits_found).'*permutea , n);           % ���ź���±�
hashed_to = fix(permuted_index .* cur_B /n );				% hash function �ֵ���Ͱ�±�
dist = fix(rem(permuted_index , (n ./ cur_B)));				% offest function Ͱ��ƫ����
% ��ƫ���������˲�����ͨ����Ƶ��Ҫ�ֵ���һ��Ͱ��
pand = dist > (n./cur_B)/2 ;
hashed_to(pand) = rem((hashed_to(pand) + 1),cur_B(pand));
dist(pand) =dist(pand)- n./cur_B(pand);
%% ����Ƶ��ֵ
dist = rem((n - dist) , n);              % ��Ͱʱ��Ӧ���˲����±�
dist1=dist(:,1:location_loops);
dist2=dist(:,location_loops+1:loops);
filter_value_1=[filter_freq(dist1+1).';filter_est_freq(dist2+1).'].';
% out_valy(1,j+1) =real(x_samp(j+1,hashed_to(i+1,j+1)+1) / filter_value(i+1,j+1)); % Ƶ��ֵʵ��
out_val=zeros(2,hits_found,loops);
filter_value=reshape(filter_value_1,1,[]);
hashed_to1=reshape(hashed_to,1,[]);
out_val1=real(x_samp(:,hashed_to1+1) ./ filter_value); % Ƶ��ֵʵ��
out_val2=imag(x_samp(:,hashed_to1+1) ./ filter_value); % Ƶ��ֵ�鲿
for jj=1:loops
    out_val(1,:,jj)=out_val1(jj,(jj-1)*hits_found+1:jj*hits_found);
    out_val(2,:,jj)=out_val2(jj,(jj-1)*hits_found+1:jj*hits_found);
end
%% �����й�ֵ��ȡ��ֵ
out_val=sort(out_val,3); % ��λ��i�����й�ֵ����������
realv = out_val(1,:,location+1);                 % ʵ��ȡ��ֵ
imagv = out_val(2,:,location+1);                 % �鲿ȡ��ֵ
out(hits((1:hits_found))+1) = realv + 1i*imagv;  % ��i����ֵ���յ�Ƶ�ʹ�ֵ
ton= hits~=0;  
I(ton)=hits(ton);                                % ��i����ֵ��λ��
I=sort(I);                                       % ��λ�ð���������
out_I=out(I+1);                                  % ������Ӧ��Ƶ��ֵ
I=[0,I];
out_I=[0,out_I];
end



