%% �ۼƴ�ֵλ�ü���ÿ��ѭ������λ�õĲ�����
function [hits, hits_found,score_out]=inner_loop_filter_regular(J, n, num,...
        cur_B, a,loop_threshold, score,hits_ini, hits_found_ini)
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
%     % �ҷֵ���ͰJ(i)��ԭʼλ�ü� �����б��
%     low  = mod(ceil((J(1:num)-0.5)*n/cur_B + n), n);	    % ceil����:�����������ȡ�� 
%     high = mod(ceil((J(1:num)+0.5)*n/cur_B + n), n);		% �ֵ�ͰJ(i)��Ƶ��λ�ü�[low,high]
%     loc = timesmod(low, a, n);                          % ���ź�λ��low ��Ӧԭʼ�źŵ�λ��loc=a*low-a*b
%     len = mod(high-low,n);
% %     for i=1:num
% %         score(sub2ind(size(score),[mod(loc(i)+(0:len)*a,n)+1]))=score(sub2ind(size(score),[mod(loc(i)+(0:len)*a,n)+1]))+1;      
% %     end        
%   for i=1:num  
%     j = low(i);
%     while  j ~= high(i)                                    % �Էֵ�ͰJ(i)��ԭʼƵ��λ�ñ��һ��
%         score(loc(i)+1)=score(loc(i)+1)+1;                   
%         if  score(loc(i)+1)==loop_threshold                % ������Ǵ����ﵽ��ֵ������λ�ô���hits
%             hits_found =hits_found+1;
%             hits(hits_found)=loc(i);
%         end
%         loc(i) = mod((loc(i)+ a),n);
%         j = mod((j + 1),n);
%     end
%   end

score_out=score;

end


