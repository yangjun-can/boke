%% �����������2��n����
function [out] = floor_to_pow2 (in)
%     out = 1;
%     while out <= in
%         out = out*2;
%     end
%     out=out/2;

% out=2^(nextpow2(in)-1);
out=2^ floor(log2(in));
end