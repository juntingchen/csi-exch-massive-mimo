function V = power_truncate(V0, t)
% v = power_truncate(v0, t)
% 
% INPUT
%   v0      a vector
%   t       a threshold, keep the total output power just above the
%           threshold t * sum(v0)

[VS, I] = sort(V0, 'descend');
VS_sum = cumsum(VS);
I_cutoff = find(VS_sum > t * sum(V0), 1);
V = zeros(size(V0));
V(I(1:I_cutoff)) = V0(I(1:I_cutoff));