function [B, Bs] = bit_allocate_equal(Btot, K)
% Bs = bit_allocate_equal(Btot, K)
%
% Author: Junting Chen (eejtchen@connect.ust.hk)
% Date: Feb 21, 2017.

N = K * (K - 1);
b_vec = ones(1, N) * floor(Btot / N);

r = Btot - sum(b_vec);
I = randperm(N);
b_vec(I(1:r)) = floor(Btot / N) + 1;

cnt = 0;
Bs = zeros(K);
for i = 1:K
    for j = 1:K
        if i == j
            continue;
        end
        cnt = cnt + 1;
        Bs(i, j) = b_vec(cnt);
    end
end

B = max(Bs, [], 2);