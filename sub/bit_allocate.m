function Bs = bit_allocate(Alpha, Dims, Btot)
% Bs = bit_allocate(Alpha, Dims, Btot)
%
% Simulation code of the algorithm developed in the following paper
%
% J. Chen, H. Yin, L. Cottatellucci, and D. Gesbert, "Feedback Mechanisms 
% for FDD Massive MIMO with D2D-based Limited CSI Sharing", IEEE Trans. 
% Wireless Commun., 2017.
%
% Author: Junting Chen (eejtchen@connect.ust.hk)
% Date: Feb 21, 2017.

K = size(Alpha, 1);

bmax = Btot;
bmin = - Btot;

b = mean([bmax, bmin]);
Bs = zeros(K);

maxloop = 100;
i = 0;
while bmax - bmin > 1e-4 && i < maxloop
    i = i + 1;
    for k = 1:K
        for j = 1:K
            if k == j
                continue;
            end
        
            if Dims(k, j) <= 1
                Bs(k, j) = 0;
            else
                Bs(k, j) = (Dims(k, j) - 1) ...
                    * max(0, b + log2(Alpha(k, j) * log(2) / (Dims(k, j) - 1)));
            end
        end
    end
    Bs_tot = sum(Bs(:));
    
    if Bs_tot > Btot
        bmax = b;
    else
        bmin = b;
    end
    
    b = mean([bmax, bmin]);
end