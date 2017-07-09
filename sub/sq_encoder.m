function S = sq_encoder(R, B, M)
% Vector quantization using transformed coding and entropy-coded uniform
% scalar quantizer. 
%
% Author: Junting Chen (eejtchen@connect.ust.hk)
% Date: Feb 21, 2017.
% 
% INPUT
%   B       Required total entropy
%   R       Covariance matrix of the input source
%   M       (Optional) The M-most-significant dimension to quantize
%
% OUTPUT
%   S.type  'complex' source or 'real' source
%   S.D     total distortion
%   S.U     Transformation matrix, y = U'x
%   S.Bd    Cell array, each entry contains an array for the decision
%           boundary of a complex dimension
%   S.Xc    Cell array, each entry contains an array for the reconstruction
%           point per complex dimension

N = size(R, 1);
if nargin < 3
    M = N;
else
    M = max(1, min(N, M));
end

if norm(R' - R.', 'fro') < 1e-10
    S.type = 'real';
else
    S.type = 'complex';
end

[U, D] = eig(R);
d = real(diag(D));
[~, I] = sort(d, 'descend');
Sigma2 = d(I);
U = U(:, I);

S.U = U;
S.Sigma2 = Sigma2;

eps = B * 1e-2;
if strcmp(S.type, 'complex') || 1
    
%     % Method 1: Reverse waterfilling bit allocation
%     mu_max = log(2) * Sigma2(1);
%     mu_min = 0;
%     bits = zeros(1, N);
%     
%     MAXLOOP = 100;
%     cnt = 0;
%     
%     while abs(sum(bits) * 2 - B) > eps && cnt < MAXLOOP
%         cnt = cnt + 1;
%         
%         mu = (mu_min + mu_max) / 2;
%         for i = 1:N
%             bits(i) = max(0, log2(Sigma2(i) * log(2) / mu) / 2);
%         end
%         
%         if sum(bits) > B / 2
%             mu_min = mu;
%         else
%             mu_max = mu;
%         end
%     end
    
    % Method 2: Bit allocation via dimension selection
    bits2 = zeros(1, N);
    m_max = min(find(Sigma2 > 0, 1, 'last'), M);
    m_min = 1;
    while m_max - m_min > 0
        m = round((m_max + m_min) / 2);
        rho = geo_mean(Sigma2(1:m));
        bits2(1:m) = ((B / m) + log2(Sigma2(1:m) / rho)) / 2;
        if bits2(m) < 0
            m_max = m;
        else
            if m == m_max
                m_min = m;
            else
                bits2(m + 1) = ((B / m) + log2(Sigma2(m + 1) / rho)) / 2;
                if bits2(m + 1) > 0
                    m_min = m;
                else
                    m_min = m;
                    m_max = m;
                end
            end
        end 
    end
    bits2(m + 1: end) = 0;
    
    S.Bits = bits2;
    bits = bits2;
    
    S.Bd = cell(N, 1);
    S.Xc = cell(N, 1);
    D_vec = zeros(1, N);
    for i = 1:N
        if bits(i) > 0
            [S.Bd{i}, S.Xc{i}, ~, D_real] = unifsq_gauss(bits(i), Sigma2(i) / 2);
            D_vec(i) = D_real * 2;
        else
            D_vec(i) = Sigma2(i);
        end
    end
    
    S.D = sum(D_vec);
    S.Ds = D_vec;
    
else
    %
end



