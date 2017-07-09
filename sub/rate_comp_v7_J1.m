function [Rate, IntLk, Sig, Par] = rate_comp_v7_J1(B_fb, B_ms, B_ms0, ...
                                                 Sint, Uint, Rue, Gain, Cr, T)
%
% Simulation code of the algorithm developed in the following paper
%
% J. Chen, H. Yin, L. Cottatellucci, and D. Gesbert, "Feedback Mechanisms 
% for FDD Massive MIMO with D2D-based Limited CSI Sharing", IEEE Trans. 
% Wireless Commun., 2017.
%
% Author: Junting Chen (eejtchen@connect.ust.hk)
% Date: Feb 21, 2017.
%
% Version 7-J1: Quantization error awared precoding
% Version 6-J1: Use entropy-coded scalar quantization
% Version 5-J1: Consider delayed CSI exchange. Add an argument 'Cr' to
% specify the time-correlation for delayed CSI exchange. Add an argument T
% to specify the number of Monte-Carlo simulations
%
% Version 2-J1: Modification for the 1st revision of the journal paper
% titled "Efficient Feedback Mechanisms for FDD Massive MIMO under
% User-level Cooperation". The main change is to disable some unused
% schemes to speed up the simulation.
%
% Version 2: Using per-group transform quantization for D2D CSIf sharing 
% for both of the baseline (equal bit allocation) and the adaptive scheme.
%
% [Rate, IntLk, Sig] = rate_comp(B_fb, B_ms, B_ms0, Sint, Uint, Rue, Gain)
%
% Compare the sum rate over different precoding schemes:
%   - (i) MF at base station
%   - (ii) Regularized ZF at the base station
%   - (iii) SLNR precoding at the terminals
%   - (iv) SLNR precoding at the terminals using equal bit allocation
%   - (v) SLNR percoding at the terminals with perfect global CSI
% 
% INPUT
%   - B_fb      scaler, the per-user feedback bits to the BS
%   - B_ms      matrix, B_ms(k, j) the number of bits to quantize h_k to
%               user j
%   - B_ms0     Vector, B_ms(k) the number of bits to quantize h_k to
%               share to others
%   - Sint      K*K cell, Sint{k, j} containts a vector of eigenvalues of the
%               interference subspace, i.e., projecting R_k onto the subspace of R_j
%   - Uint      K*K Cell, Uint{k, j} containts the eigenvectors of the
%               interference subspace, R_k projected onto R_j
%   - Rue       K cell, Rue{k} is the covariance matrix Rk of user k
%   - Gain      K vector, Power gain vector
%
% OUTPUT
%   - Rate      N_scheme * K matrix, data rate
%   - IntLk     N_scheme * K matrix, interference leakage
%   - Sig
% 
% UPDATE
%   - Codebook design on the overlapping subspace
%   - Quantize the beamformer to be fed back to the BS

bTEST = 0;

N_schemes = 5;
% 1 - MF at BS
% 2 - RZF at BS
% 3 - SLNR D2D proposed
% 4 - SLNR D2D equal bit allocation
% 5 - SLNR D2D perfect CSI

K = size(Rue, 1);
Nt = size(Rue{1}, 1);

% Bmax = max(B_ms(:));
% T = min([10000, 2^Bmax, 2^B_fb]);
% T = 1000;

% TEST
if bTEST
    W31_stat = zeros(Nt);   % Statistics of UE 1's precoder under the proposed scheme
    W11_stat = zeros(Nt);   % Statistics of UE 1's precoder under MF
    W01_ideal = zeros(Nt);   % The ideal statistics of UE 1's precoder under SLNR  
end
% END TEST

Rm_cell = cell(1, K);
Amp = zeros(1, K);
for k = 1:K
    R = Rue{k};
    Rm = sqrtm(R);
    Rm_cell{k} = Rm;
    
    for j = 1:K
        if k == j
            continue;
        end

    end
    Amp(k) = sqrt(Gain(k));
end

% Codebook to feedback the channel / beamformer to the BS
% The codebook adapts to the covariance of the direct link channel
Cbs = cell(K, 1);
Mbf = 2^B_fb;
% R0m = sqrtm(R0);  % This is the case for using the same beamformer
                    % codebook for all the users
for k = 1:K
    C = zeros(Mbf, Nt);
    for i = 1:Mbf
        h = (randn(Nt, 1) + 1i * randn(Nt, 1)) / sqrt(2);
        h = Rm_cell{k} * h;
%         h = R0m * h;
        h = h / norm(h);
        C(i, :) = h.';
    end
    Cbs{k} = C;
end

% Channel codebook for D2D CSI sharing under equal bit allocation
% (Baseline)
Cm2m0 = cell(K, 1);
for k = 1:K
    bk = B_ms0(k);
    % Cm2m0{k} = gencodebookT(Rue{k}, bk);
    Cm2m0{k} = sq_encoder(Rue{k}, bk);
end

% Channel Codebook for D2D CSI sharing (prospoed)
Cm2m = cell(K);
for k = 1:K
    for j = 1:K
        if B_ms(k, j) == 0
            continue;
        end
        b_kj = B_ms(k, j);
        % Cm2m{k, j} = gencodebookT(diag(Sint{k, j}), b_kj);
        Cm2m{k, j} = sq_encoder(diag(Sint{k, j}), b_kj);
    end
end

% Matrices for quantization errors
QU = cell(K, 1);
for k = 1:K
    Qumat = zeros(Nt);
    for j = 1:K
        if B_ms(j, k) == 0
            continue;
        end
        % Qumat = Qumat + Qmat(j, k) * Uint{j, k} * Uint{j, k}';
        Qumat = Qumat + Gain(j) * Uint{j, k} * ...
                        (Cm2m{j, k}.U * diag(Cm2m{j, k}.Ds) * Cm2m{j, k}.U) * ...
                        Uint{j, k}';
    end
    QU{k} = Qumat;
end

% Pre-compute the regularization due to quantization error 
Cq = cell(K, 1);
for k = 1:K
    Cqmat = zeros(1, Mbf);
    for i = 1:Mbf
        Cqmat(i) = real(Cbs{k}(i, :) * QU{k} * Cbs{k}(i, :)');
    end
    Cq{k} = Cqmat;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- Evaluations ---------------------------------
Sig = zeros(N_schemes, K);
IntLk = zeros(N_schemes, K);
Rate = zeros(N_schemes, K);
for t = 1:T
    P = 1 / K;  % Per-user power (power coefficient + pathloss goes to the channel)
    
    % Channel
    H = zeros(K, Nt);
    H_delay = zeros(K, Nt);
    for k = 1:K
        h = (randn(Nt, 1) + 1i * randn(Nt, 1)) / sqrt(2);
        h = sqrt(Gain(k)) * Rm_cell{k} * h;
        H(k, :) = h.';
        
        h_delay = Cr * h + sqrt(1 - Cr^2) * Rm_cell{k} * randn(Nt, 1);
        H_delay(k, :) = h_delay.';
    end
    
    % --------------- Precoding at the BS (Baselines) --------------------
    % Channel quantization 
    Hq = zeros(K, Nt);
    for k = 1:K
        hk = H(k, :);
        gk = hk / norm(hk);
        Cbs_k = abs(Cbs{k} * gk');
        [~, I] = max(Cbs_k);
        gkq = Cbs{k}(I, :);
        Hq(k, :) = gkq * norm(hk);
    end
%     % - Match filtering
%     W10 = Hq';
%     powers = diag(W10' * W10);
%     W1 = W10 * diag(1 ./ sqrt(powers));
    % Regularized ZF precoding
    W20 = Hq' * inv(Hq * Hq' + eye(K) * 1 / P);
    powers = diag(W20' * W20);
    W2 = W20 * diag(1 ./ sqrt(powers));

    % --------------- Precoding at the User Terminals --------------------
    % Precoding at the User Terminals adaptive bit allocation
    % Quantization for DELAYED D2D CSI sharing (proposed)
    Hq_cell = cell(K);
    for k = 1:K
        hk = H_delay(k, :).';   % < ---- Delayed CSI for D2D
        for j = 1:K
            if isempty(Cm2m{k, j})
                continue;
            end
            U_kj = Uint{k, j};
            h_kj = U_kj' * hk;
           
            % h_kj_hat = quantizationT(h_kj, Cm2m{k, j});
            h_kj_hat = sq_decoder(h_kj / Amp(k), Cm2m{k, j}) * Amp(k);
            % h_hat = sq_decoder(h_kj / Amp(k), Cm2m{k, j});
            % h_kj_hat = h_hat / norm(h_hat) * norm(h_kj);
            
            Hq_cell{k, j} = U_kj * h_kj_hat;
        end 
    end
    % SLNR precoding
    W3 = zeros(Nt, K);
    for k = 1:K
        % Reconstructing the H
        Hqk = zeros(K, Nt);
        Hqk(1, :) = H(k, :);    % The 1st row is the direct link channel
        cnt = 1;
        for j = 1:K
            if ~isempty(Hq_cell{j, k})
                cnt = cnt + 1;
                Hqk(cnt, :) = Hq_cell{j, k}.';
            end
        end
        Hqk = Hqk(1:cnt, :);
        
        Cbs_Hk = abs(Hqk * Cbs{k}').^2;
        sum_Cbs_Hk = sum(Cbs_Hk, 1);
        Intf_k = sum_Cbs_Hk - Cbs_Hk(1, :);
        Sig_k = Cbs_Hk(1, :);
        Slnr_k = Sig_k ./ (1 / P + Intf_k + Cq{k});
        [~, I] = max(Slnr_k);
        W3(:, k) = Cbs{k}(I, :)';
    end
    
    % Precoding at the user terminals (equal bit allocation)
	% Quantization (Baseline)
    Hq0 = zeros(K, Nt);
    for k = 1:K
        hk = H_delay(k, :).';   % < ---- Delayed CSI for D2D
        % hk_hat = quantizationT(hk, Cm2m0{k});
        hk_hat = sq_decoder(hk / Amp(k), Cm2m0{k}) * Amp(k);
        % h_hat = sq_decoder(hk / Amp(k), Cm2m0{k});
        % hk_hat = h_hat / norm(h_hat) * norm(hk);
        Hq0(k, :) = hk_hat.';
    end
    % SLNR precoding?
    W4 = zeros(Nt, K);
    for k = 1:K
        Hqk0 = Hq0;
        Hqk0(k, :) = H(k, :);    % Direct link channel perfect
        
        Cbs_Hk = abs(Hqk0 * Cbs{k}').^2;
        sum_Cbs_Hk = sum(Cbs_Hk, 1);
        Sig_k = Cbs_Hk(k, :);
        Intf_k = sum_Cbs_Hk - Sig_k;
        Slnr_k = Sig_k ./ (1 / P + Intf_k);
        [~, I] = max(Slnr_k);
        W4(:, k) = Cbs{k}(I, :)';
    end
    
%     % Precoding at the user terminals (perfect CSI)
% 	% Perfect CSI
%     % SLNR precoding
%     W5 = zeros(Nt, K);
%     for k = 1:K
%         Cbs_Hk = abs(H * Cbs{k}').^2;
%         sum_Cbs_Hk = sum(Cbs_Hk, 1);
%         Sig_k = Cbs_Hk(k, :);
%         Intf_k = sum_Cbs_Hk - Sig_k;
%         Slnr_k = Sig_k ./ (1 / P + Intf_k);
%         [~, I] = max(Slnr_k);
%         W5(:, k) = Cbs{k}(I, :)';
%     end
    
    % TEST
    if bTEST
        W11_stat = W11_stat + W1(:, 1) * W1(:, 1)' / T;
        W31_stat = W31_stat + W3(:, 1) * W3(:, 1)' / T;
        
        Rin1 = eye(Nt) * (1 / P);
        for k = 2:K
            hk = (randn(Nt, 1) + 1i * randn(Nt, 1)) / sqrt(2);
            hk = sqrt(Gain(k)) * Rm_cell{k} * hk;
            hk1 = Uint{k, 1} * Uint{k, 1}' * hk;
            Rin1 = Rin1 + hk1 * hk1';
        end
        h1 = H(1, :).';
        Rins1 = Rin1 \ (h1 * h1');
        [U_Rins1, D_Rins1] = eig(Rins1);
        [~, I_Rins1] = max(real(diag(D_Rins1)));
        w1_ideal = U_Rins1(:, I_Rins1);
        W01_ideal = W01_ideal + w1_ideal * w1_ideal' / T;
    end
    % END TEST
    
    % -------------------------- Evaluations ---------------------------- %
    % SINR and rate evaluation
%     Y1 = abs(H * W1 * sqrt(P)).^2;
    Y2 = abs(H * W2 * sqrt(P)).^2;
    Y3 = abs(H * W3 * sqrt(P)).^2;
    Y4 = abs(H * W4 * sqrt(P)).^2;    
%     Y5 = abs(H * W5 * sqrt(P)).^2;
    
%     % BS MF
%     S1 = diag(Y1);
%     I1 = sum(Y1, 2) - S1;
%     sinr1 = S1 ./ (1 + I1);
%     r1 = log2(1 + sinr1);
    
    % BS RZF
    S2 = diag(Y2);
    I2 = sum(Y2, 2) - S2;
    sinr2 = S2 ./ (1 + I2);
    r2 = log2(1 + sinr2);
    
    % D2D proposed
    S3 = diag(Y3);
    I3 = sum(Y3, 2) - S3;
    sinr3 = S3 ./ (1 + I3);
    r3 = log2(1 + sinr3);
%     r3 = log2(1 + S3);  % TEST
    
    % D2D naive (equal bit)
    S4 = diag(Y4);
    I4 = sum(Y4, 2) - S4;
    sinr4 = S4 ./ (1 + I4);
    r4 = log2(1 + sinr4);
    
%     % D2D perfect CSI
%     S5 = diag(Y5);
%     I5 = sum(Y5, 2) - S5;
%     sinr5 = S5 ./ (1 + I5);
%     r5 = log2(1 + sinr5);
% %     r5 = log2(1 + S5);  % TEST

%     Sig(1, :) = Sig(1, :) + S1.' / T;
    Sig(2, :) = Sig(2, :) + S2.' / T;
    Sig(3, :) = Sig(3, :) + S3.' / T;
    Sig(4, :) = Sig(4, :) + S4.' / T;
%     Sig(5, :) = Sig(5, :) + S5.' / T;
    
%     IntLk(1, :) = IntLk(1, :) + I1.' / T;
    IntLk(2, :) = IntLk(2, :) + I2.' / T;
    IntLk(3, :) = IntLk(3, :) + I3.' / T;
    IntLk(4, :) = IntLk(4, :) + I4.' / T;
%     IntLk(5, :) = IntLk(5, :) + I5.' / T;
    
%     Rate(1, :) = Rate(1, :) + r1.' / T;
    Rate(2, :) = Rate(2, :) + r2.' / T;
    Rate(3, :) = Rate(3, :) + r3.' / T;
    Rate(4, :) = Rate(4, :) + r4.' / T;
%     Rate(5, :) = Rate(5, :) + r5.' / T;
    
end

if bTEST
    Par.W11_stat = W11_stat;
    Par.W31_stat = W31_stat;
    Par.W01_ideal = W01_ideal;
else
    Par = 0;
end







          