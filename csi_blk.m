% This script generates simulation results depicted in Fig. 2 of the paper
%
% J. Chen, H. Yin, L. Cottatellucci, and D. Gesbert, "Feedback Mechanisms 
% for FDD Massive MIMO with D2D-based Limited CSI Sharing", IEEE Trans. 
% Wireless Commun., 2017.
%
% Author: Junting Chen (juntingc@usc.edu)
% Date: Feb 21, 2017.

clear
close all
addpath sub

Nmc = 1000000;    % Number of Monte-Carlo simulation for ecah configuration
% ========================= [Geographic Configuration] ====================
Pos_ue = [0 0];           % in meters
Distance_ms2bs = 60;       % in meters
AngularSpread_deg = 10;

% ========================= [Basic Configuration] =========================
Nt = 60;
outer_dim = 20;            % outer dimension for the two-layer precoding
    
SNR_dB = [];
SNR0_dB = 20;

Blockage_dB = 0:5:30;
Blockage0_dB = 0;

Bfb_vec = []; % # of bits per user to feedback to BS
Bfb0 = 6;

Bms0 = 80;
Bms_vec = [];   % Total # of bits for D2D CSI sharing
    
Vel_km = 5;                % user velocity
D2D_delay_ms = 10;         % delay for CSI sharing
fc = 2e9;                  % carrier frequency

N_schemes = 5;             % number of schemes (proposed + 4 baselines)
% -------------------------- TIME CORRRELATIONS  ------------------------ %
fd = Vel_km / 3.6 / 3e8 * fc;   % Maximum doppler frequency
Cr = besselj(0, 2 * pi * fd * D2D_delay_ms * 1e-3);

% ----------------------------- SUBSPACE PROCESSING --------------------- %
Nusers = length(Pos_ue);

Aod_ue_deg = Pos_ue / Distance_ms2bs * 180 / pi;
% Preprocessing
Rue0 = cell(Nusers, 1);
R0all = zeros(Nt);
for k = 1:Nusers
    Rue0{k} = genChannCorr(Nt, Aod_ue_deg(k), AngularSpread_deg);
    R0all = R0all + Rue0{k};
end

% The joint dominate subspace
[U, D] = eig(R0all);
[Ds, Ids] = sort(real(diag(D)), 'descend');
Ut = U(:, Ids(1:outer_dim));
Rue = cell(Nusers, 1);
for k = 1:Nusers
    Rue{k} = Ut' * Rue0{k} * Ut; % Equivalent channel
end
R0 = Ut' * R0all * Ut;

%%
% The pair-wise interference subspace
power_cutoff0 = 0.99;    % The power cutoff parameters are for implementation 
power_cutoff1 = 0.99;    % interest. Note that, the approximation works well
Uint = cell(Nusers);     % under well-conditioned covariance matrices
Rint = cell(Nusers);
Sint = cell(Nusers);
Sigma0 = zeros(1, Nusers);
Sigma_cell = cell(Nusers);
Dims = zeros(Nusers);
for j = 1:Nusers 
    [Uj0, Dj] = eig(Rue{j});
    Dj_vec = real(diag(Dj));
    Dj_vec_cutoff = power_truncate(Dj_vec, power_cutoff0);
    Uj = Uj0(:, Dj_vec_cutoff > 0);
    Pj = Uj * Uj';     % Projection matrix onto the subspace spanned by Uj
    
    R_j = zeros(outer_dim);
    for k = 1:Nusers
        if j == k
            continue;
        end
        Rkj = Pj * Rue{k} * Pj';    % Project Rk onto (dominate) subspace of Rj
        [Ukj0, Dkj] = eig(Rkj);
        Dkj_vec = real(diag(Dkj));
        Dkj_vec_cutoff = power_truncate(Dkj_vec, power_cutoff1);
        Ukj = Ukj0(:, Dkj_vec_cutoff > 0);
        Uint{k, j} = Ukj;
        Rint{k, j} = Rkj;
        Sint{k, j} = real(diag(Ukj' * Rkj * Ukj));
        Dims(k, j) = size(Ukj, 2);
        
        R_j = R_j + Rue{k} / trace(Rue{k}) * outer_dim;
    end
    
    [U_j0, D_j] = eig(R_j);
    D_j_vec = real(diag(D_j));
    D_j_vec_cutoff = power_truncate(D_j_vec, power_cutoff1);
    U_j = U_j0(:, D_j_vec_cutoff > 0);
    P_j = eye(outer_dim) - U_j * U_j';
    Rj0 = P_j * Rue{j} * P_j;
    Sigma0(j) = real(trace(Rj0)) / real(trace(Rue{j})) * outer_dim;
end

% ------------------------- END SUBSPACE PROCESSING --------------------- %
Nsnr = length(SNR_dB);
Nblk = length(Blockage_dB);
Nbbf = length(Bfb_vec);
Nbms = length(Bms_vec);

if Nusers ~= 2
    Nbms = 0;
end

Ncase = Nsnr + Nblk + Nbbf + Nbms;

SNRs = [10.^(SNR_dB / 10), ones(1, Nblk + Nbbf + Nbms) * 10^(SNR0_dB / 10)];
Blks = [ones(1, Nsnr) * 10^(- Blockage0_dB / 10), ...
        10.^(- Blockage_dB / 10), ...
        ones(1, Nbbf + Nbms) * 10^(- Blockage0_dB / 10)];
Bfs = [Bfb0 * ones(1, Nsnr + Nblk), Bfb_vec, Bfb0 * ones(1, Nbms)];
Bmss = [Bms0 * ones(1, Nsnr + Nblk + Nbbf), Bms_vec];

Rsum = zeros(N_schemes, Ncase);
Rsum_all = zeros(N_schemes * Nusers, Ncase);
Intsum = zeros(N_schemes, Ncase);
Bit1_vec = zeros(1, Ncase);

parfor i_case = 1:Ncase
    fprintf('Running case %d/%d ...\n', i_case, Ncase);
    
    Gain = SNRs(i_case) * [ones(1, Nusers - 1), Blks(i_case)];
    B_fb = Bfs(i_case);
    Btot_ms_i = Bmss(i_case);
    
    % Bit allocation coefficient
    % Bits(k, j) - how many bits to quantize h_k to user j (through
    %               projecting h_k onto the subspace of user j
    Alpha = zeros(Nusers);
    for k = 1:Nusers
        for j = 1:Nusers
            if k == j
                continue;
            end
            if Dims(k, j) > 0
                Alpha(k, j) = real(prod(Sint{k, j}))^(1 / Dims(k, j)) ...
                                * (1 - Sigma0(k) / outer_dim)...
                                * Gain(k);
            else
                Alpha(k, j) = 0;
            end
        end
    end

    % Bit allocation
  
    Bits_preset = [0,                      Bmss(i_case)
                   Bfb0 - Bmss(i_case),           0];    
    Bits_eq = bit_allocate_equal(Btot_ms_i, Nusers);
    Bits_prop = bit_allocate(Alpha, Dims, Btot_ms_i);
    Bit1_vec(i_case) = sum(Bits_prop(1, :));
    
    % Throughput evaluation
    [Rate, IntLk, ~, Par] = rate_comp_v7_J1(B_fb, round(Bits_prop), Bits_eq, ...
                                            Sint, Uint, Rue, Gain, Cr, Nmc);
    Rsum(:, i_case) = sum(Rate, 2);
    Rsum_all(:, i_case) = Rate(:);
    Intsum(:, i_case) = sum(IntLk, 2);
    
end

%% Figures
my_markers = {'d', 's', '^', 'o', '<', '>', 'v', '+', 'p', '.', '*'}.';
Alg_scheme_name = {
    'BS-MRC'
    'CSI Feedback'
    'Proposed'
    'Naive CSI Exch'
    'D2D Perfect CSI'
};

schemes_to_show = [3 4 2]; 
N_scheme_to_show = length(schemes_to_show);

% Throughput over SNRs
fig_i = 6;
if Nsnr > 1
    fig_i = fig_i + 1;
    figure(fig_i),
    p_handle = plot(SNR_dB, Rsum(schemes_to_show, 1:Nsnr), ...
                '-', ...
                'LineWidth', 1.5,...
                'MarkerSize', 9);
    set(p_handle, {'marker'}, my_markers(1:N_scheme_to_show));
    set(gca, 'FontSize', 14);
    ylim([4 14]);
    set(gca, 'YTick', 4:2:14);
    xlabel('Total TX Power (dB)');
    ylabel('Sum Rate (bps/Hz)');
    title(sprintf('%d Users, B_{fb} = %d per UE, B_{tot,m2m} = %d', Nusers, Bfb0, Bms0));
    legend(Alg_scheme_name{schemes_to_show},'Location', 'southeast');
    tune_figure,
end

% Throughput over blockage of the Second user
if Nblk > 1
    fig_i = fig_i + 1;
    figure(fig_i),
    p_handle = plot(Blockage_dB, Rsum(schemes_to_show, Nsnr + 1 : Nsnr + Nblk), ...
                '-', ...
                'LineWidth', 1.5,...
                'MarkerSize', 9);
    set(p_handle, {'marker'}, my_markers(1:N_scheme_to_show));
    set(gca, 'YTick', 5:9);
    xlabel('Blockage of the 2nd user (dB)');
    ylabel('Sum Rate (bps/Hz)');
    title(sprintf('Total TX Power = %2.0f dB', SNR0_dB));
    legend(Alg_scheme_name{schemes_to_show}, 'Location', 'southwest');
    tune_figure,
end

% Throughput over # of bits per user to feedback to the BS
if Nbbf > 0
    fig_i = fig_i + 1;
    figure(fig_i),
    p_handle = plot(Bfb_vec, Rsum(schemes_to_show, Nsnr + Nblk + 1 : Nsnr + Nblk + Nbbf), ...
                '-', ...
                'LineWidth', 1.5,...
                'MarkerSize', 9);
    set(p_handle, {'marker'}, my_markers(1:N_scheme_to_show));
    set(gca, 'FontSize', 14);
    ylim([7 13]);
    set(gca, 'XTick', 4:2:12);
    xlabel('Bits per UE to feedback to BS');
    ylabel('Sum Rate (bps/Hz)');
    title(sprintf('Total TX Power = %2.0f dB', SNR0_dB));
    legend(Alg_scheme_name{schemes_to_show}, 'Location', 'southeast');  
    tune_figure,
end

% % Blockage versus bit allocation on user 1
% if Nblk > 1
%     fig_i = fig_i + 1;
%     figure(fig_i),
%     p_handle = plot(Blockage_dB, Bit1_vec(Nsnr + 1 : Nsnr + Nblk), ...
%                 'r^-', ...
%                 'LineWidth', 1,...
%                 'MarkerSize', 8);
%     xlabel('Blockage of the 2nd user (dB)');
%     ylabel('Bits of user 1');
%     title(sprintf('Total TX Power = %2.0f dB', SNR0_dB));
%     legend('Adaptive bit allocation', 'Location', 'southwest');
% 
% end

% Throughput versus the total number of bits for D2D CSI sharing
if Nbms > 1
    fig_i = fig_i + 1;
    figure(fig_i),
    p_handle = plot(Bms_vec, Rsum(schemes_to_show, Nsnr + Nblk + Nbbf + 1:end), ...
                '-', ...
                'LineWidth', 1.5,...
                'MarkerSize', 9);
    set(p_handle, {'marker'}, my_markers(1:N_scheme_to_show));
    set(gca, 'FontSize', 14);
    ylim([7 11]);
    set(gca, 'YTick', 7:11);
    xlabel('Total number of bits for D2D CSI sharing');
    ylabel('Sum Rate (bps/Hz)');
    title(sprintf('%d Users, Total TX Power = %2.0f dB, B_{fb} = %d', Nusers, SNR0_dB, Bfb0));
    legend(Alg_scheme_name{schemes_to_show}, 'Location', 'southeast');
    tune_figure,
end
















