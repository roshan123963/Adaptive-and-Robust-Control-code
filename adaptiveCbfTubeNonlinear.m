% adaptiveCbfTubeNonlinear: Adaptive PI servo-control with barrier augmentation on the reference model and a Lyapunov tube,
% for the lateral-directional aircraft model with nonlinear matched uncertainty Theta' * Phi(x_p).
%
% The script follows the notation of TCST_PI_note/TCST_PI_Sections.pdf (Sections II-VI). Matrix and signal names are the
% paper's symbols (A_ref, H_pi, DeltaH_max, ...); configuration quantities with physical units carry unit suffixes.
% All states and signals are in SI units (rad, rad/s); lateral acceleration N_y is in g.
%
% Cases:
%   1) Proposed: adaptation + barrier on the reference model + tube (Theorem 4)
%   2) No tube: same, with delta = 0
%   3) Infeasible tube: adaptation rate too small, so the feasibility condition (F) fails
%   4) Safety only: barrier of [2] on the plant with the nominal model, no adaptation
%   5) Adaptation only: MRAC with pi = 0, no barrier
%
%=================================
%                                  INPUTS
%=================================
% None (script). Edit the "Configuration" section to change the uncertainty, limits, gains, or cases.
%=================================
% OUTPUTS
%=================================
% Command-window table of the theorem checks for every case, and figures of the constrained quantities, the tube,
% the Lyapunov bound, the feasibility margin, tracking, and the adaptive parameters, plus comparison figures across the
% cases: peak usage of each limit, usage over time per case, tube feasibility (Proposed vs Infeasible tube), and errors.
%
%=================================
% EXAMPLE USAGE
%=================================
%
% adaptiveCbfTubeNonlinear

clc; clear; close all;

%% Configuration
dt_s = 1e-4;                          % RK4 step
finalTime_s = 40;
logEvery_nd = 10;                     % log every logEvery_nd steps

% Design model (Sections II-III): x_p = [beta; p_s; r_s], u = [delta_ail; delta_rud]
A_p = [-0.11794  0.00085  -1.0001;
       -7.0113  -1.4492    0.22059;
        6.3035   0.06511  -0.41172];
B_p = [ 0        0.015257;
       -7.9662   2.6875;
        0.60926 -2.3577];
C_reg = [ 0       1         0;
         -2.6049  0.018724  0.067695];           % y_reg = [p_s; N_y]
D_reg = [0  0;
         0  0.33698];
C_lim = [0 1 0;
         1 0 0];                                 % y_lim = [p_s; beta]

n_p = 3;  m = 2;  n = n_p + m;  q = 2;

A_nom = [zeros(m) C_reg; zeros(n_p, m) A_p];     % nominal design model
B = [D_reg; B_p];
B_ref = [-eye(m); zeros(n_p, m)];
B_v = [eye(m); zeros(n_p, m)];
B_pi = [B_v B];
E_p = [zeros(n_p, m) eye(n_p)];

K_lqr = lqr(A_nom, B, diag([1.025 1.0289 0 0 1.6021]), diag([1 0.49129]));
A_ref = A_nom - B * K_lqr;                       % known, Hurwitz
Q = eye(n);
P = lyap(A_ref', Q);                             % A_ref' * P + P * A_ref = -Q
P_inv = inv(P);

% True plant (unknown to the controller): A = A_ref + B * K_star (Assumption 1), plus Theta' * Phi(x_p) (Section VI)
K_star = K_lqr + 0.4 * [zeros(m) K_lqr(:, 3:5)];   % uncertainty does not change K_I, so K_I* = K_lqr,I
A = A_ref + B * K_star;
betaScale_rad = deg2rad(0.5);
yawRateScale_radps = deg2rad(5);
Phi = @(x_p) [betaScale_rad * tanh(x_p(1) / betaScale_rad); yawRateScale_radps * tanh(x_p(3) / yawRateScale_radps)];
J_Phi = @(x_p) [sech(x_p(1) / betaScale_rad)^2 0 0; 0 0 sech(x_p(3) / yawRateScale_radps)^2];
L_Phi = 1;                                       % ||Phi(a) - Phi(b)|| <= ||a - b||
Theta = [0.4  0.5;
         0.2 -0.4];                              % q x m; column i acts on input i

Omega_star = [K_star Theta'];
Omega_hat_0 = [K_lqr zeros(m, q)];
theta_max = 1.4;                                 % Assumption 2 (Section VI form)

% Operational limits (Section II)
u_max = deg2rad([3; 2]);                         % aileron, rudder
y_max = [deg2rad(18); deg2rad(0.5)];             % p_s, beta
z_max = [u_max; y_max];
z_min = -z_max;
alpha = diag([80 8 40 4]);
epsilon_nu = 1e-3;                               % smoothing in nu_i = sqrt(||theta_hat_i||^2 + epsilon_nu^2)

% Command: p_s pulses of 15 dps and N_y pulses. The N_y amplitude needs more sideslip and rudder than the limits allow,
% so the beta and rudder constraints stay active while it is commanded (0.005 g keeps every constraint inactive).
rollRateCmd_radps = deg2rad(15);
lateralAccelCmd_nd = 0.02;                       % N_y in g
y_cmd = @(t_s) [rollRateCmd_radps; lateralAccelCmd_nd] * ((t_s >= 0 && t_s < 10) - (t_s >= 20 && t_s < 30));

% Cases
gammaProposed_nd = 1e5;
cases = struct( ...
    'name',        {'Proposed', 'No tube', 'Infeasible tube', 'Safety only', 'Adaptation only'}, ...
    'isAdaptive',  {true,       true,      true,              false,         true}, ...
    'cbfMode',     {'ref',      'ref',     'ref',             'plant',       'none'}, ...
    'isTube',      {true,       false,     true,              false,         false}, ...
    'gamma',       {gammaProposed_nd, gammaProposed_nd, 1e4,  0,             gammaProposed_nd});
caseColor = {[0 0.25 0.95], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.49 0.18 0.56], [0.5 0.5 0.5]};

%% Offline checks of the assumptions
fprintf('Assumption 1: eig(A_ref) = %s\n', mat2str(eig(A_ref).', 3));
fprintf('Open-loop true plant (linear part) eig(A_p + B_p*Delta) = %s\n', mat2str(eig(A(m+1:n, m+1:n)).', 3));
priorBound_nd = norm(Omega_star - Omega_hat_0, 'fro');
fprintf('Assumption 2: ||Omega* - Omega_hat(0)||_F = %.3f < theta_max = %.3f : %s\n', priorBound_nd, theta_max, ...
    passFail(priorBound_nd < theta_max));
H_w = C_lim * B_p;
kappa = min(svd(K_lqr(:, 1:m)));                 % K_hat_I is frozen at K_lqr,I (Remark 5.1)
fprintf('Assumption 4: sigma_min(C_lim*B_p) = %.4f, kappa = sigma_min(K_I) = %.4f\n', min(svd(H_w)), kappa);
N = null(B_p')';                                 % N * B_p = 0
T = inv([C_lim; N]);
T_y = T(:, 1:m);  T_eta = T(:, m+1:end);
A_ref_pp = A_ref(m+1:n, m+1:n);
A_0 = N * A_ref_pp * T_eta;
fprintf('Assumption 6: eig(A_0) = %s, tzero(A_p, B_p, C_lim) = %s : %s\n', mat2str(eig(A_0).', 4), ...
    mat2str(tzero(ss(A_p, B_p, C_lim, 0)).', 4), passFail(all(real(eig(A_0)) < 0)));
p_bar = max(eig(E_p * P_inv * E_p'));

% A priori bound on x_p_ref from Lemma 7 (scalar A_0 here: c_0 = 1, lambda_0 = -A_0)
lambda_0 = -max(real(eig(A_0)));
r_ref = norm(T_y) * norm(y_max) + norm(T_eta) * (norm(N * A_ref_pp * T_y) * norm(y_max) / lambda_0);

%% Simulation
nSteps = round(finalTime_s / dt_s);
nLog = floor(nSteps / logEvery_nd) + 1;
results = cell(numel(cases), 1);

for c = 1:numel(cases)
    cfg = cases(c);
    gammas = [0 * ones(1, m), cfg.gamma * ones(1, n_p), cfg.gamma * ones(1, q)];   % Gamma_I = 0: K_hat_I frozen
    Gamma = diag(gammas);
    isAdapted = gammas > 0;
    if any(isAdapted)
        V_bar_0 = theta_max^2 / min(gammas(isAdapted));                         % theta_max^2 * lambda_max(Gamma^-1)
    else
        V_bar_0 = 0;
    end

    par = struct('A', A, 'B', B, 'B_ref', B_ref, 'B_v', B_v, 'B_pi', B_pi, 'A_ref', A_ref, 'P', P, 'P_inv', P_inv, ...
        'Q', Q, 'Gamma', Gamma, 'C_lim', C_lim, 'E_p', E_p, 'alpha', alpha, 'z_min', z_min, 'z_max', z_max, ...
        'Theta', Theta, 'Phi', Phi, 'J_Phi', J_Phi, 'L_Phi', L_Phi, 'p_bar', p_bar, 'epsilon_nu', epsilon_nu, ...
        'y_cmd', y_cmd, 'm', m, 'n', n, 'q', q);

    chi = [zeros(n, 1); zeros(n, 1); Omega_hat_0(:); V_bar_0];   % x(0) = x_ref(0) = 0 (Assumption 3)

    logData = struct('t_s', zeros(1, nLog), 'z', zeros(2*m, nLog), 'z_ref', zeros(2*m, nLog), 'delta', zeros(2*m, nLog), ...
        'u', zeros(m, nLog), 'pi_aug', zeros(2*m, nLog), 'x', zeros(n, nLog), 'x_ref', zeros(n, nLog), ...
        'V', zeros(1, nLog), 'V_bar', zeros(1, nLog), 'F', zeros(2*m, nLog), 'y_reg', zeros(m, nLog), ...
        'Omega_hat', zeros(m, n + q, nLog));
    worst = struct('refTight', -inf, 'plantBox', -inf, 'tube', -inf, 'lyap', -inf, 'F', inf, 'xpRef', 0, 'x', 0);
    divergedAt_s = NaN;
    iLog = 0;

    for k = 0:nSteps
        t_s = k * dt_s;
        [k1, d] = closedLoopRhs(t_s, chi, par, cfg);

        % Per-step checks of the theorem claims (continuous-time quantities at the start of each step)
        e = chi(1:n) - chi(n+1:2*n);
        Delta_Omega = reshape(chi(2*n+1:2*n+m*(n+q)), m, n + q) - Omega_star;
        V = e' * P * e + sum(sum((Delta_Omega(:, isAdapted) ./ gammas(isAdapted)) .* Delta_Omega(:, isAdapted)));
        V_bar = chi(end);
        if cfg.isTube
            worst.refTight = max(worst.refTight, max([d.z_ref - (z_max - d.delta); (z_min + d.delta) - d.z_ref]));
            worst.tube = max(worst.tube, max(abs(d.z - d.z_ref) - d.delta));
            worst.F = min(worst.F, min(d.F));
        end
        if cfg.isAdaptive
            worst.lyap = max(worst.lyap, V - V_bar);
        end
        worst.plantBox = max(worst.plantBox, max(abs(d.z) ./ z_max));
        worst.xpRef = max(worst.xpRef, norm(chi(n+m+1:2*n)));
        worst.x = max(worst.x, norm(chi(1:n)));

        if mod(k, logEvery_nd) == 0
            iLog = iLog + 1;
            logData.t_s(iLog) = t_s;  logData.z(:, iLog) = d.z;  logData.z_ref(:, iLog) = d.z_ref;  logData.delta(:, iLog) = d.delta;
            logData.u(:, iLog) = d.u;  logData.pi_aug(:, iLog) = d.pi_aug;  logData.x(:, iLog) = chi(1:n);  logData.x_ref(:, iLog) = chi(n+1:2*n);
            logData.V(iLog) = V;  logData.V_bar(iLog) = V_bar;  logData.F(:, iLog) = d.F;
            logData.y_reg(:, iLog) = C_reg * chi(m+1:n) + D_reg * (d.u + Theta' * Phi(chi(m+1:n)));
            logData.Omega_hat(:, :, iLog) = reshape(chi(2*n+1:2*n+m*(n+q)), m, n + q);
        end
        if k == nSteps
            break
        end

        k2 = closedLoopRhs(t_s + dt_s / 2, chi + dt_s / 2 * k1, par, cfg);
        k3 = closedLoopRhs(t_s + dt_s / 2, chi + dt_s / 2 * k2, par, cfg);
        k4 = closedLoopRhs(t_s + dt_s, chi + dt_s * k3, par, cfg);
        chi = chi + dt_s / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

        if any(~isfinite(chi)) || norm(chi(1:n)) > 1e3
            divergedAt_s = t_s;
            break
        end
    end
    fn = fieldnames(logData);
    for f = 1:numel(fn)
        if strcmp(fn{f}, 'Omega_hat')
            logData.Omega_hat = logData.Omega_hat(:, :, 1:iLog);
        else
            logData.(fn{f}) = logData.(fn{f})(:, 1:iLog);
        end
    end
    results{c} = struct('log', logData, 'worst', worst, 'divergedAt_s', divergedAt_s, 'V_bar_0', V_bar_0);
end

%% Theorem checks
fprintf('\nLimits: |u_bl| <= [%.1f %.1f] deg, |p_s| <= %.1f dps, |beta| <= %.2f deg\n', rad2deg(z_max));
fprintf('Lemma 7 a priori bound: ||x_p_ref|| <= r_ref = %.4f\n\n', r_ref);
for c = 1:numel(cases)
    r = results{c};  wc = r.worst;  lg = r.log;
    fprintf('=== Case %d: %s ===\n', c, cases(c).name);
    if ~isnan(r.divergedAt_s)
        fprintf('  DIVERGED at t = %.2f s\n', r.divergedAt_s);
    end
    delta0 = lg.delta(:, 1);
    fprintf('  delta(0) / half-width             : %s\n', mat2str((delta0 ./ z_max)', 3));
    fprintf('  Plant peak |z_i| / limit           : %s  (Thm 4.2: <= 1)  %s\n', ...
        mat2str(max(abs(lg.z), [], 2)' ./ z_max', 4), passFail(wc.plantBox <= 1 + 1e-9));
    fprintf('  Total command peak |u| / limit     : %s  (not guaranteed, Remark 1)\n', mat2str(max(abs(lg.u), [], 2)' ./ u_max', 4));
    if cases(c).isTube
        fprintf('  Ref. model in tightened box (Thm 2.1): max violation %.2e  %s\n', wc.refTight, passFail(wc.refTight <= 1e-9));
        fprintf('  Tube |z - z_ref| <= delta (Lem 3/5) : max violation %.2e  %s\n', wc.tube, passFail(wc.tube <= 1e-12));
        fprintf('  Feasibility (F) margin, min over t  : %.4g  %s\n', wc.F, passFail(wc.F > 0));
    end
    if cases(c).isAdaptive
        fprintf('  V <= V_bar (Lemma 5)               : max(V - V_bar) = %.2e  %s\n', wc.lyap, passFail(wc.lyap <= 1e-12));
    end
    if strcmp(cases(c).cbfMode, 'ref')
        fprintf('  ||x_p_ref|| <= r_ref (Lemma 7)     : %.4f <= %.4f  %s\n', wc.xpRef, r_ref, passFail(wc.xpRef <= r_ref));
    end
    eNorm = vecnorm(lg.x - lg.x_ref);
    fprintf('  max ||x|| = %.4f, ||e|| peak = %.2e, ||e(T)|| = %.2e  (Thm 4.3-4.4)\n\n', wc.x, max(eNorm), eNorm(end));
end

%% Figures
labels = {'u_{bl,ail}, deg', 'u_{bl,rud}, deg', 'p_s, dps', '\beta, deg'};
figure('Name', 'Constrained quantities z = [u_bl; y_lim]');
for i = 1:2*m
    subplot(2*m, 1, i); hold on; grid on;
    for c = 1:numel(cases)
        plot(results{c}.log.t_s, rad2deg(results{c}.log.z(i, :)), 'Color', caseColor{c}, 'LineWidth', 1.1);
    end
    yline(rad2deg([z_min(i) z_max(i)]), 'k-.');
    ylabel(labels{i});
    ylim(rad2deg(2.5 * [z_min(i) z_max(i)]));
end
xlabel('Time, s');
subplot(2*m, 1, 1); legend({cases.name}, 'Location', 'eastoutside');

figure('Name', 'Total command u = u_bl + w');
for i = 1:m
    subplot(m, 1, i); hold on; grid on;
    for c = 1:numel(cases)
        plot(results{c}.log.t_s, rad2deg(results{c}.log.u(i, :)), 'Color', caseColor{c}, 'LineWidth', 1.1);
    end
    yline(rad2deg([-u_max(i) u_max(i)]), 'k-.');
    ylabel(labels{i}(1:end-5));
    ylim(rad2deg(3 * [-u_max(i) u_max(i)]));
end
legend({cases.name}, 'Location', 'eastoutside'); xlabel('Time, s');

lg = results{1}.log;
figure('Name', 'Proposed: tube and reference model');
for i = 1:2*m
    subplot(2*m, 1, i); hold on; grid on;
    plot(lg.t_s, rad2deg(lg.z(i, :)), 'b', lg.t_s, rad2deg(lg.z_ref(i, :)), 'b--', 'LineWidth', 1.1);
    plot(lg.t_s, rad2deg(z_max(i) - lg.delta(i, :)), 'r--', lg.t_s, rad2deg(z_min(i) + lg.delta(i, :)), 'r--');
    yline(rad2deg([z_min(i) z_max(i)]), 'k-.');
    ylabel(labels{i});
end
xlabel('Time, s');
subplot(2*m, 1, 1); legend('z (plant)', 'z_{ref}', 'tightened limits', '', 'limits', 'Location', 'eastoutside');

figure('Name', 'Proposed: theorem quantities');
subplot(3, 1, 1);
plot(lg.t_s, rad2deg(abs(lg.z - lg.z_ref)), '-', 'LineWidth', 1.1); hold on; grid on;
set(gca, 'ColorOrderIndex', 1);
plot(lg.t_s, rad2deg(lg.delta), '--', 'LineWidth', 1.1);
ylabel('deg (dps)'); title('|z_i - z_{ref,i}| (solid) <= \delta_i (dashed)');
subplot(3, 1, 2);
semilogy(lg.t_s, lg.V, 'b', lg.t_s, lg.V_bar, 'r--', 'LineWidth', 1.1); grid on;
legend('V', 'V_{bar}'); title('Lemma 5: V <= V_{bar}');
subplot(3, 1, 3);
plot(lg.t_s, lg.F, 'LineWidth', 1.1); grid on; hold on; yline(0, 'k');
ylabel('(F) margin'); title('\alpha_i(z_{max}-z_{min}-2\delta_i) - 2 d\delta_i/dt > 0'); xlabel('Time, s');

figure('Name', 'Tracking');
for i = 1:m
    subplot(m, 1, i); hold on; grid on;
    for c = 1:numel(cases)
        scale = [180 / pi, 1];
        plot(results{c}.log.t_s, scale(i) * results{c}.log.y_reg(i, :), 'Color', caseColor{c}, 'LineWidth', 1.1);
    end
    cmd = arrayfun(@(t) y_cmd(t), lg.t_s, 'UniformOutput', false);
    cmd = [cmd{:}];
    plot(lg.t_s, scale(i) * cmd(i, :), 'k--');
end
subplot(m, 1, 1); ylabel('p_s, dps'); legend([{cases.name}, {'command'}], 'Location', 'eastoutside');
subplot(m, 1, 2); ylabel('N_y, g'); xlabel('Time, s');

figure('Name', 'Proposed: adaptive parameters');
Omega_log = results{1}.log.Omega_hat;
for i = 1:m
    for j = m+1:n+q
        subplot(m, n + q - m, (i - 1) * (n + q - m) + j - m); hold on; grid on;
        plot(lg.t_s, squeeze(Omega_log(i, j, :)), 'b', 'LineWidth', 1.1);
        yline(Omega_star(i, j), 'r--');
        title(sprintf('\\Omega_{%d,%d}', i, j));
    end
end
sgtitle('Adaptive parameters (red dashed: \Omega^*)');

%% Comparison figures
caseNames = {cases.name};
nCases = numel(cases);
shortLabels = {'u_{bl,ail}', 'u_{bl,rud}', 'p_s', '\beta'};
constraintColor = lines(2*m);
halfWidth = (z_max - z_min) / 2;
iProposed = find(strcmp(caseNames, 'Proposed'));
iInfeasible = find(strcmp(caseNames, 'Infeasible tube'));

% 1) Peak usage of each limit (Theorem 4.2 guarantees <= 1 for z; Remark 1: no guarantee for the total command u)
peakZ_nd = zeros(2*m, nCases);
peakU_nd = zeros(m, nCases);
for c = 1:nCases
    peakZ_nd(:, c) = max(abs(results{c}.log.z), [], 2) ./ z_max;
    peakU_nd(:, c) = max(abs(results{c}.log.u), [], 2) ./ u_max;
end
figure('Name', 'Comparison: peak usage of each limit');
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile;
plotPeakBars(peakZ_nd, caseColor, shortLabels, 'Constrained quantities z = [u_{bl}; y_{lim}] (guaranteed \leq 1)');
nexttile;
plotPeakBars(peakU_nd, caseColor, {'u_{ail}', 'u_{rud}'}, 'Total command u = u_{bl} + w (not guaranteed)');
legendHandle = legend(caseNames, 'Orientation', 'horizontal');
legendHandle.Layout.Tile = 'south';

% 2) Usage of each limit over time, one panel per case; the shaded band is outside the limit
figure('Name', 'Comparison: |z_i| / limit over time');
layoutHandle = tiledlayout(nCases, 1, 'TileSpacing', 'compact');
for c = 1:nCases
    nexttile; hold on; grid on;
    tCase_s = results{c}.log.t_s;
    usage_nd = abs(results{c}.log.z) ./ z_max;
    yTop_nd = max(1.5, 1.05 * max(usage_nd(:)));
    patch([tCase_s(1) tCase_s(end) tCase_s(end) tCase_s(1)], [1 1 yTop_nd yTop_nd], [1 0.85 0.85], 'EdgeColor', 'none');
    for i = 1:2*m
        plot(tCase_s, usage_nd(i, :), 'Color', constraintColor(i, :), 'LineWidth', 1.1);
    end
    yline(1, 'k-.');
    ylim([0 yTop_nd]);
    title(sprintf('%s: peak %.2f', caseNames{c}, max(usage_nd(:))), 'Color', caseColor{c});
end
legendHandle = legend([{'outside limit'}, shortLabels]);
legendHandle.Layout.Tile = 'east';
xlabel(layoutHandle, 'Time, s');
ylabel(layoutHandle, '|z_i| / z_{max,i}');

% 3) Tube feasibility: Proposed vs Infeasible tube
tubeCases = [iProposed iInfeasible];
tubeStyle = {'-', '--'};
figure('Name', 'Comparison: tube feasibility, Proposed vs Infeasible tube');
layoutHandle = tiledlayout(4, 1, 'TileSpacing', 'compact');
nexttile; hold on; grid on;
for iTube = 1:2
    lgTube = results{tubeCases(iTube)}.log;
    for i = 1:2*m
        plot(lgTube.t_s, lgTube.delta(i, :) / halfWidth(i), tubeStyle{iTube}, 'Color', constraintColor(i, :), 'LineWidth', 1.1);
    end
end
yline(1, 'k-.');
ylabel('\delta_i / half-width');
title('Tube width (solid: Proposed, dashed: Infeasible tube); \geq 1 leaves no room for z_{ref}');
legend(shortLabels, 'Location', 'eastoutside');

nexttile; hold on; grid on;
for iTube = 1:2
    lgTube = results{tubeCases(iTube)}.log;
    marginF_nd = min(lgTube.F ./ (diag(alpha) .* (z_max - z_min)), [], 1);
    plot(lgTube.t_s, marginF_nd, 'Color', caseColor{tubeCases(iTube)}, 'LineWidth', 1.1);
end
yline(0, 'k-.');
ylabel('min_i F_i / (\alpha_i w_i)');
title('Feasibility condition (F): margin must stay > 0 (w_i = z_{max,i} - z_{min,i})');
legend(caseNames(tubeCases), 'Location', 'eastoutside');

nexttile; hold on; grid on;
for iTube = 1:2
    lgTube = results{tubeCases(iTube)}.log;
    refViolation_nd = max(max(lgTube.z_ref - (z_max - lgTube.delta), (z_min + lgTube.delta) - lgTube.z_ref) ./ halfWidth, [], 1);
    plot(lgTube.t_s, refViolation_nd, 'Color', caseColor{tubeCases(iTube)}, 'LineWidth', 1.1);
end
yline(0, 'k-.');
ylabel('max_i dist / half-width');
title('Reference model vs tightened box: \leq 0 means inside (Theorem 2.1)');
legend(caseNames(tubeCases), 'Location', 'eastoutside');

nexttile; hold on; grid on;
for iTube = 1:2
    lgTube = results{tubeCases(iTube)}.log;
    plot(lgTube.t_s, max(abs(lgTube.z) ./ z_max, [], 1), 'Color', caseColor{tubeCases(iTube)}, 'LineWidth', 1.1);
end
yline(1, 'k-.');
ylabel('max_i |z_i| / z_{max,i}');
title('Plant vs original box: \leq 1 is safe (guaranteed only when (F) holds)');
legend(caseNames(tubeCases), 'Location', 'eastoutside');
xlabel(layoutHandle, 'Time, s');

% 4) Model-following error and parameter error
figure('Name', 'Comparison: model-following and parameter errors');
layoutHandle = tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile; hold on; grid on;
for c = 1:nCases
    eNorm = vecnorm(results{c}.log.x - results{c}.log.x_ref);
    plot(results{c}.log.t_s, max(eNorm, 1e-12), 'Color', caseColor{c}, 'LineWidth', 1.1);
end
set(gca, 'YScale', 'log');
ylim([1e-8 1]);
ylabel('||e|| = ||x - x_{ref}||');
title('Model-following error (Theorem 4.4: e \rightarrow 0 for the adaptive cases)');
nexttile; hold on; grid on;
for c = 1:nCases
    Omega_err = squeeze(sqrt(sum((results{c}.log.Omega_hat - Omega_star).^2, [1 2])));
    plot(results{c}.log.t_s, Omega_err, 'Color', caseColor{c}, 'LineWidth', 1.1);
end
yline(theta_max, 'k-.', '\theta_{max}');
yLimits = ylim();
ylim([yLimits(1) 1.05 * theta_max]);              % keep the theta_max label inside the axes
ylabel('||\Omega_{hat} - \Omega^*||_F');
title('Parameter error (bounded by Theorem 1; convergence to zero is not claimed)');
legendHandle = legend(caseNames);
legendHandle.Layout.Tile = 'east';
xlabel(layoutHandle, 'Time, s');

%% Local functions
function [chiDot, d] = closedLoopRhs(t_s, chi, par, cfg)
% closedLoopRhs: Closed-loop right-hand side (plant, reference model, adaptive law, V_bar) with the barrier augmentation
%
%=================================
%                                  INPUTS
%=================================
% t_s: time
% chi: closed-loop state [x; x_ref; vec(Omega_hat); V_bar]
% par: struct of model, controller, and limit parameters
% cfg: case configuration (isAdaptive, cbfMode = 'ref' | 'plant' | 'none', isTube)
%=================================
% OUTPUTS
%=================================
% chiDot: time derivative of chi
% d: diagnostics (z, z_ref, delta, F margin, u, pi_aug)
%
%=================================
% EXAMPLE USAGE
%=================================
%
% chiDot = closedLoopRhs(0, chi, par, cfg)
    m = par.m;  n = par.n;  q = par.q;
    x = chi(1:n);
    x_ref = chi(n+1:2*n);
    Omega_hat = reshape(chi(2*n+1:2*n+m*(n+q)), m, n + q);
    V_bar = chi(end);
    x_p = x(m+1:n);
    x_p_ref = x_ref(m+1:n);

    psi = [x; par.Phi(x_p)];
    psi_ref = [x_ref; par.Phi(x_p_ref)];
    e = x - x_ref;
    u_bl = -Omega_hat * psi;                                     % (48)
    if cfg.isAdaptive
        Omega_hat_dot = par.B' * par.P * e * psi' * par.Gamma;    % (49)
    else
        Omega_hat_dot = zeros(m, n + q);
    end
    eQe = e' * par.Q * e;
    V_bar_dot = -eQe;
    pi_0 = [-par.y_cmd(t_s); zeros(m, 1)];
    K_hat = Omega_hat(:, 1:n);
    Theta_hat = Omega_hat(:, n+1:n+q)';                           % q x m

    delta = zeros(2*m, 1);
    deltaDot = zeros(2*m, 1);
    switch cfg.cbfMode
        case 'none'
            pi_aug = zeros(2*m, 1);
        otherwise
            if strcmp(cfg.cbfMode, 'ref')
                xc = x_ref;  psic = psi_ref;                      % barrier on the reference model
            else
                xc = x;  psic = psi;                              % barrier of [2] on the plant with the nominal model
            end
            x_pc = xc(m+1:n);
            zeta = [-Omega_hat * psic; par.C_lim * x_pc];         % (52)
            G = [-K_hat - Theta_hat' * par.J_Phi(x_pc) * par.E_p; par.C_lim * par.E_p];   % (53)
            dzeta_dt = [-Omega_hat_dot * psic; zeros(m, 1)];      % (54)
            H_pi = G * par.B_pi;                                  % (55)
            ell = dzeta_dt + G * par.A_ref * xc + G * par.B_pi * pi_0 + par.alpha * zeta;   % (56)

            if cfg.isTube
                c_y = [zeros(m) par.C_lim];
                s_y = sum((c_y * par.P_inv) .* c_y, 2);
                delta_y = sqrt(s_y * V_bar);                      % (29)
                deltaDot_y = -s_y * eQe ./ (2 * delta_y);         % (30) with c_i constant

                K_hat_dot = Omega_hat_dot(:, 1:n);
                Theta_hat_dot = Omega_hat_dot(:, n+1:n+q)';
                s_u = sum((K_hat * par.P_inv) .* K_hat, 2);
                sDot_u = 2 * sum((K_hat_dot * par.P_inv) .* K_hat, 2);
                nu = sqrt(sum(Theta_hat.^2, 1)' + par.epsilon_nu^2);
                nuDot = sum(Theta_hat .* Theta_hat_dot, 1)' ./ nu;
                delta_u = sqrt(s_u * V_bar) + par.L_Phi * nu * sqrt(par.p_bar * V_bar);     % (57)
                deltaDot_u = (sDot_u * V_bar + s_u * V_bar_dot) ./ (2 * sqrt(s_u * V_bar)) ...
                    + par.L_Phi * sqrt(par.p_bar) * (nuDot * sqrt(V_bar) + nu * V_bar_dot / (2 * sqrt(V_bar)));   % (58)
                delta = [delta_u; delta_y];
                deltaDot = [deltaDot_u; deltaDot_y];
            end
            DeltaH_max = ell + deltaDot - par.alpha * (par.z_max - delta);   % (37)
            DeltaH_min = -ell + deltaDot + par.alpha * (par.z_min + delta);  % (38)
            pi_aug = H_pi \ (max(0, DeltaH_min) - max(0, DeltaH_max));      % (40)
    end
    v = pi_aug(1:m);
    w = pi_aug(m+1:2*m);
    u = u_bl + w;

    xDot = par.A * x + par.B * (u + par.Theta' * par.Phi(x_p)) + par.B_ref * par.y_cmd(t_s) + par.B_v * v;   % (47)
    x_refDot = par.A_ref * x_ref + par.B_pi * (pi_0 + pi_aug);                                        % (17)
    chiDot = [xDot; x_refDot; Omega_hat_dot(:); V_bar_dot];

    if nargout > 1
        d.z = [u_bl; par.C_lim * x_p];
        d.z_ref = [-Omega_hat * psi_ref; par.C_lim * x_p_ref];
        d.delta = delta;
        d.F = diag(par.alpha) .* (par.z_max - par.z_min - 2 * delta) - 2 * deltaDot;
        d.u = u;
        d.pi_aug = pi_aug;
    end
end

function str = passFail(isPass)
% passFail: Returns 'PASS' or 'FAIL' for a logical check
%
%=================================
%                                  INPUTS
%=================================
% isPass: logical result of a check
%=================================
% OUTPUTS
%=================================
% str: 'PASS' if isPass is true, otherwise 'FAIL'
%
%=================================
% EXAMPLE USAGE
%=================================
%
% passFail(1 < 2)
    if isPass
        str = 'PASS';
    else
        str = 'FAIL';
    end
end

function plotPeakBars(peak_nd, caseColor, tickLabels, titleText)
% plotPeakBars: Grouped bar chart of peak |signal| / limit, one group per limit and one bar per case, with the limit at 1
%
%=================================
%                                  INPUTS
%=================================
% peak_nd: nLimits x nCases peak ratios (nLimits >= 2)
% caseColor: 1 x nCases cell of RGB triplets
% tickLabels: 1 x nLimits cell of group labels
% titleText: axes title
%=================================
% OUTPUTS
%=================================
% None (draws into the current axes)
%
%=================================
% EXAMPLE USAGE
%=================================
%
% plotPeakBars([0.9 1.2; 0.3 0.5], {[0 0 1], [1 0 0]}, {'a', 'b'}, 'Peak usage')
    hold on; grid on;
    barHandle = bar(peak_nd, 'grouped');
    for c = 1:numel(barHandle)
        barHandle(c).FaceColor = caseColor{c};
        text(barHandle(c).XEndPoints, barHandle(c).YEndPoints, compose('%.2f', peak_nd(:, c)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
    yline(1, 'k-.', 'limit', 'LineWidth', 1.2);
    set(gca, 'XTick', 1:size(peak_nd, 1), 'XTickLabel', tickLabels);
    ylabel('peak |\cdot| / limit');
    title(titleText);
end
