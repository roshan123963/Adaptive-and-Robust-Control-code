% adaptiveCbfTubeVsNoTube: Tube vs no tube for the adaptive PI servo-controller with barrier augmentation on the
% reference model, for the lateral-directional aircraft model with nonlinear matched uncertainty Theta' * Phi(x_p).
%
% Both cases use the same plant, adaptive law, and barrier (TCST_PI_note/TCST_PI_Sections.pdf, Sections III-VI).
% The only difference is the tube delta_i in (37)-(38):
%   Tube:    the reference model is kept in the box tightened by delta_i (Theorem 4)
%   No tube: delta_i = 0, so the reference model is kept in the original box
% During adaptation the plant and the reference model differ by e = x - x_ref, so z_i = z_ref,i + c_i * e. Without the
% tube, z_ref,i rides the original limit and the plant overshoots it by c_i * e. With the tube, |c_i * e| <= delta_i
% (Lemmas 5 and 9) and z_ref,i stays delta_i inside the limit, so the plant stays inside.
%
% The barrier rates alpha are faster than in adaptiveCbfTubeNonlinear (200 on both input rows, 40 on beta), so the
% reference model reaches the rudder limit during the command transients, when e is largest.
%
%=================================
%                                  INPUTS
%=================================
% None (script). Edit the "Configuration" section to change the uncertainty, limits, or gains.
%=================================
% OUTPUTS
%=================================
% Command-window summary of the peaks, the no-tube violations, and the theorem checks for the tube case, and figures:
% one per case, no tube vs tube, a zoom on the worst violation, the tube bound, and command tracking.
%
%=================================
% EXAMPLE USAGE
%=================================
%
% adaptiveCbfTubeVsNoTube

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
K_star = K_lqr + 0.4 * [zeros(m) K_lqr(:, 3:5)];
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

% Adaptation: Gamma_I = 0 keeps K_hat_I frozen (Remark 5.1)
gamma_nd = 1e5;
Gamma = diag([zeros(1, m), gamma_nd * ones(1, n_p + q)]);
V_bar_0 = theta_max^2 / gamma_nd;                % theta_max^2 * lambda_max(Gamma^-1) over the adapted entries

% Operational limits (Section II)
u_max = deg2rad([3; 2]);                         % aileron, rudder
y_max = [deg2rad(18); deg2rad(0.5)];             % p_s, beta
z_max = [u_max; y_max];
z_min = -z_max;
alpha = diag([200 200 40 40]);
epsilon_nu = 1e-3;                               % smoothing in nu_i = sqrt(||theta_hat_i||^2 + epsilon_nu^2)

% Command: p_s and N_y pulses; the N_y amplitude keeps the beta constraint active while it is commanded
rollRateCmd_radps = deg2rad(15);
lateralAccelCmd_nd = 0.02;                       % N_y in g
y_cmd = @(t_s) [rollRateCmd_radps; lateralAccelCmd_nd] * ((t_s >= 0 && t_s < 10) - (t_s >= 20 && t_s < 30));

caseNames = {'No tube', 'Tube'};
isTubeCase = [false true];
caseColor = {[0.6 0.6 0.6], [0 0 1]};           % no tube gray, tube blue

%% Offline checks of the assumptions
priorBound_nd = norm(Omega_star - Omega_hat_0, 'fro');
fprintf('Assumption 2: ||Omega* - Omega_hat(0)||_F = %.3f < theta_max = %.3f\n', priorBound_nd, theta_max);
N = null(B_p')';                                 % N * B_p = 0
M_inv = inv([C_lim; N]);
A_0 = N * A_ref(m+1:n, m+1:n) * M_inv(:, m+1:end);
fprintf('Assumption 6: eig(A_0) = %s (minimum phase)\n', mat2str(eig(A_0).', 4));
p_bar = max(eig(E_p * P_inv * E_p'));

%% Simulation
par = struct('A', A, 'B', B, 'B_ref', B_ref, 'B_v', B_v, 'B_pi', B_pi, 'A_ref', A_ref, 'P', P, 'P_inv', P_inv, ...
    'Q', Q, 'Gamma', Gamma, 'C_lim', C_lim, 'E_p', E_p, 'alpha', alpha, 'z_min', z_min, 'z_max', z_max, ...
    'Theta', Theta, 'Phi', Phi, 'J_Phi', J_Phi, 'L_Phi', L_Phi, 'p_bar', p_bar, 'epsilon_nu', epsilon_nu, ...
    'y_cmd', y_cmd, 'm', m, 'n', n, 'q', q);
isAdapted = diag(Gamma)' > 0;
nSteps = round(finalTime_s / dt_s);
nLog = floor(nSteps / logEvery_nd) + 1;
results = cell(1, 2);

for c = 1:2
    par.isTube = isTubeCase(c);
    chi = [zeros(n, 1); zeros(n, 1); Omega_hat_0(:); V_bar_0];   % x(0) = x_ref(0) = 0 (Assumption 3)
    logData = struct('t_s', zeros(1, nLog), 'z', zeros(2*m, nLog), 'z_ref', zeros(2*m, nLog), 'delta', zeros(2*m, nLog), ...
        'u', zeros(m, nLog), 'x', zeros(n, nLog), 'x_ref', zeros(n, nLog), 'V', zeros(1, nLog), 'V_bar', zeros(1, nLog), ...
        'F', zeros(2*m, nLog), 'y_reg', zeros(m, nLog), 'y_cmd', zeros(m, nLog), 'OmegaErr', zeros(1, nLog));
    iLog = 0;
    for k = 0:nSteps
        t_s = k * dt_s;
        [k1, d] = closedLoopRhs(t_s, chi, par);
        if mod(k, logEvery_nd) == 0
            iLog = iLog + 1;
            e = chi(1:n) - chi(n+1:2*n);
            Delta_Omega = reshape(chi(2*n+1:2*n+m*(n+q)), m, n + q) - Omega_star;
            logData.t_s(iLog) = t_s;
            logData.z(:, iLog) = d.z;
            logData.z_ref(:, iLog) = d.z_ref;
            logData.delta(:, iLog) = d.delta;
            logData.u(:, iLog) = d.u;
            logData.x(:, iLog) = chi(1:n);
            logData.x_ref(:, iLog) = chi(n+1:2*n);
            logData.V(iLog) = e' * P * e + sum(sum(Delta_Omega(:, isAdapted).^2 ./ diag(Gamma(isAdapted, isAdapted))'));
            logData.V_bar(iLog) = chi(end);
            logData.F(:, iLog) = d.F;
            logData.y_reg(:, iLog) = C_reg * chi(m+1:n) + D_reg * (d.u + Theta' * Phi(chi(m+1:n)));
            logData.y_cmd(:, iLog) = y_cmd(t_s);
            logData.OmegaErr(iLog) = norm(Delta_Omega, 'fro');
        end
        if k == nSteps
            break
        end
        k2 = closedLoopRhs(t_s + dt_s / 2, chi + dt_s / 2 * k1, par);
        k3 = closedLoopRhs(t_s + dt_s / 2, chi + dt_s / 2 * k2, par);
        k4 = closedLoopRhs(t_s + dt_s, chi + dt_s * k3, par);
        chi = chi + dt_s / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    end
    results{c} = logData;
end
lgNo = results{1};
lgTube = results{2};
t_s = lgTube.t_s;

%% Summary
zMax_deg = rad2deg(z_max);
fprintf('\n%-30s %8s %8s %8s %8s\n', 'peak |z_i|', 'dela', 'delr', 'ps', 'beta');
fprintf('%-30s %8.3f %8.3f %8.3f %8.3f\n', 'limit (deg, dps)', zMax_deg);
for c = 1:2
    fprintf('%-30s %8.3f %8.3f %8.3f %8.3f\n', caseNames{c}, rad2deg(max(abs(results{c}.z), [], 2)));
    fprintf('%-30s %8.3f %8.3f %8.3f %8.3f\n', '  peak / limit', max(abs(results{c}.z), [], 2) ./ z_max);
end

usageNo_nd = abs(lgNo.z) ./ z_max;
[peakNo_nd, iPeak] = max(usageNo_nd, [], 2);
[worstUsage_nd, iWorst] = max(peakNo_nd);
tWorst_s = t_s(iPeak(iWorst));
rowNames = {'dela', 'delr', 'ps', 'beta'};
fprintf('\nNo tube: worst violation %.1f%% of the %s limit at t = %.2f s\n', 100 * (worstUsage_nd - 1), rowNames{iWorst}, tWorst_s);
fprintf('  at that time ||Omega_hat - Omega*||_F = %.3f (initial %.3f): the adaptation has not converged\n', ...
    lgNo.OmegaErr(iPeak(iWorst)), lgNo.OmegaErr(1));
fprintf('  z_ref,i = %.3f (limit %.3f), c_i e = %+.3f (deg or dps)\n', rad2deg(lgNo.z_ref(iWorst, iPeak(iWorst))), ...
    zMax_deg(iWorst), rad2deg(lgNo.z(iWorst, iPeak(iWorst)) - lgNo.z_ref(iWorst, iPeak(iWorst))));
for i = 1:2*m
    isViolated = usageNo_nd(i, :) > 1;
    if any(isViolated)
        fprintf('  %-4s outside its limit for %.3f s in total\n', rowNames{i}, sum(isViolated) * dt_s * logEvery_nd);
    end
end

refDist = max(lgTube.z_ref - (z_max - lgTube.delta), (z_min + lgTube.delta) - lgTube.z_ref);
checkTol_nd = 1e-9;                              % numerical tolerance on the <= 0 checks
checkNames = {'Theorem 2.1, x_ref in the tightened box', 'Lemmas 5, 9, |c_i e| - delta_i', 'Lemma 5, V - V_bar', ...
    'Theorem 2.2, |z_i| - z_max,i (plant)', '(F), minus the margin'};
checkValues = [max(refDist(:)), max(max(abs(lgTube.z - lgTube.z_ref) - lgTube.delta)), max(lgTube.V - lgTube.V_bar), ...
    max(max(abs(lgTube.z) - z_max)), -min(lgTube.F(:))];
fprintf('\nTube case theorem checks (each <= 0):\n');
for iCheck = 1:numel(checkNames)
    checkStatus = 'FAIL';
    if checkValues(iCheck) <= checkTol_nd
        checkStatus = 'ok';
    end
    fprintf('  %-42s %-4s  %.2e\n', checkNames{iCheck}, checkStatus, checkValues(iCheck));
end

%% Figures
grayColor = caseColor{1};
ylab = {'dela, deg', 'delr, deg', 'ps, dps', 'beta, deg'};

% One figure per case, in the layout of Robust_Adaptive_TK.m
for c = 1:2
    lg = results{c};
    beta_deg = rad2deg(lg.x(m+1, :));
    ps_dps = rad2deg(lg.x(m+2, :));
    figure('Name', caseNames{c});
    subplot(2, 2, 1);
    plot(rad2deg(lg.u(1, :)), rad2deg(lg.u(2, :)), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); hold on; grid on;
    plot(rad2deg(lg.z(1, :)), rad2deg(lg.z(2, :)), 'b', 'LineWidth', 0.5);
    rectangle('Position', [-zMax_deg(1) -zMax_deg(2) 2 * zMax_deg(1) 2 * zMax_deg(2)], 'EdgeColor', 'r', 'LineWidth', 1);
    xlabel('dela (deg)'); ylabel('delr (deg)'); title('Control Operational Bounds');
    axis([-3 3 -3 3] .* zMax_deg([1 1 2 2])');  % u_bl + w can be much larger (Remark 1); keep the box readable
    legend('u_{bl} + w', 'u_{bl} (constrained)', 'Location', 'best');

    subplot(2, 2, 3);
    plot(ps_dps, beta_deg, 'LineWidth', 0.5); hold on; grid on;
    rectangle('Position', [-zMax_deg(3) -zMax_deg(4) 2 * zMax_deg(3) 2 * zMax_deg(4)], 'EdgeColor', 'r', 'LineWidth', 1);
    xlabel('ps (dps)'); ylabel('beta (deg)'); title('Limited Output Operational Bounds');

    rowsToPlot = [4 2 3];                        % beta, delr, ps
    for j = 1:3
        i = rowsToPlot(j);
        subplot(4, 2, 2 * j); hold on; grid on;
        if i == 3
            plot(t_s, rad2deg(lg.y_cmd(1, :)), 'k--', 'LineWidth', 1.0);
        end
        plot(t_s, rad2deg(lg.z(i, :)), 'b', t_s, rad2deg(lg.z_ref(i, :)), 'c--', 'LineWidth', 1.2);
        yline([-zMax_deg(i) zMax_deg(i)], 'k-.');
        if isTubeCase(c)
            plot(t_s, rad2deg(z_max(i) - lg.delta(i, :)), 'r--', t_s, rad2deg(z_min(i) + lg.delta(i, :)), 'r--');
        end
        isViolated = abs(lg.z(i, :)) > z_max(i);
        plot(t_s(isViolated), rad2deg(lg.z(i, isViolated)), 'r.', 'MarkerSize', 6);
        ylabel(ylab{i});
        if j == 1
            legend('z (plant)', 'z_{ref}', 'limit', 'Location', 'eastoutside');
        end
    end

    subplot(4, 2, 8);
    plot(t_s, vecnorm(lg.x - lg.x_ref), 'b', 'LineWidth', 1.2); grid on; ylabel('||e||');
    yyaxis right; plot(t_s, lg.OmegaErr, 'LineWidth', 1.0); ylabel('||\Omega_{hat} - \Omega^*||_F');
    xlabel('Time, sec');
    sgtitle(sprintf('MRAC + CBF, %s (red dots: outside the limit)', lower(caseNames{c})));
end

% No tube vs tube: constrained quantities z = [u_bl; y_lim]
figure('Name', 'No tube vs tube');
for i = 1:2*m
    subplot(2*m, 1, i); hold on; grid on;
    hNo = plot(t_s, rad2deg(lgNo.z(i, :)), 'Color', grayColor, 'LineWidth', 1.2);
    hTube = plot(t_s, rad2deg(lgTube.z(i, :)), 'b', 'LineWidth', 1.2);
    hLim = yline([-zMax_deg(i) zMax_deg(i)], 'r-.');
    hTight = plot(t_s, rad2deg(z_max(i) - lgTube.delta(i, :)), 'r--', t_s, rad2deg(z_min(i) + lgTube.delta(i, :)), 'r--');
    isViolated = usageNo_nd(i, :) > 1;
    hViol = plot(t_s(isViolated), rad2deg(lgNo.z(i, isViolated)), 'r.', 'MarkerSize', 8);
    ylabel(ylab{i});
    if i == 1
        legend([hNo hTube hLim(1) hTight(1) hViol], 'No tube', 'Tube', 'limit', 'tightened limit (tube)', ...
            'no tube outside limit', 'Location', 'eastoutside');
    end
end
xlabel('Time, sec');
sgtitle('Constrained quantities z = [u_{bl}; y_{lim}]: no tube vs tube');

% Zoom on the worst no-tube violation: z = z_ref + c_i e
zoomWindow_s = tWorst_s + [-0.5 1.5];
isZoom = t_s >= zoomWindow_s(1) & t_s <= zoomWindow_s(2);
tz_s = t_s(isZoom);
figure('Name', 'Zoom on the worst violation');
subplot(2, 1, 1); hold on; grid on;
plot(tz_s, rad2deg(lgNo.z(iWorst, isZoom)), 'Color', grayColor, 'LineWidth', 1.5);
plot(tz_s, rad2deg(lgNo.z_ref(iWorst, isZoom)), '--', 'Color', grayColor, 'LineWidth', 1.2);
plot(tz_s, rad2deg(lgTube.z(iWorst, isZoom)), 'b', 'LineWidth', 1.5);
plot(tz_s, rad2deg(lgTube.z_ref(iWorst, isZoom)), 'b--', 'LineWidth', 1.2);
yline([-zMax_deg(iWorst) zMax_deg(iWorst)], 'r-.');
plot(tz_s, rad2deg(z_max(iWorst) - lgTube.delta(iWorst, isZoom)), 'r--', tz_s, rad2deg(z_min(iWorst) + lgTube.delta(iWorst, isZoom)), 'r--');
ylabel(ylab{iWorst});
legend('z, no tube', 'z_{ref}, no tube', 'z, tube', 'z_{ref}, tube', 'limit', '', 'tightened limit', 'Location', 'eastoutside');
title('No tube: z_{ref} rides the limit, so z = z_{ref} + c_i e goes past it. Tube: z_{ref} stays \delta_i inside.');
subplot(2, 1, 2); hold on; grid on;
plot(tz_s, rad2deg(abs(lgNo.z(iWorst, isZoom) - lgNo.z_ref(iWorst, isZoom))), 'Color', grayColor, 'LineWidth', 1.5);
plot(tz_s, rad2deg(abs(lgTube.z(iWorst, isZoom) - lgTube.z_ref(iWorst, isZoom))), 'b', 'LineWidth', 1.5);
plot(tz_s, rad2deg(lgTube.delta(iWorst, isZoom)), 'r--', 'LineWidth', 1.2);
ylabel(['|c_i e|, ' ylab{iWorst}(end-2:end)]);
legend('|c_i e|, no tube', '|c_i e|, tube', '\delta_i (tube)', 'Location', 'eastoutside');
xlabel('Time, sec');
subplot(2, 1, 1); xlim(zoomWindow_s);
subplot(2, 1, 2); xlim(zoomWindow_s);

% Tube bound |c_i e| <= delta_i and V <= V_bar (as in the "Safety margin" figure of Robust_Adaptive_TK.m)
figure('Name', 'Tube bound');
for i = 1:2*m
    subplot(2*m + 1, 1, i); hold on; grid on;
    plot(t_s, rad2deg(abs(lgNo.z(i, :) - lgNo.z_ref(i, :))), 'Color', grayColor, 'LineWidth', 1.0);
    plot(t_s, rad2deg(abs(lgTube.z(i, :) - lgTube.z_ref(i, :))), 'b', t_s, rad2deg(lgTube.delta(i, :)), 'r--', 'LineWidth', 1.2);
    ylabel(ylab{i});
end
subplot(2*m + 1, 1, 1); legend('|c_i e|, no tube', '|c_i e|, tube', '\delta_i', 'Location', 'eastoutside');
title('Tube bound |c_i e| \leq \delta_i (Lemmas 5 and 9)');
subplot(2*m + 1, 1, 2*m + 1);
semilogy(t_s, lgTube.V, 'b', t_s, lgTube.V_bar, 'r--', 'LineWidth', 1.2); grid on;
ylabel('V'); legend('V', 'V_{bar}', 'Location', 'eastoutside'); xlabel('Time, sec');

% Command tracking (as in Robust_Adaptive_TK.m)
figure('Name', 'Command tracking');
subplot(2, 1, 1);
plot(t_s, rad2deg(lgTube.y_cmd(1, :)), 'k--', t_s, rad2deg(lgNo.y_reg(1, :)), 'Color', grayColor, 'LineWidth', 1.0); hold on; grid on;
plot(t_s, rad2deg(lgTube.y_reg(1, :)), 'b', 'LineWidth', 1.2);
ylabel('ps, dps'); legend('Cmd', caseNames{1}, caseNames{2}, 'Location', 'best');
subplot(2, 1, 2);
plot(t_s, lgTube.y_cmd(2, :), 'k--', t_s, lgNo.y_reg(2, :), 'Color', grayColor, 'LineWidth', 1.0); hold on; grid on;
plot(t_s, lgTube.y_reg(2, :), 'b', 'LineWidth', 1.2);
ylabel('N_y, g'); xlabel('Time, sec');
sgtitle('Command tracking (N_y is limited by the \beta and rudder constraints)');

%% Local functions
function [chiDot, d] = closedLoopRhs(t_s, chi, par)
% closedLoopRhs: Closed-loop right-hand side (plant, reference model, adaptive law, V_bar) with the barrier on the
% reference model, with or without the tube
%
%=================================
%                                  INPUTS
%=================================
% t_s: time
% chi: closed-loop state [x; x_ref; vec(Omega_hat); V_bar]
% par: struct of model, controller, and limit parameters; par.isTube selects the tube (true) or delta = 0 (false)
%=================================
% OUTPUTS
%=================================
% chiDot: time derivative of chi
% d: diagnostics (z, z_ref, delta, F margin, u)
%
%=================================
% EXAMPLE USAGE
%=================================
%
% chiDot = closedLoopRhs(0, chi, par)
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
    Omega_hat_dot = par.B' * par.P * e * psi' * par.Gamma;        % (49)
    eQe = e' * par.Q * e;
    V_bar_dot = -eQe;
    pi_0 = [-par.y_cmd(t_s); zeros(m, 1)];
    K_hat = Omega_hat(:, 1:n);
    Theta_hat = Omega_hat(:, n+1:n+q)';                           % q x m

    zeta = [-Omega_hat * psi_ref; par.C_lim * x_p_ref];           % (52)
    G = [-K_hat - Theta_hat' * par.J_Phi(x_p_ref) * par.E_p; par.C_lim * par.E_p];   % (53)
    dzeta_dt = [-Omega_hat_dot * psi_ref; zeros(m, 1)];           % (54)
    H_pi = G * par.B_pi;                                          % (55)
    ell = dzeta_dt + G * par.A_ref * x_ref + H_pi * pi_0 + par.alpha * zeta;   % (56)

    delta = zeros(2*m, 1);
    deltaDot = zeros(2*m, 1);
    if par.isTube
        c_y = [zeros(m) par.C_lim];
        s_y = sum((c_y * par.P_inv) .* c_y, 2);
        delta_y = sqrt(s_y * V_bar);                              % (29)
        deltaDot_y = -s_y * eQe ./ (2 * delta_y);                 % (30) with c_i constant

        K_hat_dot = Omega_hat_dot(:, 1:n);
        Theta_hat_dot = Omega_hat_dot(:, n+1:n+q)';
        s_u = sum((K_hat * par.P_inv) .* K_hat, 2);
        sDot_u = 2 * sum((K_hat_dot * par.P_inv) .* K_hat, 2);
        nu = sqrt(sum(Theta_hat.^2, 1)' + par.epsilon_nu^2);
        nuDot = sum(Theta_hat .* Theta_hat_dot, 1)' ./ nu;
        delta_u = sqrt(s_u * V_bar) + par.L_Phi * nu * sqrt(par.p_bar * V_bar);             % (57)
        deltaDot_u = (sDot_u * V_bar + s_u * V_bar_dot) ./ (2 * sqrt(s_u * V_bar)) ...
            + par.L_Phi * sqrt(par.p_bar) * (nuDot * sqrt(V_bar) + nu * V_bar_dot / (2 * sqrt(V_bar)));   % (58)
        delta = [delta_u; delta_y];
        deltaDot = [deltaDot_u; deltaDot_y];
    end
    DeltaH_max = ell + deltaDot - par.alpha * (par.z_max - delta);   % (37)
    DeltaH_min = -ell + deltaDot + par.alpha * (par.z_min + delta);  % (38)
    pi_aug = H_pi \ (max(0, DeltaH_min) - max(0, DeltaH_max));      % (40)
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
    end
end
