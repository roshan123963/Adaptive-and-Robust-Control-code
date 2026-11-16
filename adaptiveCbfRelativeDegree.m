% adaptiveCbfRelativeDegree: Tube vs no tube for the adaptive PI servo-controller with barrier augmentation on the
% reference model and nonlinear matched uncertainty Theta' * Phi(x_p), with limited outputs of vector relative degree
% r = (r_1, ..., r_m).
%
% Theory: TCST_PI_note/TCST_PI_Sections.pdf, Sections II-V (vector relative degree) and VI (nonlinear uncertainty).
% Each y_lim row of relative degree r_j is constrained through q_j(d/dt) applied to the tightened constraint function,
% with q_j(s) = (s + alpha) * prod_k (s + lambda_k) of degree r_j, and its tube width is passed through r_j - 1
% first-order filters of rate mu, so that its derivatives are computable. The u_bl rows have relative degree 1.
%
% Model: lateral-directional, y_lim = [p_s; beta], with the rudder side force Y_delr / V = 0, so r = (1, 2): beta is
% moved by the rudder through the yaw rate, and its decoupling row is C_beta * A_ref,pp * B_p.
% Both cases use the same plant, adaptive law, and barrier; the only difference is the tube delta_i:
%   Tube:    the reference model is kept in the box tightened by delta_i (Theorem 2)
%   No tube: delta_i = 0, so the reference model is kept in the original box, and the plant, z = z_ref + c_i * e,
%            can leave it during adaptation
% The total command u_bl + w is also constrained (Section VII): w solves a QP that enforces the output-row barrier
% conditions (Y) and the tightened total-command box (U) of the reference model, and v follows in closed form from the
% u_bl rows. Theorem 5 then bounds the total command of the plant, provided (Y) and (U) stay compatible (Assumption 8),
% which the QP solver checks at every step. If the total command leaves its box, the script prints the terms of Lemma 8.
%
%=================================
%                                  INPUTS
%=================================
% None (script). Edit the "Configuration" section to change the relative degrees, rates, uncertainty, or limits.
%=================================
% OUTPUTS
%=================================
% Command-window summary (assumption checks, peaks of z and of u_bl + w, no-tube violations, theorem and proposition
% checks for the tube case, Lemma 8 terms) and figures: one per case, no tube vs tube, total command, tube bound, the beta
% cascade, and command tracking.
%
%=================================
% EXAMPLE USAGE
%=================================
%
% adaptiveCbfRelativeDegree

clc; clear; close all;

%% Configuration
dt_s = 1e-4;                          % RK4 step
finalTime_s = 40;
logEvery_nd = 10;                     % log every logEvery_nd steps

% Design model (Sections II-III): x_p = [beta; p_s; r_s], u = [delta_ail; delta_rud], Y_delr / V = 0
A_p = [-0.11794  0.00085  -1.0001;
       -7.0113  -1.4492    0.22059;
        6.3035   0.06511  -0.41172];
B_p = [ 0        0;
       -7.9662   2.6875;
        0.60926 -2.3577];
C_reg = [ 0       1         0;
         -2.6049  0.018724  0.067695];           % y_reg = [p_s; N_y]
D_reg = [0  0;
         0  0.33698];
C_lim = [0 1 0;
         1 0 0];                                 % y_lim = [p_s; beta]
n_p = 3;  m = 2;  n = n_p + m;  q = 2;
relDeg = [1; 2];                                 % r_j of y_lim = [p_s; beta]

% Nonlinear matched uncertainty (Section VI)
betaScale_rad = deg2rad(0.5);
yawRateScale_radps = deg2rad(5);
Phi = @(x_p) [betaScale_rad * tanh(x_p(1) / betaScale_rad); yawRateScale_radps * tanh(x_p(3) / yawRateScale_radps)];
J_Phi = @(x_p) [sech(x_p(1) / betaScale_rad)^2 0 0; 0 0 sech(x_p(3) / yawRateScale_radps)^2];
L_Phi = 1;
Theta = [0.4  0.5;
         0.2 -0.4];
theta_max = 1.4;

% Adaptation: Gamma_I = 0 keeps K_hat_I frozen (Remark 5.1)
gamma_nd = 2e5;                                  % the tube scales as theta_max / sqrt(gamma_nd); a narrower tube widens (U)
Gamma = diag([zeros(1, m), gamma_nd * ones(1, n_p + q)]);
V_bar_0 = theta_max^2 / gamma_nd;

% Operational limits and barrier rates. alpha holds the last-stage rate of every row; lambdaRates{j} holds the r_j - 1
% cascade rates of y_lim row j
u_max = deg2rad([3; 2]);                         % aileron, rudder
y_max = [deg2rad(18); deg2rad(0.5)];             % p_s, beta
z_max = [u_max; y_max];
z_min = -z_max;
alpha = [200; 200; 40; 3];                       % slow beta rates widen (Y), so (Y) and (U) stay compatible (Assumption 8)
lambdaRates = {[], 1};
mu = 4;                                          % tube filter rate; b_{i,l} >= 0 makes (F) static (Lemma 6)
epsilon_nu = 1e-3;

% Command: p_s and N_y pulses
rollRateCmd_radps = deg2rad(15);
lateralAccelCmd_nd = 0.02;
y_cmd = @(t_s) [rollRateCmd_radps; lateralAccelCmd_nd] * ((t_s >= 0 && t_s < 10) - (t_s >= 20 && t_s < 30));

caseNames = {'No tube', 'Tube'};
isTubeCase = [false true];
caseColor = {[0.6 0.6 0.6], [0 0 1]};           % no tube gray, tube blue

%% Design and offline checks of the assumptions
A_nom = [zeros(m) C_reg; zeros(n_p, m) A_p];
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
Omega_star = [K_star Theta'];
Omega_hat_0 = [K_lqr zeros(m, q)];

% Assumption 4: C_j * A_ref,pp^k * B_p = 0 for k <= r_j - 2, and D_y nonsingular
A_refpp = A_ref(m+1:n, m+1:n);
D_y = zeros(m);
for j = 1:m
    for k = 0:relDeg(j) - 2
        if norm(C_lim(j, :) * A_refpp^k * B_p) > 0
            error('Row %d does not have relative degree %d: C_lim(%d, :) * A_ref,pp^%d * B_p ~= 0.', j, relDeg(j), j, k);
        end
    end
    D_y(j, :) = C_lim(j, :) * A_refpp^(relDeg(j) - 1) * B_p;
end
fprintf('Assumption 2: ||Omega* - Omega_hat(0)||_F = %.3f < theta_max = %.3f\n', norm(Omega_star - Omega_hat_0, 'fro'), theta_max);
fprintf('Assumption 4: r = %s, D_y = %s, cond %.1f, kappa = sigma_min(K_I) = %.3f\n', mat2str(relDeg'), mat2str(D_y, 4), ...
    cond(D_y), min(svd(K_lqr(:, 1:m))));
if sum(relDeg) < n_p
    fprintf('Assumption 6: transmission zeros %s\n', mat2str(tzero(ss(A_p, B_p, C_lim, 0)).', 4));
else
    fprintf('Assumption 6: sum of relative degrees = n_p = %d, no zero dynamics\n', n_p);
end
fprintf('Simulated (true) plant: eig(A_p) = %s\n', mat2str(eig(A(m+1:n, m+1:n)).', 3));

% Rows of the barrier: q_i, b_{i,l}, cascade polynomials p_{i,k}, H_pi row, and c_i * q_i(A_ref)
rowData = struct('rho', {}, 'c', {}, 'q0', {}, 'b', {}, 'cqA', {}, 'H', {}, 'p0', {}, 'cpA', {}, 'bp', {}, 'filtIdx', {}, 's', {});
iFilt = 0;
for j = 1:m
    i = m + j;
    rho = relDeg(j);
    lam = lambdaRates{j};
    if numel(lam) ~= rho - 1
        error('Row %d needs %d cascade rates.', j, rho - 1);
    end
    c_i = [zeros(1, m) C_lim(j, :)];
    qPoly = poly(-[alpha(i) lam(:)']);           % descending coefficients of q_i
    row.rho = rho;
    row.c = c_i;
    row.q0 = qPoly(end);
    row.b = shiftCoefficients(qPoly, mu);
    row.cqA = c_i * polyvalm(qPoly, A_ref);
    row.H = c_i * A_ref^(rho - 1) * B_pi;
    row.p0 = zeros(1, rho - 1);
    row.cpA = zeros(rho - 1, n);
    row.bp = cell(1, rho - 1);
    for k = 1:rho - 1
        pPoly = poly(-lam(1:k));
        row.p0(k) = pPoly(end);
        row.cpA(k, :) = c_i * polyvalm(pPoly, A_ref);
        row.bp{k} = shiftCoefficients(pPoly, mu);
    end
    row.filtIdx = iFilt + (1:rho - 1);
    iFilt = iFilt + rho - 1;
    row.s = c_i * P_inv * c_i';
    rowData(j) = row;
    if rho >= 2
        fprintf('Lemma 6, y_lim row %d (r = %d): b = %s (>= 0 needed), 2 * delta(0) / (y_max - y_min) = %.3f (< 1 needed)\n', ...
            j, rho, mat2str(row.b, 4), 2 * sqrt(row.s * V_bar_0) / (2 * y_max(j)));
    end
end

%% Simulation
par = struct('A', A, 'B', B, 'B_ref', B_ref, 'B_v', B_v, 'B_pi', B_pi, 'A_ref', A_ref, 'P', P, 'P_inv', P_inv, 'Q', Q, ...
    'Gamma', Gamma, 'C_lim', C_lim, 'E_p', E_p, 'alpha', alpha, 'z_min', z_min, 'z_max', z_max, 'Theta', Theta, ...
    'Phi', Phi, 'J_Phi', J_Phi, 'L_Phi', L_Phi, 'p_bar', max(eig(E_p * P_inv * E_p')), 'epsilon_nu', epsilon_nu, ...
    'y_cmd', y_cmd, 'm', m, 'n', n, 'q', q, 'mu', mu, 'rows', rowData, 'nFilt', iFilt);
isAdapted = diag(Gamma)' > 0;
nSteps = round(finalTime_s / dt_s);
nLog = floor(nSteps / logEvery_nd) + 1;
results = cell(1, 2);
for c = 1:2
    par.isTube = isTubeCase(c);
    filt0 = zeros(par.nFilt, 1);                 % delta_{i,k}(0) = delta_{i,0}(0)
    for j = 1:m
        filt0(rowData(j).filtIdx) = sqrt(rowData(j).s * V_bar_0);
    end
    chi = [zeros(n, 1); zeros(n, 1); Omega_hat_0(:); V_bar_0; filt0];   % x(0) = x_ref(0) = 0 (Assumption 3)
    logData = struct('t_s', zeros(1, nLog), 'z', zeros(2*m, nLog), 'z_ref', zeros(2*m, nLog), 'delta', zeros(2*m, nLog), ...
        'u', zeros(m, nLog), 'w', zeros(m, nLog), 'x', zeros(n, nLog), 'x_ref', zeros(n, nLog), 'V', zeros(1, nLog), ...
        'V_bar', zeros(1, nLog), 'F', zeros(2*m, nLog), 'psiMax', zeros(2*m, nLog), 'stageOrder', zeros(1, nLog), ...
        'y_reg', zeros(m, nLog), 'y_cmd', zeros(m, nLog), 'OmegaErr', zeros(1, nLog), 'K_hat', zeros(m, n, nLog));
    iLog = 0;
    logData.nInfeasible = 0;                     % RK4 steps at which the total-command QP is infeasible
    logData.nTotalActive = 0;                    % RK4 steps at which the total-command box changes w
    for k = 0:nSteps
        t_s = k * dt_s;
        [k1, d] = closedLoopRhs(t_s, chi, par);
        logData.nInfeasible = logData.nInfeasible + ~d.isQpFeasible;
        logData.nTotalActive = logData.nTotalActive + d.isTotalActive;
        if mod(k, logEvery_nd) == 0
            iLog = iLog + 1;
            e = chi(1:n) - chi(n+1:2*n);
            Omega_hat = reshape(chi(2*n+1:2*n+m*(n+q)), m, n + q);
            Delta_Omega = Omega_hat - Omega_star;
            logData.t_s(iLog) = t_s;
            logData.z(:, iLog) = d.z;
            logData.z_ref(:, iLog) = d.z_ref;
            logData.delta(:, iLog) = d.delta;
            logData.u(:, iLog) = d.u;
            logData.w(:, iLog) = d.w;
            logData.x(:, iLog) = chi(1:n);
            logData.x_ref(:, iLog) = chi(n+1:2*n);
            logData.V(iLog) = e' * P * e + sum(sum(Delta_Omega(:, isAdapted).^2 ./ diag(Gamma(isAdapted, isAdapted))'));
            logData.V_bar(iLog) = chi(2*n+m*(n+q)+1);
            logData.F(:, iLog) = d.F;
            logData.psiMax(:, iLog) = d.psiMax;
            logData.stageOrder(iLog) = d.stageOrder;
            logData.y_reg(:, iLog) = C_reg * chi(m+1:n) + D_reg * (d.u + Theta' * Phi(chi(m+1:n)));
            logData.y_cmd(:, iLog) = y_cmd(t_s);
            logData.OmegaErr(iLog) = norm(Delta_Omega, 'fro');
            logData.K_hat(:, :, iLog) = Omega_hat(:, 1:n);
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
dtLog_s = dt_s * logEvery_nd;

%% Summary
zMax_deg = rad2deg(z_max);
rowNames = {'dela', 'delr', 'ps', 'beta'};
fprintf('\n%-30s %8s %8s %8s %8s %10s %10s\n', 'peak / limit', 'dela', 'delr', 'ps', 'beta', 'u_ail', 'u_rud');
for c = 1:2
    lg = results{c};
    fprintf('%-30s %8.3f %8.3f %8.3f %8.3f %10.3f %10.3f\n', caseNames{c}, max(abs(lg.z), [], 2) ./ z_max, ...
        max(abs(lg.u), [], 2) ./ u_max);
end
fprintf('(dela, delr: the baseline u_bl; u_ail, u_rud: the total command u_bl + w; both are constrained)\n');
for c = 1:2
    fprintf('%s: total-command box active for %.3f s, QP infeasible for %.4f s\n', caseNames{c}, ...
        results{c}.nTotalActive * dt_s, results{c}.nInfeasible * dt_s);
end

usageNo_nd = abs(lgNo.z) ./ z_max;
[peakNo_nd, iPeak] = max(usageNo_nd, [], 2);
[worstUsage_nd, iWorst] = max(peakNo_nd);
tWorst_s = t_s(iPeak(iWorst));
fprintf('\nNo tube: worst violation %.1f%% of the %s limit at t = %.2f s\n', 100 * (worstUsage_nd - 1), rowNames{iWorst}, tWorst_s);
fprintf('  at that time ||Omega_hat - Omega*||_F = %.3f (initial %.3f): the adaptation has not converged\n', ...
    lgNo.OmegaErr(iPeak(iWorst)), lgNo.OmegaErr(1));
fprintf('  z_ref,i = %.3f (limit %.3f), c_i e = %+.3f (deg or dps)\n', rad2deg(lgNo.z_ref(iWorst, iPeak(iWorst))), ...
    zMax_deg(iWorst), rad2deg(lgNo.z(iWorst, iPeak(iWorst)) - lgNo.z_ref(iWorst, iPeak(iWorst))));
for i = 1:2*m
    isViolated = usageNo_nd(i, :) > 1;
    if any(isViolated)
        fprintf('  %-4s outside its limit for %.3f s in total\n', rowNames{i}, sum(isViolated) * dtLog_s);
    end
end

checkTol_nd = 1e-9;                              % numerical tolerance on the <= 0 checks
uRefTot = lgTube.z_ref(1:m, :) + lgTube.w;       % u_ref^tot = u_ref + w (computable)
refDist = max(lgTube.z_ref - (z_max - lgTube.delta), (z_min + lgTube.delta) - lgTube.z_ref);
checkNames = {'Theorem 2.1, x_ref in the tightened box', 'Theorem 2.1, every cascade stage psi <= 0', ...
    'Lemmas 5, 6, 9, |c_i e| - delta_i', 'Lemma 5, V - V_bar', 'Lemma 6, filter stages ordered', ...
    'Theorem 2.2, |z_i| - z_max,i (plant)', '(F), minus the margin', 'Proposition 1, |u_i - u_ref,i^tot| - delta_i', ...
    'Assumption 8, QP infeasible (s)', 'Theorem 5, |u_i| - u_max,i (plant total command)'};
checkValues = [max(refDist(:)), max(lgTube.psiMax(:)), max(max(abs(lgTube.z - lgTube.z_ref) - lgTube.delta)), ...
    max(lgTube.V - lgTube.V_bar), max(lgTube.stageOrder), max(max(abs(lgTube.z) - z_max)), -min(lgTube.F(:)), ...
    max(max(abs(lgTube.u - uRefTot) - lgTube.delta(1:m, :))), lgTube.nInfeasible * dt_s, max(max(abs(lgTube.u) - u_max))];
fprintf('\nTube case theorem checks (each <= 0):\n');
for iCheck = 1:numel(checkNames)
    checkStatus = 'FAIL';
    if checkValues(iCheck) <= checkTol_nd
        checkStatus = 'ok';
    end
    fprintf('  %-46s %-4s  %.2e\n', checkNames{iCheck}, checkStatus, checkValues(iCheck));
end

% Total command of the tube case: where it leaves the actuator box, and the terms of Lemma 8 on the beta row
isOutTube = any(abs(lgTube.u) > u_max, 1);
fprintf('\nTube case: total command u_bl + w outside its limits for %.3f s (Theorem 5 requires 0)\n', ...
    sum(isOutTube) * dtLog_s);
if any(isOutTube)
    [~, kWorst] = max(max(abs(lgTube.u) ./ u_max, [], 1));
    jBeta = 2;
    rowBeta = rowData(jBeta);
    x_refW = lgTube.x_ref(:, kWorst);
    yDdot_ref = rowBeta.c * A_ref^2 * x_refW + D_y(jBeta, :) * lgTube.w(:, kWorst);    % Lemma 2, (22)
    drift = C_lim(jBeta, :) * A(m+1:n, m+1:n)^2 * x_refW(m+1:n);                         % C_j A_p^r x_p,ref (true A_p)
    gainErr = D_y(jBeta, :) * (-lgTube.z_ref(1:m, kWorst) - K_star * x_refW);         % d_j (Omega_hat psi_ref - K* x_ref)
    fprintf('  largest excursion %.2f x limit at t = %.3f s: u_bl = %s deg, w = %s deg\n', ...
        max(abs(lgTube.u(:, kWorst)) ./ u_max), t_s(kWorst), mat2str(rad2deg(lgTube.z(1:m, kWorst))', 3), ...
        mat2str(rad2deg(lgTube.w(:, kWorst))', 3));
    fprintf('  Lemma 8, beta row (rad/s^2): d_j u_ref^tot = %+.4f = beta_ref'''' (%+.4f) - drift (%+.4f) - gain error (%+.4f)\n', ...
        D_y(jBeta, :) * uRefTot(:, kWorst), yDdot_ref, drift, gainErr);
    fprintf('  (gain error includes Theta_hat'' * Phi(x_p,ref) of u_ref; identity residual %.1e)\n', ...
        D_y(jBeta, :) * uRefTot(:, kWorst) - (yDdot_ref - drift - gainErr));
end

%% Figures
ylab = {'dela, deg', 'delr, deg', 'ps, dps', 'beta, deg'};
grayColor = caseColor{1};

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
    axis([-2 2 -2 2] .* zMax_deg([1 1 2 2])');
    legend('u_{bl} + w', 'u_{bl} (constrained)', 'Location', 'best');

    subplot(2, 2, 3);
    plot(ps_dps, beta_deg, 'LineWidth', 0.5); hold on; grid on;
    rectangle('Position', [-zMax_deg(3) -zMax_deg(4) 2 * zMax_deg(3) 2 * zMax_deg(4)], 'EdgeColor', 'r', 'LineWidth', 1);
    xlabel('ps (dps)'); ylabel('beta (deg)'); title('Limited Output Operational Bounds');

    rowsToPlot = [4 2 3];                        % beta, delr, ps
    for jj = 1:3
        i = rowsToPlot(jj);
        subplot(4, 2, 2 * jj); hold on; grid on;
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
        if jj == 1
            legend('z (plant)', 'z_{ref}', 'limit', 'Location', 'eastoutside');
        end
    end

    subplot(4, 2, 8);
    plot(t_s, vecnorm(lg.x - lg.x_ref), 'b', 'LineWidth', 1.2); grid on; ylabel('||e||');
    yyaxis right; plot(t_s, lg.OmegaErr, 'LineWidth', 1.0); ylabel('||\Omega_{hat} - \Omega^*||_F');
    xlabel('Time, sec');
    sgtitle(sprintf('MRAC + CBF, r = (1, 2), %s (red dots: outside the limit)', lower(caseNames{c})));
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
    hViol = plot([NaN t_s(isViolated)], [NaN rad2deg(lgNo.z(i, isViolated))], 'r.', 'MarkerSize', 8);   % NaN keeps a handle
    ylabel(ylab{i});
    if i == 1
        legend([hNo hTube hLim(1) hTight(1) hViol], 'No tube', 'Tube', 'limit', 'tightened limit (tube)', ...
            'no tube outside limit', 'Location', 'eastoutside');
    end
end
xlabel('Time, sec');
sgtitle('Constrained quantities z = [u_{bl}; y_{lim}], r = (1, 2): no tube vs tube');

% Total command u = u_bl + w against the actuator limits
figure('Name', 'Total command');
uLab = {'u_{ail} = u_{bl,1} + w_1, deg', 'u_{rud} = u_{bl,2} + w_2, deg'};
for i = 1:m
    subplot(m, 1, i); hold on; grid on;
    hCase = gobjects(1, 2);
    for c = 1:2
        hCase(c) = plot(t_s, rad2deg(results{c}.u(i, :)), 'Color', caseColor{c}, 'LineWidth', 1.2);
    end
    hLim = yline([-zMax_deg(i) zMax_deg(i)], 'r-.');
    isOut = abs(lgTube.u(i, :)) > u_max(i);
    hOut = plot([NaN t_s(isOut)], [NaN rad2deg(lgTube.u(i, isOut))], 'r.', 'MarkerSize', 8);   % NaN keeps a handle if empty
    ylabel(uLab{i});
    if i == 1
        legend([hCase hLim(1) hOut], caseNames{1}, caseNames{2}, 'limit', 'tube case outside limit', 'Location', 'eastoutside');
    end
end
xlabel('Time, sec');
sgtitle('Total command u = u_{bl} + w, constrained through (U) (Theorem 5)');

% Tube bound |c_i e| <= delta_i and V <= V_bar (as in the "Safety margin" figure of Robust_Adaptive_TK.m)
figure('Name', 'Tube bound');
for i = 1:2*m
    subplot(2*m + 1, 1, i); hold on; grid on;
    plot(t_s, rad2deg(abs(lgNo.z(i, :) - lgNo.z_ref(i, :))), 'Color', grayColor, 'LineWidth', 1.0);
    plot(t_s, rad2deg(abs(lgTube.z(i, :) - lgTube.z_ref(i, :))), 'b', t_s, rad2deg(lgTube.delta(i, :)), 'r--', 'LineWidth', 1.2);
    ylabel(ylab{i});
end
subplot(2*m + 1, 1, 1); legend('|c_i e|, no tube', '|c_i e|, tube', '\delta_i', 'Location', 'eastoutside');
title('Tube bound |c_i e| \leq \delta_i (Lemmas 5, 6, and 9)');
subplot(2*m + 1, 1, 2*m + 1);
semilogy(t_s, lgTube.V, 'b', t_s, lgTube.V_bar, 'r--', 'LineWidth', 1.2); grid on;
ylabel('V'); legend('V', 'V_{bar}', 'Location', 'eastoutside'); xlabel('Time, sec');

% The beta cascade of the tube case: h <= 0 and every stage psi_k <= 0
iBeta = m + 2;
hMax = lgTube.z_ref(iBeta, :) - z_max(iBeta) + lgTube.delta(iBeta, :);
hMin = z_min(iBeta) + lgTube.delta(iBeta, :) - lgTube.z_ref(iBeta, :);
figure('Name', 'Beta cascade');
subplot(2, 1, 1); hold on; grid on;
plot(t_s, rad2deg(max(hMax, hMin)), 'b', 'LineWidth', 1.2); yline(0, 'r-.');
ylabel('max(h_{max}, h_{min}), deg'); title('Tube case: reference model in the tightened \beta box');
subplot(2, 1, 2); hold on; grid on;
plot(t_s, lgTube.psiMax(iBeta, :), 'b', 'LineWidth', 1.2); yline(0, 'r-.');
ylabel('max_k \psi_k'); title('Every cascade stage \psi_k = p_k(d/dt) h \leq 0'); xlabel('Time, sec');

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
function b = shiftCoefficients(pPoly, mu)
% shiftCoefficients: Coefficients b_0, ..., b_{d-1} (ascending) of (p(s) - p(0)) / s written in powers of (s + mu)
%
%=================================
%                                  INPUTS
%=================================
% pPoly: descending coefficients of a monic polynomial p of degree d >= 1
% mu: shift (filter rate)
%=================================
% OUTPUTS
%=================================
% b: 1 x d row, (p(s) - p(0)) / s = sum_l b(l + 1) * (s + mu)^l
%
%=================================
% EXAMPLE USAGE
%=================================
%
% b = shiftCoefficients(poly([-2 -5]), 7)       % returns [0 1]
    a = fliplr(pPoly(1:end-1));                  % ascending coefficients of (p(s) - p(0)) / s
    d = numel(a);
    b = zeros(1, d);
    for l = 0:d-1                                % s^l = (sigma - mu)^l with sigma = s + mu
        for k = 0:l
            b(k + 1) = b(k + 1) + a(l + 1) * nchoosek(l, k) * (-mu)^(l - k);
        end
    end
end

function [chiDot, d] = closedLoopRhs(t_s, chi, par)
% closedLoopRhs: Closed-loop right-hand side (plant, reference model, adaptive law, V_bar, tube filters) with the barrier
% on the reference model, with or without the tube; y_lim row j has relative degree par.rows(j).rho >= 1
%
%=================================
%                                  INPUTS
%=================================
% t_s: time
% chi: closed-loop state [x; x_ref; vec(Omega_hat); V_bar; tube filter states delta_{i,k}, k >= 1]
% par: struct of model, controller, limit, and barrier-row parameters; par.isTube selects the tube or delta = 0
%=================================
% OUTPUTS
%=================================
% chiDot: time derivative of chi
% d: diagnostics (z, z_ref, delta, F margin, largest cascade stage, filter ordering, u, w)
%
%=================================
% EXAMPLE USAGE
%=================================
%
% chiDot = closedLoopRhs(0, chi, par)
    m = par.m;  n = par.n;  q = par.q;
    x = chi(1:n);
    x_ref = chi(n+1:2*n);
    iOmega = 2*n + (1:m*(n+q));
    Omega_hat = reshape(chi(iOmega), m, n + q);
    V_bar = chi(iOmega(end) + 1);
    filt = chi(iOmega(end) + 2:end);
    x_p = x(m+1:n);
    x_p_ref = x_ref(m+1:n);

    psiReg = [x; par.Phi(x_p)];
    psiReg_ref = [x_ref; par.Phi(x_p_ref)];
    e = x - x_ref;
    u_bl = -Omega_hat * psiReg;                                   % (48)
    Omega_hat_dot = par.B' * par.P * e * psiReg' * par.Gamma;      % (49)
    eQe = e' * par.Q * e;
    V_bar_dot = -eQe;
    pi_0 = [-par.y_cmd(t_s); zeros(m, 1)];
    K_hat = Omega_hat(:, 1:n);
    Theta_hat = Omega_hat(:, n+1:n+q)';

    % Input rows: relative degree 1, nonlinear tube (57)-(58)
    alpha_u = par.alpha(1:m);
    if par.isTube
        K_hat_dot = Omega_hat_dot(:, 1:n);
        Theta_hat_dot = Omega_hat_dot(:, n+1:n+q)';
        s_u = sum((K_hat * par.P_inv) .* K_hat, 2);
        sDot_u = 2 * sum((K_hat_dot * par.P_inv) .* K_hat, 2);
        nu = sqrt(sum(Theta_hat.^2, 1)' + par.epsilon_nu^2);
        nuDot = sum(Theta_hat .* Theta_hat_dot, 1)' ./ nu;
        delta_u = sqrt(s_u * V_bar) + par.L_Phi * nu * sqrt(par.p_bar * V_bar);
        deltaDot_u = (sDot_u * V_bar + s_u * V_bar_dot) ./ (2 * sqrt(s_u * V_bar)) ...
            + par.L_Phi * sqrt(par.p_bar) * (nuDot * sqrt(V_bar) + nu * V_bar_dot / (2 * sqrt(V_bar)));
    else
        delta_u = zeros(m, 1);
        deltaDot_u = zeros(m, 1);
    end
    G_u = -K_hat - Theta_hat' * par.J_Phi(x_p_ref) * par.E_p;       % (53)
    H_pi = zeros(2*m);
    H_pi(1:m, :) = G_u * par.B_pi;                                 % (55)
    ell_u = -Omega_hat_dot * psiReg_ref + G_u * par.A_ref * x_ref + H_pi(1:m, :) * pi_0 - alpha_u .* (Omega_hat * psiReg_ref);
    qDelta_u = deltaDot_u + alpha_u .* delta_u;                    % (33) with rho = 1
    DeltaH_max = [ell_u + qDelta_u - alpha_u .* par.z_max(1:m); zeros(m, 1)];   % (37)
    DeltaH_min = [-ell_u + qDelta_u + alpha_u .* par.z_min(1:m); zeros(m, 1)];  % (38)
    delta = [delta_u; zeros(m, 1)];
    F = [alpha_u .* (par.z_max(1:m) - par.z_min(1:m)) - 2 * qDelta_u; zeros(m, 1)];
    psiMax = -inf(2*m, 1);
    filtDot = zeros(size(filt));
    stageOrder = -inf;

    % Output rows of relative degree rho: q(d/dt) applied to the tightened constraint, with rho - 1 tube filters (31)
    for j = 1:m
        i = m + j;
        row = par.rows(j);
        rho = row.rho;
        if par.isTube
            delta0 = sqrt(row.s * V_bar);                          % (29)
            stages = [delta0; filt(row.filtIdx)];                  % delta_{i,0}, ..., delta_{i,rho-1}
            epsStage = [-row.s * eQe / (2 * delta0); par.mu * (stages(1:end-1) - stages(2:end))];   % (30), (31)
            filtDot(row.filtIdx) = epsStage(2:end);
            stageOrder = max([stageOrder; stages(1:end-1) - stages(2:end)]);
        else
            stages = zeros(rho, 1);
            epsStage = zeros(rho, 1);
        end
        delta_j = stages(end);
        qDelta = row.q0 * delta_j + sum(row.b .* par.mu.^(0:rho-1) .* epsStage(rho:-1:1)');            % (33)
        H_pi(i, :) = row.H;                                                                              % (23)
        ell_j = row.cqA * x_ref + row.H * pi_0;                                                          % (23)
        DeltaH_max(i) = ell_j + qDelta - row.q0 * par.z_max(i);                                          % (37)
        DeltaH_min(i) = -ell_j + qDelta + row.q0 * par.z_min(i);                                         % (38)
        F(i) = row.q0 * (par.z_max(i) - par.z_min(i)) - 2 * qDelta;                                      % (F)
        delta(i) = delta_j;
        if rho >= 2
            hMax = row.c * x_ref - par.z_max(i) + delta_j;
            hMin = par.z_min(i) + delta_j - row.c * x_ref;
            psiMax(i) = max(hMax, hMin);
            for k = 1:rho - 1                                      % psi_k = p_k(d/dt) h, which does not contain pi (36)
                pDelta = row.p0(k) * delta_j + sum(row.bp{k} .* par.mu.^(0:k-1) .* epsStage(rho:-1:rho-k+1)');
                psiMax(i) = max([psiMax(i), row.cpA(k, :) * x_ref + pDelta - row.p0(k) * par.z_max(i), ...
                    -row.cpA(k, :) * x_ref + pDelta + row.p0(k) * par.z_min(i)]);
            end
        end
    end

    % Augmentation with the total-command constraint (Section VII): w from the QP in nu = D_y * w, then v in closed form
    D_y = H_pi(m+1:2*m, m+1:2*m);
    u_ref = -Omega_hat * psiReg_ref;
    w_lo = par.z_min(1:m) + delta_u - u_ref;                       % (U): u_min + delta <= u_ref + w <= u_max - delta
    w_hi = par.z_max(1:m) - delta_u - u_ref;
    [nu, isQpFeasible, isTotalActive] = solveTotalCommandQp(D_y, DeltaH_min(m+1:2*m), -DeltaH_max(m+1:2*m), w_lo, w_hi);
    w = D_y \ nu;
    s_u = max(0, DeltaH_min(1:m)) - max(0, DeltaH_max(1:m));       % input rows of (H_pi * pi), as in (40)
    v = H_pi(1:m, 1:m) \ (s_u - H_pi(1:m, m+1:2*m) * w);
    pi_aug = [v; w];
    u = u_bl + w;

    xDot = par.A * x + par.B * (u + par.Theta' * par.Phi(x_p)) + par.B_ref * par.y_cmd(t_s) + par.B_v * v;   % (47)
    x_refDot = par.A_ref * x_ref + par.B_pi * (pi_0 + pi_aug);                                        % (17)
    chiDot = [xDot; x_refDot; Omega_hat_dot(:); V_bar_dot; filtDot];

    if nargout > 1
        d.z = [u_bl; par.C_lim * x_p];
        d.z_ref = [-Omega_hat * psiReg_ref; par.C_lim * x_p_ref];
        d.delta = delta;
        d.F = F;
        d.psiMax = psiMax;
        d.stageOrder = stageOrder;
        d.u = u;
        d.w = w;
        d.isQpFeasible = isQpFeasible;
        d.isTotalActive = isTotalActive;
    end
end

function [nu, isFeasible, isTotalActive] = solveTotalCommandQp(D_y, lo_y, hi_y, w_lo, w_hi)
% solveTotalCommandQp: Minimum-norm nu = D_y * w subject to the output-row barrier conditions lo_y <= nu <= hi_y and the
% total-command box w_lo <= D_y^-1 * nu <= w_hi, by enumeration of the faces of the feasible polytope (small m)
%
%=================================
%                                  INPUTS
%=================================
% D_y: m x m decoupling matrix
% lo_y, hi_y: bounds on nu from the output rows, DeltaH_min and -DeltaH_max of the y_lim rows
% w_lo, w_hi: bounds on w from the tightened total-command box
%=================================
% OUTPUTS
%=================================
% nu: minimizer; if the QP is infeasible, the output-row solution is kept and the total-command box is dropped
% isFeasible: true if the QP is feasible
% isTotalActive: true if the total-command box changes the solution
%
%=================================
% EXAMPLE USAGE
%=================================
%
% nu = solveTotalCommandQp(eye(2), [-1; -1], [1; 1], [-2; -2], [2; 2])     % returns [0; 0]
    m = numel(lo_y);
    D_inv = D_y \ eye(m);
    A = [eye(m); -eye(m); D_inv; -D_inv];
    b = [hi_y; -lo_y; w_hi; -w_lo];
    tol = 1e-12 * max(1, max(abs(b)));
    nu = min(max(0, lo_y), hi_y);                % output rows alone: the closed form (40)
    isTotalActive = false;
    isFeasible = true;
    if all(A * nu <= b + tol)
        return
    end
    isTotalActive = true;
    isFeasible = false;
    bestNorm = inf;
    for k = 1:m                                  % the minimizer is the min-norm point of the face that contains it
        activeSets = nchoosek(1:4*m, k);
        for s = 1:size(activeSets, 1)
            A_S = A(activeSets(s, :), :);
            gram = A_S * A_S';
            if rcond(gram) < 1e-12
                continue
            end
            candidate = A_S' * (gram \ b(activeSets(s, :)));
            if all(A * candidate <= b + tol) && norm(candidate) < bestNorm
                bestNorm = norm(candidate);
                nu = candidate;
                isFeasible = true;
            end
        end
    end
end
