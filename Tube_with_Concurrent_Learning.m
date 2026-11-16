%% MRAC + CBF: no tube  vs  tube  vs  tube + CONCURRENT LEARNING
% Case 1: MRAC + CBF, no tube          (original baseline)
% Case 2: MRAC + CBF + tube            (Theorem 2)
% Case 3: MRAC + CBF + tube + CL       (concurrent learning)
%   K_hat_dot = B'*P*e*x'*Gamma - kce * sum_j eps_j*x_j' * Gamma
%   eps_j = K_hat*x_j - K_lqr*x_j - y_j = -K_tilde*x_j,  y_j = Theta*x_j (from data)
%   V_dot = -e'*Q_tilde*e - 2*kce*sum_j ||eps_j||^2   (exact, computable)
clc; clear; close all;

g = 32.174;                          % gravitational acceleration (ft/s^2)

Yb_V0 = -0.11794;
Yps_V_0 = 0.00085;
Yrs_V0_m1 = -1.0001;
L_b = -7.0113;
L_ps = -1.4492;
L_rs = 0.22059;
N_b = 6.3035;
N_ps = 0.06511;
N_rs = -0.41172;

Yail_V0 = 0;
Yrud_V0 = 0.015257;
L_ail = -7.9662;
L_rud = 2.6875;
N_ail = 0.60926;
N_rud = -2.3577;

Yb_g = -2.6049;
Yps_g = 0.018724;
Yrs_g = 0.067695;

Yail_g = 0;
Yrud_g = 0.33698;

A_p = [Yb_V0  Yps_V_0    Yrs_V0_m1;
       L_b      L_ps        L_rs;
       N_b      N_ps        N_rs];

B_p = [Yail_V0    Yrud_V0;
      L_ail         L_rud;
      N_ail          N_rud];

C_p_reg = [0    1   0;
          Yb_g    Yps_g   Yrs_g];

D_p_reg = [0    0;
          Yail_g     Yrud_g];

% Extended Plant (nominal model, known to the controller)
A = [zeros(2, 2) C_p_reg;
    zeros(3, 2) A_p];

B = [D_p_reg;
    B_p];

B_ref = [-eye(2);
    zeros(3,2)];

B_v = [eye(2);
    zeros(3,2)];

% LQR
Q_lqr = diag([1.025 1.0289 0 0 1.6021]);
R_lqr = diag([1 0.49129]);
K_lqr = lqr(A,B,Q_lqr,R_lqr);

A_ref = A - B * K_lqr;

% True Plant (unknown to the controller, matched uncertainty)
Theta  = 0.5 * [zeros(2,2) K_lqr(:,3:5)];
A_true = A + B*Theta;
K_star = K_lqr + Theta;              % ideal gains: A_true - B*K_star = A_ref

% Observer-like reference model
nu  = 0.1;
Q_v = eye(5) + ((nu + 1)/nu) * eye(5);
R_v = (nu/(nu + 1)) * eye(5);
P_v = icare(A_ref', eye(5), Q_v, R_v);
L_v = P_v / R_v;
P   = inv(P_v);

% Adaptive law
gamma   = 1e5;
Gamma_I = gamma*eye(2);
Gamma_P = gamma*eye(3);
Gamma   = blkdiag(Gamma_I, Gamma_P);

% Tube (Theorem 2)
theta_max = 1.5;                        % get from "norm(Theta)", under root of sum of squares of all entries of theta
V_bar0    = theta_max^2 / gamma;
Q_tilde   = -((A_ref - L_v)'*P + P*(A_ref - L_v));
if norm(K_star - K_lqr, 'fro') > theta_max
    warning('theta_max = %.3f does not bound ||Theta||_F = %.3f: tube is NOT certified.', ...
            theta_max, norm(K_star - K_lqr, 'fro'));
end

% Concurrent learning parameters
kc      = 0.5;       % CL gain; normalized so kce*gamma*lambda_max < kc*gamma (Euler: kc*gamma*h < 2)
N_max   = 20;        % history stack size
%nov_tol = 0.1;       % relative novelty threshold for recording a point
nov_tol = 1e-3;      % epsilon in eq. 15 (paper's example uses 0.08 in its own units)
%x_min   = 1e-3;      % ignore near-zero states
B_pinv  = pinv(B);

% Input and Output Constraints
Y_max = deg2rad([2; 1; 18; 0.5]);       % aileron, rudder, roll rate, sideslip angle                % with looser limits (large value of Ymax) the CBF rarely hits them, so the CBF does almost nothing
Y_min = -Y_max;
C_plim = [0 1 0;
          1 0 0];
H_w = C_plim * B_p;
alpha = diag([80 8 40 4]);

% Command Signal
dt = 1e-4; T = 40;
n_sub = 10;
t = 0:dt:T; nT = numel(t);
pulse = zeros(1,nT);
pulse(t >= 0  & t < 10) =  1;
pulse(t >= 20 & t < 30) = -1;
%pulse(t >= 40  & t < 50) =  1;
%pulse(t >= 60 & t < 70) = -1;
Y_cmd = [deg2rad(15)*pulse; 0.0312*pulse];   % ps_cmd Ny_cmd            % what the pilot ask the aircraft to do at each time step

% Cases
case_names = {'MRAC + CBF, no tube', 'MRAC + CBF + tube', 'MRAC + CBF + tube + CL'};
adapt_on   = [1 1 1];
tube_on    = [0 1 1];
cl_on      = [0 0 1];
nC         = numel(case_names);
col        = {[0.65 0.65 0.65], [0.85 0.33 0.10], [0 0.25 0.95]};   % grey, orange, blue
lw         = [1.0 1.2 1.2];

X = zeros(5,nT,nC);  X_ref = zeros(5,nT,nC);  U = zeros(2,nT,nC);  VW = zeros(4,nT,nC);
K_hist = zeros(2,5,nT,nC);
Y_lim = zeros(4,nT,nC);  Y_lim_ref = zeros(4,nT,nC);  D_hist = zeros(4,nT,nC);  Ddot_hist = zeros(4,nT,nC);
F_min = inf(4,nC);  Vbar_hist = zeros(nT,nC);  sig_hist = zeros(nT,nC);

for c = 1:nC
    x = zeros(5,1);  x_ref = zeros(5,1);
    K_I = K_lqr(:,1:2);  K_P = K_lqr(:,3:5);
    V_bar = V_bar0;
    Xs = zeros(5,0);  Ys = zeros(2,0);  kce = 0;  sig_cur = 0;
    %x_pre = x;  u_pre = zeros(2,1);
    
    x_pre = x;  e_pre = zeros(5,1);  K_pre = [K_I K_P];  u_pre = zeros(2,1);

    for k = 1:nT

        e_yi = x(1:2);
        x_p = x(3:5);
        e = x - x_ref;

        % Adaptive laws (MRAC + CL)
        [Kdot_cl, diss] = cl_term([K_I K_P], Xs, Ys, K_lqr, kce, Gamma);
        K_I_dot = adapt_on(c) * B' * P * e * e_yi' * Gamma_I + Kdot_cl(:,1:2);
        K_P_dot = adapt_on(c) * B' * P * e * x_p' * Gamma_P + Kdot_cl(:,3:5);

        % Baseline with (adaptive) gains
        U_bl = - K_I * e_yi - K_P * x_p;

        % CBF matrices
        C_lim = [-K_I   -K_P;
            zeros(2,2)           C_plim];
        H_x = C_lim * A_ref + alpha * C_lim;
        H_u = [-K_I        -K_I*D_p_reg - K_P*B_p;
               zeros(2,2)   H_w];
        E     = [-[K_I_dot K_P_dot]*x_ref; zeros(2,1)] + C_lim*L_v*e;       % E is built from x_ref instead of x (for tube case)

        % Tube and its time derivative (V_bar_dot = -eQe - diss); zero when tube_on = 0
        eQe   = e' * Q_tilde * e;
        s     = diag(C_lim * P_v * C_lim');
        s_dot = 2*diag([-[K_I_dot K_P_dot]; zeros(2,5)] * P_v * C_lim');
        delta     = tube_on(c) * sqrt(s*V_bar);
        delta_dot = tube_on(c) * (s_dot*V_bar - s*(eQe + diss)) ./ (2*sqrt(s*V_bar));

        % CBF on the reference model with tightened limits
        q = [-Y_cmd(:,k); zeros(2,1)];
        H_max_delta =  H_x*x_ref + H_u*q + E + delta_dot - alpha*(Y_max - delta);       % x is replaced with x_ref (for tube case)
        H_min_delta = -H_x*x_ref - H_u*q - E + delta_dot + alpha*(Y_min + delta);
        vw = H_u \ (max(0, H_min_delta) - max(0, H_max_delta));
        V = vw(1:2);  W = vw(3:4);

        u = U_bl + W;

        % Log
        X(:,k,c) = x;  X_ref(:,k,c) = x_ref;  U(:,k,c) = u;  VW(:,k,c) = vw;
        K_hist(:,:,k,c) = [K_I K_P];
        Y_lim(:,k,c) = C_lim*x;  Y_lim_ref(:,k,c) = C_lim*x_ref;
        D_hist(:,k,c) = delta;  Ddot_hist(:,k,c) = delta_dot;
        Vbar_hist(k,c) = V_bar;  sig_hist(k,c) = sig_cur;
        F_min(:,c) = min(F_min(:,c), diag(alpha).*(Y_max - Y_min - 2*delta) - 2*delta_dot);

        % Dynamics: n_sub Euler sub-steps (V, W held over dt)
        h = dt/n_sub;
        for j = 1:n_sub
            e = x - x_ref;
            [Kdot_cl, diss] = cl_term([K_I K_P], Xs, Ys, K_lqr, kce, Gamma);
            K_I_dot = adapt_on(c) * B' * P * e * x(1:2)' * Gamma_I + Kdot_cl(:,1:2);
            K_P_dot = adapt_on(c) * B' * P * e * x(3:5)' * Gamma_P + Kdot_cl(:,3:5);
            u = - K_I * x(1:2) - K_P * x(3:5) + W;
            %x_pre = x;  u_pre = u;
            x_pre = x;  e_pre = e;  K_pre = [K_I K_P];
            V_bar = V_bar - (e' * Q_tilde * e + diss)*h;
            x     = x + (A_true*x + B*u + B_ref*Y_cmd(:,k) + B_v*V)*h;
            x_ref = x_ref + (A_ref*x_ref + B*W + B_ref*Y_cmd(:,k) + B_v*V + L_v*e)*h;
            K_I = K_I + K_I_dot*h;
            K_P = K_P + K_P_dot*h;
        end

        % Record data: y = Theta*x_pre (xdot from measured states)
        if cl_on(c)
            %xdot = (x - x_pre)/h;
            %y = B_pinv*(xdot - A*x_pre - B*u_pre - B_ref*Y_cmd(:,k) - B_v*V);
            edot = ((x - x_ref) - e_pre)/h;
            y = K_pre*x_pre + B_pinv*(edot - (A_ref - L_v)*e_pre);   % = K_star*x_pre
            [Xs, Ys, changed] = cl_record(Xs, Ys, x_pre, y, N_max, nov_tol);
            if changed
            %     kce = kc / (1 + max(eig(Xs*Xs')));
                 if size(Xs,2) >= 5, sig_cur = min(svd(Xs)); end
             end
            % if changed && gamma*h*max(eig(Xs*Xs')) >= 2
            %     warning('CL Euler step unstable at t = %.2f s: increase n_sub.', t(k));
            % end
        end
    end
end

% True Lyapunov function (simulation only: uses K_star)
Kt_true = zeros(nT,nC);  V_true = zeros(nT,nC);
for c = 1:nC
    for k = 1:nT
        Kt = K_star - K_hist(:,:,k,c);
        ek = X(:,k,c) - X_ref(:,k,c);
        Kt_true(k,c) = norm(Kt, 'fro');
        V_true(k,c)  = ek'*P*ek + trace(Kt*Kt')/gamma;
    end
end

% Summary
lim = rad2deg(Y_max);
fprintf('||Theta||_F = %.4f,  theta_max = %.4f\n\n', norm(K_star - K_lqr,'fro'), theta_max);
fprintf('%-26s %8s %8s %8s %8s\n', '', 'dela', 'delr', 'ps', 'beta');
fprintf('%-26s %8.3f %8.3f %8.3f %8.3f\n', 'limit (deg, dps)', lim);
for c = 1:nC
    pk = [rad2deg(max(abs(U(:,:,c)),[],2)); rad2deg(max(abs(X(4,:,c)))); rad2deg(max(abs(X(3,:,c))))];
    fprintf('%-26s %8.3f %8.3f %8.3f %8.3f\n', case_names{c}, pk);
end
for c = 1:nC
    fprintf('\nChecks (%s):\n', case_names{c});
    fprintf('  max |c_i*x| / limit      : %s  (<= 1: plant inside original limits)\n', mat2str(max(abs(Y_lim(:,:,c)),[],2)'./Y_max', 4));
    fprintf('  ||K_tilde||_F  t=0 -> T  : %.4f -> %.4f\n', Kt_true(1,c), Kt_true(end,c));
    if tube_on(c)
        fprintf('  max(|c_i*e| - delta_i)   : %s deg  (<= 0: tube bound holds)\n', mat2str(rad2deg(max(abs(Y_lim(:,:,c) - Y_lim_ref(:,:,c)) - D_hist(:,:,c),[],2))', 3));
        fprintf('  min feasibility margin   : %s  (> 0: condition (F) holds)\n', mat2str(F_min(:,c)', 3));
        fprintf('  V_bar          t=0 -> T  : %.3e -> %.3e\n', Vbar_hist(1,c), Vbar_hist(end,c));
        fprintf('  V_bar - V      t=0 -> T  : %.3e -> %.3e  (constant)\n', Vbar_hist(1,c)-V_true(1,c), Vbar_hist(end,c)-V_true(end,c));
    end
end

% Plots (tube + CL case)
c = nC;
beta = rad2deg(X(3,:,c));  ps = rad2deg(X(4,:,c));  rs = rad2deg(X(5,:,c));

figure('Name', case_names{c});
subplot(2,2,1);
plot(rad2deg(U(1,:,c)), rad2deg(U(2,:,c)), 'b', 'LineWidth', 0.5); hold on; grid on;
rectangle('Position', [-lim(1) -lim(2) 2*lim(1) 2*lim(2)], 'EdgeColor', 'r', 'LineWidth', 1);
xlabel('dela (deg)'); ylabel('delr (deg)'); title('Control Operational Bounds');

subplot(2,2,3);
plot(ps, beta, 'LineWidth', 0.5); hold on; grid on;
rectangle('Position', [-lim(3) -lim(4) 2*lim(3) 2*lim(4)], 'EdgeColor', 'r', 'LineWidth', 1);
xlabel('ps (dps)'); ylabel('beta (deg)'); title('Limited Output Operational Bounds');

subplot(4,2,2); plot(t, beta, 'LineWidth', 1.5); hold on; grid on;
yline([-lim(4) lim(4)], 'k-.'); ylabel('beta, deg');
plot(t, lim(4) - rad2deg(D_hist(4,:,c)), 'r--', t, -lim(4) + rad2deg(D_hist(4,:,c)), 'r--');

subplot(4,2,4); h = plot(t, rad2deg(Y_cmd(1,:)), t, ps, t, rad2deg(X_ref(4,:,c)), '--', 'LineWidth', 1.5); grid on;
yline([-lim(3) lim(3)], 'k-.'); ylabel('ps, dps'); legend(h, 'Cmd', 'ps', 'ps_{ref}');

subplot(4,2,6); plot(t, rs, 'LineWidth', 1.5); grid on; ylabel('rs, dps');

subplot(4,2,8); h1 = plot(t, X(1,:,c), t, X(2,:,c), 'LineWidth', 1.5); grid on; ylabel('eyI');
yyaxis right; h2 = plot(t, vecnorm(X(:,:,c) - X_ref(:,:,c)), 'LineWidth', 1); ylabel('||e||');
legend([h1; h2], 'ePsI', 'eNyI', '||x - x_{ref}||'); xlabel('Time, sec');
sgtitle(case_names{c});

% K_hat(t)
figure('Name', 'Adaptive gains');
for i = 1:2
    for j = 1:5
        subplot(2,5,(i-1)*5 + j); hold on; grid on;
        for cc = 1:nC
            plot(t, squeeze(K_hist(i,j,:,cc)), 'Color', col{cc}, 'LineWidth', lw(cc));
        end
        if i == 1 && j == 1, legend(case_names, 'AutoUpdate', 'off', 'Location', 'best'); end
        yline(K_star(i,j), 'r--');
        title(sprintf('K_{hat}(%d,%d)', i, j));
        if i == 2, xlabel('Time, sec'); end
    end
end
sgtitle('Adaptive gains (red dashed = K^*)');

% Aileron and rudder
figure('Name', 'Aileron and rudder');
ylab = {'dela, deg', 'delr, deg'};
for i = 1:2
    subplot(2,1,i); hold on; grid on;
    for cc = 1:nC
        plot(t, rad2deg(U(i,:,cc)), 'Color', col{cc}, 'LineWidth', lw(cc));
    end
    yline([-lim(i) lim(i)], 'r-.');
    ylabel(ylab{i});
end
legend(case_names); xlabel('Time, sec');

% Case comparison (same uncertain plant)
figure('Name', 'Case comparison');
ylab = {'beta, deg', 'ps, dps', 'rs, dps'};
for i = 1:3
    subplot(4,1,i); hold on; grid on;
    for cc = 1:nC
        plot(t, rad2deg(X(2+i,:,cc)), 'Color', col{cc}, 'LineWidth', lw(cc));
    end
    ylabel(ylab{i});
end
subplot(4,1,1); legend(case_names, 'AutoUpdate', 'off'); yline([-lim(4) lim(4)], 'r-.');
subplot(4,1,2); plot(t, rad2deg(Y_cmd(1,:)), 'k--');
subplot(4,1,4); plot(t, VW(1:2,:,nC), 'LineWidth', 1.2); grid on;
ylabel('v (AW)'); legend('v_{ps}', 'v_{Ny}'); xlabel('Time, sec');
title(['v (AW): ' case_names{nC}]);

% Tube: |c_i*e| <= delta_i (tube + CL case)
figure('Name', 'Tube');
ylab = {'dela, deg', 'delr, deg', 'ps, dps', 'beta, deg'};
for i = 1:4
    subplot(4,1,i);
    plot(t, rad2deg(abs(Y_lim(i,:,nC) - Y_lim_ref(i,:,nC))), 'b', t, rad2deg(D_hist(i,:,nC)), 'r--', 'LineWidth', 1.2); grid on;
    ylabel(ylab{i});
end
subplot(4,1,1); legend('|c_i e|', '\delta_i'); title(['Tube bound |c_i e| \leq \delta_i (' case_names{nC} ')']);
subplot(4,1,4); xlabel('Time, sec');

% delta and delta_dot (tube cases only)
figure('Name', 'delta and delta_dot');
ylab = {'dela', 'delr', 'ps', 'beta'};
unit = {'deg', 'deg', 'dps', 'deg'};
tc = find(tube_on);
for i = 1:4
    subplot(4,2,2*i-1); hold on; grid on;
    for cc = tc
        plot(t, rad2deg(D_hist(i,:,cc)), 'Color', col{cc}, 'LineWidth', lw(cc));
    end
    ylabel(sprintf('\\delta_{%s}, %s', ylab{i}, unit{i}));
    if i == 1, title('\delta_i(t)'); legend(case_names(tc)); end
    subplot(4,2,2*i); hold on; grid on;
    for cc = tc
        plot(t, rad2deg(Ddot_hist(i,:,cc)), 'Color', col{cc}, 'LineWidth', lw(cc));
    end
    ylabel(sprintf('d\\delta_{%s}/dt, %s/s', ylab{i}, unit{i}));
    if i == 1, title('d\delta_i/dt'); end
end
subplot(4,2,7); xlabel('Time, sec');
subplot(4,2,8); xlabel('Time, sec');
sgtitle('Tube radius and its rate');

% Lateral acceleration
Ny = zeros(nC,nT);
for cc = 1:nC
    Ny(cc,:) = C_p_reg(2,:)*X(3:5,:,cc) + D_p_reg(2,:)*U(:,:,cc);
end
figure('Name', 'Lateral acceleration'); hold on; grid on;
for cc = 1:nC
    plot(t, Ny(cc,:), 'Color', col{cc}, 'LineWidth', lw(cc));
end
plot(t, Y_cmd(2,:), 'k--', 'LineWidth', 1);
xlabel('Time, sec'); ylabel('N_y, g'); legend([case_names, {'Cmd'}]);
title('Lateral acceleration');

% CL diagnostics
figure('Name', 'Concurrent learning');
subplot(3,1,1); hold on; grid on;
for cc = 1:nC
    plot(t, Kt_true(:,cc), 'Color', col{cc}, 'LineWidth', lw(cc));
end
legend(case_names, 'AutoUpdate', 'off'); yline(theta_max, 'k-.');
ylabel('||K^* - K_{hat}||_F');
subplot(3,1,2);
semilogy(t, Vbar_hist(:,2), 'Color', col{2}, 'LineWidth', 1.2); hold on; grid on;
semilogy(t, Vbar_hist(:,nC), 'r', t, V_true(:,nC), 'b', t, Vbar_hist(:,nC) - V_true(:,nC), 'k--', 'LineWidth', 1.2);
ylabel('V'); legend('V_{bar} (tube)', 'V_{bar} (tube + CL)', 'V true (tube + CL)', 'V_{bar} - V (tube + CL)');
subplot(3,1,3);
plot(t, sig_hist(:,nC), 'b', 'LineWidth', 1.2); grid on;
ylabel('\sigma_{min}(X)'); xlabel('Time, sec'); title('History-stack rank condition (> 0 needed)');
sgtitle('Concurrent learning');


%% Local functions
function [Kdot_cl, diss] = cl_term(K, Xs, Ys, K_lqr, kce, Gam)
% CL update and its exact Lyapunov dissipation (V_dot gets -diss)
    if isempty(Xs) 
        Kdot_cl = zeros(size(K));  diss = 0;  return;
    end
    %E       = K*Xs - K_lqr*Xs - Ys;          % eps_j = -K_tilde*x_j
    E       = K*Xs - Ys;                     % eps_j = K_hat*x_j - K_star*x_j = -K_tilde*x_j            % Remark 5 of paper
    Kdot_cl = -(E*Xs') * Gam;             % second part of equation 18                 
    diss    = 2*sum(E(:).^2);               % Second term of equation 21. The paper uses it only to prove convergence. This code uses it to shrink the tube
end

function [Xs, Ys, changed] = cl_record(Xs, Ys, x, y, N_max, nov_tol)
% Add novel points; when full, replace the point that maximizes sigma_min(X)
    changed = false;
    %if norm(x) < x_min, return; end         % Skip near zero states
    %if ~isempty(Xs) && norm(x - Xs(:,end))/norm(x) < nov_tol, return; end
    if ~isempty(Xs) && norm(x - Xs(:,end))^2/norm(x) < nov_tol, return; end   % eq. 15
    if size(Xs,2) < N_max
        Xs = [Xs x];  Ys = [Ys y];  changed = true;  return;
    end
    best = min(svd(Xs));  jb = 0;
    for j = 1:N_max
        Xt = Xs;  Xt(:,j) = x;
        sj = min(svd(Xt));
        if sj > best, best = sj; jb = j; end
    end
    if jb > 0
        Xs(:,jb) = x;  Ys(:,jb) = y;  changed = true;
    end
end