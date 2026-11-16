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

% Design model. The true A (and A_p, C_p_reg) is unknown to the controller. This
% nominal model is used only at design time, to choose A_ref and the initial gains
% K_hat(0) = K_lqr. The controller itself uses only A_ref, B, B_ref, B_v and P.
A_nom = [zeros(2, 2) C_p_reg;
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
K_lqr = lqr(A_nom,B,Q_lqr,R_lqr);

A_ref = A_nom - B * K_lqr;                   % reference model (Hurwitz), known

% True Plant: A_true is unknown. The only assumption is matching: there is an
% unknown K_star with A_true - B*K_star = A_ref, and ||K_star - K_hat(0)||_F <= theta_max.
% Every A_true admitted by matching has the form A_ref + B*K_star. K_star below picks
% one such plant for the simulation; the controller never uses A_true or K_star.
K_star = K_lqr + 0.5 * [zeros(2,2) K_lqr(:,3:5)];
A_true = A_ref + B*K_star;


% Reference model:  x_ref_dot = A_ref*x_ref + ... + L_v*(x - x_ref)
nu  = 0.1;
Q_v = eye(5) + ((nu + 1)/nu) * eye(5);
R_v = (nu/(nu + 1)) * eye(5);
P_v = icare(A_ref', eye(5), Q_v, R_v);      % A_ref*P_v + P_v*A_ref' - P_v*R_v^-1*P_v + Q_v = 0
% L_v = P_v / R_v;  P = inv(P_v);           % observer-like reference model
L_v = zeros(5);                              % open-loop reference model
P   = lyap(A_ref', eye(5));                  % (A_ref - L_v)'*P + P*(A_ref - L_v) = -I
P_inv = inv(P);

% Adaptive law: K_hat_dot = B'*P*e*x'*Gamma
gamma   = 5e4;                               % learning rate (large enough for feasible tightened limits)
Gamma_I = gamma*eye(2);
Gamma_P = gamma*eye(3);

% Constraint tightening (Theorem 2): the error ellipsoid e'*P*e <= V_bar gives the
% safety margin |c_i*e| <= delta_i = sqrt(c_i*P^-1*c_i' * V_bar)
theta_max = 1.5;                             % known bound on ||K_star - K_hat(0)||_F
V_bar0    = theta_max^2 / gamma;             % V(0) bound (e(0) = 0)
Q_tilde   = -((A_ref - L_v)'*P + P*(A_ref - L_v));   % V_dot = -e'*Q_tilde*e


% Input and Output Constraints
Y_max = deg2rad([3; 2; 18; 0.5]);       % aileron, rudder, roll rate, sideslip angle
Y_min = -Y_max;
C_plim = [0 1 0;
          1 0 0];
H_w = C_plim * B_p;
alpha = diag([80 8 40 4]);


% Command Signal
% The N_y command must be reachable inside the limits, otherwise no controller can
% track it. 0.0312 g at ps = 15 dps needs a trim beta of -0.95 deg (limit 0.5 deg) and
% rudder of -2.87 deg (limit 2 deg). 0.005 g needs beta = -0.07 deg and rudder = -0.5 deg,
% which also fits inside the tightened limits (|beta| <= 0.5 - delta_beta ~ 0.15 deg).
Ny_amp = 0.005;
dt = 1e-4; T = 40;
t = 0:dt:T; nT = numel(t);
pulse = zeros(1,nT);
pulse(t >= 0  & t < 10) =  1;
pulse(t >= 20 & t < 30) = -1;
Y_cmd = [deg2rad(15)*pulse; Ny_amp*pulse];   % ps_cmd Ny_cmd

trim = [A_true(3:5,3:5) B_p; A_true(1:2,3:5) D_p_reg] \ [zeros(3,1); deg2rad(15); Ny_amp];   % true plant: simulation check only
fprintf('Trim for the command (true plant): beta = %.3f deg, dela = %.3f deg, delr = %.3f deg\n\n', rad2deg(trim([1 4 5])));

% Parameters for the closed-loop right-hand side
prm.A_true = A_true;  prm.A_ref = A_ref;  prm.B = B;  prm.B_ref = B_ref;  prm.B_v = B_v;
prm.L_v = L_v;  prm.P = P;  prm.P_inv = P_inv;  prm.Q_tilde = Q_tilde;
prm.Gamma = blkdiag(Gamma_I, Gamma_P);  prm.C_plim = C_plim;  prm.B_pi = [B_v B];
prm.alpha = alpha;  prm.Y_max = Y_max;  prm.Y_min = Y_min;

% Time Loop
% The whole closed loop (plant, reference model, adaptive law, V_bar and the CBF) is one
% continuous-time ODE, integrated with RK4. The CBF is evaluated inside every RK4 stage,
% so the signals (v, w) are a continuous feedback law, as in the proofs, not held
% between steps.
case_names = {'MRAC + CBF, no tightening', 'MRAC + CBF + tightening'};
adapt_on   = [1 1];
tighten_on = [0 1];

X = zeros(5,nT,2);  X_ref = zeros(5,nT,2);  U = zeros(2,nT,2);  U_bl = zeros(2,nT,2);  VW = zeros(4,nT,2);
K_hist = zeros(2,5,nT,2);  Y_out = zeros(2,nT,2);  Y_out_ref = zeros(2,nT,2);
Y_lim = zeros(4,nT,2);  Y_lim_ref = zeros(4,nT,2);  D_hist = zeros(4,nT,2);  Dd_hist = zeros(4,nT,2);
AB_hist = zeros(4,nT,2);  Vbar_hist = zeros(1,nT,2);  V_hist = zeros(1,nT,2);

for c = 1:2
    prm.adapt = adapt_on(c);  prm.tighten = tighten_on(c);
    z = [zeros(5,1); zeros(5,1); reshape(K_lqr,[],1); V_bar0];   % x(0) = x_ref(0) = 0 (A3)

    for k = 1:nT
        yc = Y_cmd(:,k);
        [k1, o] = cl_rhs(z, yc, prm);

        % Log
        K_hat = reshape(z(11:20), 2, 5);
        X(:,k,c) = o.x;  X_ref(:,k,c) = o.x_ref;  U(:,k,c) = o.u;  U_bl(:,k,c) = o.u_bl;  VW(:,k,c) = o.vw;
        K_hist(:,:,k,c) = K_hat;  Y_out(:,k,c) = o.y;  Y_out_ref(:,k,c) = o.y_ref;
        Y_lim(:,k,c) = o.C_lim*o.x;  Y_lim_ref(:,k,c) = o.C_lim*o.x_ref;
        D_hist(:,k,c) = o.delta;  Dd_hist(:,k,c) = o.delta_dot;  AB_hist(:,k,c) = o.a + o.b;
        Vbar_hist(1,k,c) = z(21);
        dK = K_hat - K_star;                                          % for verification only
        V_hist(1,k,c) = o.e'*P*o.e + trace(dK/prm.Gamma*dK');

        if k == nT, break; end
        k2 = cl_rhs(z + dt/2*k1, yc, prm);
        k3 = cl_rhs(z + dt/2*k2, yc, prm);
        k4 = cl_rhs(z + dt*k3,   yc, prm);
        z  = z + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    end
end

% Summary
lim = rad2deg(Y_max);
fprintf('%-26s %8s %8s %8s %8s\n', 'peaks', 'dela', 'delr', 'ps', 'beta');
fprintf('%-26s %8.3f %8.3f %8.3f %8.3f\n', 'limit (deg, dps)', lim);
for c = 1:2
    fprintf('%-26s %8.3f %8.3f %8.3f %8.3f   (u_bl, ps, beta)\n', case_names{c}, rad2deg(max(abs(Y_lim(:,:,c)),[],2)));
    fprintf('%-26s %8.3f %8.3f   (total command u_bl + w, not covered by Theorem 2)\n', '', rad2deg(max(abs(U(:,:,c)),[],2)));
end

c = 2;  tol = 1e-9;
ref_margin = max(Y_lim_ref(:,:,c) - (Y_max - D_hist(:,:,c)), (Y_min + D_hist(:,:,c)) - Y_lim_ref(:,:,c));
pl_margin  = max(Y_lim(:,:,c) - Y_max, Y_min - Y_lim(:,:,c));
margin_gap = abs(Y_lim(:,:,c) - Y_lim_ref(:,:,c)) - D_hist(:,:,c);
claims = {
    'Lemma 3: V <= V_bar',                         max(V_hist(1,:,c) - Vbar_hist(1,:,c)),       sprintf('%.2e', max(V_hist(1,:,c) - Vbar_hist(1,:,c)))
    'Lemma 3(a): |c_i e| <= delta_i',              max(margin_gap(:)),                          [mat2str(rad2deg(max(margin_gap,[],2))',3) ' deg']
    '(F): a_i + b_i < 0',                          max(max(AB_hist(:,:,c)))+tol,                mat2str(max(AB_hist(:,:,c),[],2)',3)
    'delta_dot <= 0 on ps, beta rows',             max(max(Dd_hist(3:4,:,c))),                  mat2str(max(Dd_hist(3:4,:,c),[],2)',3)
    'Claim 1: x_ref in tightened box',         max(ref_margin(:))-tol,                      [mat2str(rad2deg(max(ref_margin,[],2))',3) ' deg']
    'Claim 2: plant in original box',              max(pl_margin(:)),                           [mat2str(rad2deg(max(pl_margin,[],2))',3) ' deg']
    };
fprintf('\nTheorem 2 checks (%s). Each value must be <= 0:\n', case_names{c});
for i = 1:size(claims,1)
    st = 'FAIL'; if claims{i,2} <= 0, st = 'ok'; end
    fprintf('  %-34s %-4s  %s\n', claims{i,1}, st, claims{i,3});
end
fprintf('  (aileron/rudder delta_dot > 0 on %.0f%% / %.0f%% of samples: allowed, c_i depends on K_hat)\n', ...
    100*mean(Dd_hist(1,:,c) > 0), 100*mean(Dd_hist(2,:,c) > 0));
fprintf('  V_bar(T)/V_bar(0) = %.3f,  ||e|| peak %.2e -> final %.2e\n', Vbar_hist(1,end,c)/V_bar0, ...
    max(vecnorm(X(:,:,c) - X_ref(:,:,c))), norm(X(:,end,c) - X_ref(:,end,c)));

fprintf('\nTracking error at the end of each command segment (ps in dps, Ny in g):\n');
for tt = [10 20 30 40]
    k = find(t < tt - dt/2, 1, 'last');
    for c = 1:2
        fprintf('  t = %5.2f  %-22s  ps - cmd = %+8.4f   Ny - cmd = %+9.5f\n', t(k), case_names{c}, ...
            rad2deg(Y_out(1,k,c) - Y_cmd(1,k)), Y_out(2,k,c) - Y_cmd(2,k));
    end
end

%% Plots (MRAC + CBF + tightening)
c = 2;
beta = rad2deg(X(3,:,c));  ps = rad2deg(X(4,:,c));  rs = rad2deg(X(5,:,c));

figure('Name', 'MRAC + CBF');
subplot(2,2,1);
plot(rad2deg(U(1,:,c)), rad2deg(U(2,:,c)), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); hold on; grid on;
plot(rad2deg(U_bl(1,:,c)), rad2deg(U_bl(2,:,c)), 'b', 'LineWidth', 0.5);
rectangle('Position', [-lim(1) -lim(2) 2*lim(1) 2*lim(2)], 'EdgeColor', 'r', 'LineWidth', 1);
xlabel('dela (deg)'); ylabel('delr (deg)'); title('Control Operational Bounds');
legend('u_{bl} + w', 'u_{bl} (guaranteed)', 'Location', 'best');

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
sgtitle('MRAC + CBF + constraint tightening');


% Command tracking: y = [ps; Ny] of the true plant vs command and reference model
figure('Name', 'Command tracking');
subplot(3,1,1);
plot(t, rad2deg(Y_cmd(1,:)), 'k--', t, rad2deg(Y_out(1,:,1)), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0); hold on; grid on;
plot(t, rad2deg(Y_out(1,:,2)), 'b', t, rad2deg(Y_out_ref(1,:,2)), 'r:', 'LineWidth', 1.2);
ylabel('ps, dps'); legend('Cmd', case_names{1}, case_names{2}, 'ref model (tightening)', 'Location', 'best');
subplot(3,1,2);
plot(t, Y_cmd(2,:), 'k--', t, Y_out(2,:,1), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0); hold on; grid on;
plot(t, Y_out(2,:,2), 'b', t, Y_out_ref(2,:,2), 'r:', 'LineWidth', 1.2);
ylabel('N_y, g');
subplot(3,1,3);
plot(t, VW(1:2,:,2), 'LineWidth', 1.0); grid on;
ylabel('v (tightening)'); legend('v_{ps}', 'v_{Ny}'); xlabel('Time, sec');
sgtitle('Command tracking (v \neq 0 while a limit is active)');


% K_hat(t)
figure('Name', 'Adaptive gains');
for i = 1:2
    for j = 1:5
        subplot(2,5,(i-1)*5 + j);
        plot(t, squeeze(K_hist(i,j,:,2)), 'b', 'LineWidth', 1.2); hold on; grid on;
        yline(K_star(i,j), 'r--');
        title(sprintf('K_{hat}(%d,%d)', i, j));
        if i == 2, xlabel('Time, sec'); end
    end
end
sgtitle('Adaptive gains (tightening case): K_I (cols 1-2), K_P (cols 3-5), K^* dashed');


% Aileron and rudder: u_bl is guaranteed by Theorem 2, the total command u_bl + w is not
figure('Name', 'Aileron and rudder');
ylab = {'dela, deg', 'delr, deg'};
for i = 1:2
    subplot(2,1,i);
    plot(t, rad2deg(U(i,:,2)), 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on; grid on;
    plot(t, rad2deg(U_bl(i,:,2)), 'b', 'LineWidth', 1.2);
    yline([-lim(i) lim(i)], 'r-.');
    plot(t, lim(i) - rad2deg(D_hist(i,:,2)), 'r--', t, -lim(i) + rad2deg(D_hist(i,:,2)), 'r--');
    ylabel(ylab{i});
end
legend('u_{bl} + w', 'u_{bl}', 'limit', 'tightened limit'); xlabel('Time, sec');
sgtitle('Constraint tightening case');


% No tightening vs tightening (same uncertain plant)
figure('Name', 'No tightening vs tightening');
ylab = {'beta, deg', 'ps, dps', 'rs, dps'};
for i = 1:3
    subplot(4,1,i);
    plot(t, rad2deg(X(2+i,:,1)), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2); hold on; grid on;
    plot(t, rad2deg(X(2+i,:,2)), 'b', 'LineWidth', 1.2);
    ylabel(ylab{i});
end
subplot(4,1,1); yline([-lim(4) lim(4)], 'r-.'); legend(case_names);
subplot(4,1,2); plot(t, rad2deg(Y_cmd(1,:)), 'k--');
subplot(4,1,4); plot(t, VW(1:2,:,2), 'LineWidth', 1.2); grid on;
ylabel('v (AW)'); legend('v_{ps}', 'v_{Ny}'); xlabel('Time, sec');


% Safety margin: |c_i*e| <= delta_i, and V <= V_bar
figure('Name', 'Safety margin');
ylab = {'dela, deg', 'delr, deg', 'ps, dps', 'beta, deg'};
for i = 1:4
    subplot(5,1,i);
    plot(t, rad2deg(abs(Y_lim(i,:,2) - Y_lim_ref(i,:,2))), 'b', t, rad2deg(D_hist(i,:,2)), 'r--', 'LineWidth', 1.2); grid on;
    ylabel(ylab{i});
end
subplot(5,1,1); legend('|c_i e|', '\delta_i'); title('Safety margin |c_i e| \leq \delta_i');
subplot(5,1,5); plot(t, V_hist(1,:,2), 'b', t, Vbar_hist(1,:,2), 'r--', 'LineWidth', 1.2); grid on;
ylabel('V'); legend('V', 'V_{bar}'); xlabel('Time, sec');


% Margin rate per row
figure('Name', 'Margin rate');
for i = 1:4
    subplot(4,1,i);
    plot(t, rad2deg(Dd_hist(i,:,2)), 'b', 'LineWidth', 1.0); hold on; grid on; yline(0, 'k-');
    ylabel(['d\delta/dt ' ylab{i} '/s']);
end
subplot(4,1,1); title('Margin rate d\delta/dt: \leq 0 on ps and \beta rows; either sign on aileron and rudder rows (c_i = -K_{hat,i})');
subplot(4,1,4); xlabel('Time, sec');


%% Closed-loop right-hand side
function [dz, o] = cl_rhs(z, y_cmd, p)
% State z = [x; x_ref; K_hat(:); V_bar]. Returns dz/dt and the signals used for logging.
x = z(1:5);  x_ref = z(6:10);  K = reshape(z(11:20), 2, 5);  V_bar = z(21);
e = x - x_ref;

% Adaptive law: K_hat_dot = B'*P*e*x'*Gamma
K_dot = p.adapt * (p.B' * p.P * e) * (x' * p.Gamma);

% Constrained outputs y_lim = C_lim*x = [u_bl; ps; beta]
C_lim = [-K; zeros(2,2) p.C_plim];
C_dot = [-K_dot; zeros(2,5)];
H_x = C_lim*p.A_ref + p.alpha*C_lim;
H_u = C_lim*p.B_pi;

% Safety margin delta_i = sqrt(c_i*P^-1*c_i' * V_bar) and its rate
eQe = e' * p.Q_tilde * e;
if p.tighten
    s     = sum((C_lim*p.P_inv) .* C_lim, 2);
    s_dot = 2*sum((C_dot*p.P_inv) .* C_lim, 2);
    delta     = sqrt(s*V_bar);
    delta_dot = (s_dot*V_bar - s*eQe) ./ (2*delta);
else
    delta = zeros(4,1);  delta_dot = zeros(4,1);
end

% CBF on the reference model with tightened limits (Theorem 2, Lemma 1)
pi0 = [-y_cmd; zeros(2,1)];
ell = H_x*x_ref + H_u*pi0 + C_dot*x_ref + C_lim*p.L_v*e;
a =  ell + delta_dot - p.alpha*(p.Y_max - delta);
b = -ell + delta_dot + p.alpha*(p.Y_min + delta);
vw = H_u \ (max(0, b) - max(0, a));
v = vw(1:2);  w = vw(3:4);

u_bl = -K*x;
u = u_bl + w;
x_dot     = p.A_true*x + p.B*u + p.B_ref*y_cmd + p.B_v*v;
x_ref_dot = p.A_ref*x_ref + p.B*w + p.B_ref*y_cmd + p.B_v*v + p.L_v*e;
dz = [x_dot; x_ref_dot; K_dot(:); -eQe];

if nargout > 1
    o = struct('x', x, 'x_ref', x_ref, 'e', e, 'u', u, 'u_bl', u_bl, 'vw', vw, 'C_lim', C_lim, ...
        'delta', delta, 'delta_dot', delta_dot, 'a', a, 'b', b, ...
        'y', x_dot(1:2) - p.B_ref(1:2,:)*y_cmd - v, ...                % measured [ps; Ny] of the true plant
        'y_ref', x_ref_dot(1:2) - p.B_ref(1:2,:)*y_cmd - v - p.L_v(1:2,:)*e);
end
end
