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
gamma   = 5e4;                               % learning rate (large enough for a feasible tube)
Gamma_I = gamma*eye(2);
Gamma_P = gamma*eye(3);

% Tube (Theorem 2): |c_i*e| <= delta_i = sqrt(c_i*P^-1*c_i' * V_bar)
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
dt = 1e-4; T = 40;
n_sub = 10;                                  % Euler sub-steps per dt (one step diverges at gamma = 1e5)
t = 0:dt:T; nT = numel(t);
pulse = zeros(1,nT);
pulse(t >= 0  & t < 10) =  1;
pulse(t >= 20 & t < 30) = -1;
Y_cmd = [deg2rad(15)*pulse; 0.0312*pulse];   % ps_cmd Ny_cmd

% Time Loop
case_names = {'MRAC + CBF, no tube', 'MRAC + CBF + tube'};
adapt_on   = [1 1];
tube_on    = [0 1];

X = zeros(5,nT,2);  X_ref = zeros(5,nT,2);  U = zeros(2,nT,2);  VW = zeros(4,nT,2);
K_hist = zeros(2,5,nT);
Y_lim = zeros(4,nT,2);  Y_lim_ref = zeros(4,nT,2);  D_hist = zeros(4,nT,2);  F_min = inf(4,2);

for c = 1:2
    x = zeros(5,1);  x_ref = zeros(5,1);
    K_I = K_lqr(:,1:2);  K_P = K_lqr(:,3:5);
    V_bar = V_bar0;

    for k = 1:nT

        e_yi = x(1:2);
        x_p = x(3:5);
        e = x - x_ref;

        % Adaptive Laws
        K_I_dot = adapt_on(c) * B' * P * e * e_yi' * Gamma_I;
        K_P_dot = adapt_on(c) * B' * P * e * x_p' * Gamma_P;

        % Baseline with (adaptive) gains
        U_bl = - K_I * e_yi - K_P * x_p;

        % CBF matrices (certainty-equivalence closed-loop model A_ref)
        C_lim = [-K_I   -K_P;
            zeros(2,2)           C_plim];
        H_x = C_lim * A_ref + alpha * C_lim;
        H_u = [-K_I        -K_I*D_p_reg - K_P*B_p;
               zeros(2,2)   H_w];
        E     = [-[K_I_dot K_P_dot]*x_ref; zeros(2,1)] + C_lim*L_v*e;

        % Tube and its time derivative (zero when tube_on = 0)
        eQe   = e' * Q_tilde * e;
        s     = diag(C_lim * P_inv * C_lim');
        s_dot = 2*diag([-[K_I_dot K_P_dot]; zeros(2,5)] * P_inv * C_lim');
        delta     = tube_on(c) * sqrt(s*V_bar);
        delta_dot = tube_on(c) * (s_dot*V_bar - s*eQe) ./ (2*sqrt(s*V_bar));

        % CBF evaluated on the reference model with tightened limits
        q = [-Y_cmd(:,k); zeros(2,1)];
        H_max_delta =  H_x*x_ref + H_u*q + E + delta_dot - alpha*(Y_max - delta);
        H_min_delta = -H_x*x_ref - H_u*q - E + delta_dot + alpha*(Y_min + delta);
        vw = H_u \ (max(0, H_min_delta) - max(0, H_max_delta));
        V = vw(1:2);  W = vw(3:4);

        u = U_bl + W;

        % Log
        X(:,k,c) = x;  X_ref(:,k,c) = x_ref;  U(:,k,c) = u;  VW(:,k,c) = vw;
        if adapt_on(c), K_hist(:,:,k) = [K_I K_P]; end
        Y_lim(:,k,c) = C_lim*x;  Y_lim_ref(:,k,c) = C_lim*x_ref;  D_hist(:,k,c) = delta;
        F_min(:,c) = min(F_min(:,c), diag(alpha).*(Y_max - Y_min - 2*delta) - 2*delta_dot);   % (F) > 0

        % dynamics: n_sub Euler sub-steps (V, W held over dt)
        h = dt/n_sub;
        for j = 1:n_sub
            e = x - x_ref;
            K_I_dot = adapt_on(c) * B' * P * e * x(1:2)' * Gamma_I;
            K_P_dot = adapt_on(c) * B' * P * e * x(3:5)' * Gamma_P;
            u = - K_I * x(1:2) - K_P * x(3:5) + W;
            V_bar = V_bar - (e' * Q_tilde * e)*h;
            x     = x + (A_true*x + B*u + B_ref*Y_cmd(:,k) + B_v*V)*h;
            x_ref = x_ref + (A_ref*x_ref + B*W + B_ref*Y_cmd(:,k) + B_v*V + L_v*e)*h;
            K_I = K_I + K_I_dot*h;
            K_P = K_P + K_P_dot*h;
        end

    end
end

% Summary
lim = rad2deg(Y_max);
fprintf('%-26s %8s %8s %8s %8s\n', '', 'dela', 'delr', 'ps', 'beta');
fprintf('%-26s %8.3f %8.3f %8.3f %8.3f\n', 'limit (deg, dps)', lim);
for c = 1:2
    pk = [rad2deg(max(abs(U(:,:,c)),[],2)); rad2deg(max(abs(X(4,:,c)))); rad2deg(max(abs(X(3,:,c))))];
    fprintf('%-26s %8.3f %8.3f %8.3f %8.3f\n', case_names{c}, pk);
end
c = 2;
fprintf('\nTube checks (%s):\n', case_names{c});
fprintf('  max |c_i*x| / limit      : %s  (<= 1: plant inside original limits)\n', mat2str(max(abs(Y_lim(:,:,c)),[],2)'./Y_max', 4));
fprintf('  max(|c_i*e| - delta_i)   : %s deg  (<= 0: tube bound holds)\n', mat2str(rad2deg(max(abs(Y_lim(:,:,c) - Y_lim_ref(:,:,c)) - D_hist(:,:,c),[],2))', 3));
fprintf('  min feasibility margin   : %s  (> 0: condition (F) holds)\n', mat2str(F_min(:,c)', 3));

%% Plots (MRAC + CBF)
c = 2;
beta = rad2deg(X(3,:,c));  ps = rad2deg(X(4,:,c));  rs = rad2deg(X(5,:,c));

figure('Name', 'MRAC + CBF');
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
sgtitle('MRAC + CBF');


% K_hat(t)
figure('Name', 'Adaptive gains');
for i = 1:2
    for j = 1:5
        subplot(2,5,(i-1)*5 + j);
        plot(t, squeeze(K_hist(i,j,:)), 'b', 'LineWidth', 1.2); hold on; grid on;
        yline(K_star(i,j), 'r--');
        title(sprintf('K_{hat}(%d,%d)', i, j));
        if i == 2, xlabel('Time, sec'); end
    end
end
sgtitle('Adaptive gains: K_I (cols 1-2), K_P (cols 3-5), K^* dashed');


% Aileron and rudder
figure('Name', 'Aileron and rudder');
ylab = {'dela, deg', 'delr, deg'};
for i = 1:2
    subplot(2,1,i);
    plot(t, rad2deg(U(i,:,1)), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0); hold on; grid on;
    plot(t, rad2deg(U(i,:,2)), 'b', 'LineWidth', 1.2);
    yline([-lim(i) lim(i)], 'r-.');
    ylabel(ylab{i});
end
legend(case_names); xlabel('Time, sec');


% No tube vs tube (same uncertain plant)
figure('Name', 'No tube vs tube');
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


% Tube: |c_i*e| <= delta_i
figure('Name', 'Tube');
ylab = {'dela, deg', 'delr, deg', 'ps, dps', 'beta, deg'};
for i = 1:4
    subplot(4,1,i);
    plot(t, rad2deg(abs(Y_lim(i,:,2) - Y_lim_ref(i,:,2))), 'b', t, rad2deg(D_hist(i,:,2)), 'r--', 'LineWidth', 1.2); grid on;
    ylabel(ylab{i});
end
subplot(4,1,1); legend('|c_i e|', '\delta_i'); title('Tube bound |c_i e| \leq \delta_i');
subplot(4,1,4); xlabel('Time, sec');
