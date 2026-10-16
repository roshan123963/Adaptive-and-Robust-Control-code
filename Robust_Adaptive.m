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

% Extended Plant

A = [zeros(2, 2) C_p_reg;
    zeros(3, 2) A_p];

B = [D_p_reg;
    B_p];

B_ref = [-eye(2);
    zeros(3,2)];

B_v = [eye(2);
    zeros(3,2)];

B_full = [-B_ref B];



% LQR
Q_lqr = diag([1.025 1.0289 0 0 1.6021]);
R_lqr = diag([1 0.49129]);
K_lqr = lqr(A,B,Q_lqr,R_lqr);

A_ref = A - B * K_lqr;

% True Plant
Theta  = 0.5 * [zeros(2,2) K_lqr(:,3:5)];
A_true = A + B*Theta;
K_star = K_lqr + Theta;


% Adaptive law
Q       = eye(5);
P       = lyap(A_ref', Q);                   % A_ref'*P + P*A_ref = -Q
P_inv   = inv(P);
gamma   = 1e4;                               % learning rate
Gamma_I = gamma*eye(2);
Gamma_P = gamma*eye(3);

theta_max = 1.5;
V_bar0    = theta_max^2 / gamma;             % V(0) bound


% Input and Output Constraints
Y_max = deg2rad([3; 2; 18; 0.5]);       % aileron, rudder, roll rate, sideslip angle
Y_min = -Y_max;
C_plim = [0 1 0;
          1 0 0];
H_w = C_plim * B_p;
alpha = diag([80 8 40 40]); 


% Command Signal
dt = 1e-4; T = 40;
t = 0:dt:T; nT = numel(t);
pulse = zeros(1,nT);
pulse(t >= 0  & t < 10) =  1;
pulse(t >= 20 & t < 30) = -1;
Y_cmd = [deg2rad(15)*pulse; 0.0312*pulse];   % ps_cmd Ny_cmd

% Time Loop
X = zeros(5,nT);  X_ref = zeros(5,nT);  U = zeros(2,nT);
K_hist = zeros(2,5,nT);  D_hist = zeros(4,nT);  Vb_hist = zeros(1,nT);

x = zeros(5,1);  x_ref = zeros(5,1);
K_I = K_lqr(:,1:2);  K_P = K_lqr(:,3:5);  % run once
int_eQe = 0;

for k=1:nT

    e_yi = x(1:2);
    x_p = x(3:5);
    e = x - x_ref;

    % Adaptive Laws
    K_I_dot = B' * P * e * e_yi' * Gamma_I;
    K_P_dot = B' * P * e * x_p' * Gamma_P;

    % Baseline with (adaptive) gains
    U_bl = - K_I * e_yi - K_P * x_p;

    % CBF matrices
    C_lim = [-K_I   -K_P;
        zeros(2,2)           C_plim];
    H_x = C_lim * A_ref + alpha * C_lim;    % A (book) = A_ref (code)
    H_u = [-K_I        -K_I*D_p_reg - K_P*B_p;    
           zeros(2,2)   H_w];
    E     = [-[K_I_dot K_P_dot]*x_ref; zeros(2,1)];

    % Tightened Limits
    V_bar =max(V_bar0 - int_eQe, 0); %
    delta = sqrt(diag(C_lim * P_inv * C_lim') * V_bar);
    delta = min(delta, 0.9*Y_max);
    Y_max_t = Y_max - delta;
    Y_min_t = Y_min + delta;

    % CBF
    q = [-Y_cmd(:,k); zeros(2,1)];
    H_max_delta =  H_x*x_ref + H_u*q + E - alpha*Y_max_t;
    H_min_delta = -H_x*x_ref - H_u*q - E + alpha*Y_min_t;
    vw = H_u \ (max(0, H_min_delta) - max(0, H_max_delta));
    V = vw(1:2);  W = vw(3:4);

    u = U_bl + W;

    % Log
    X(:,k) = x;  X_ref(:,k) = x_ref;  U(:,k) = u;
    K_hist(:,:,k) = [K_I K_P];  D_hist(:,k) = delta;  Vb_hist(k) = V_bar;



    % dynamics
    x     = x + (A_true*x + B*u + B_ref*Y_cmd(:,k) + B_v*V)*dt;
    x_ref = x_ref + (A_ref*x_ref + B*W + B_ref*Y_cmd(:,k) + B_v*V)*dt;
    K_I = K_I + K_I_dot*dt;
    K_P = K_P + K_P_dot*dt;
    int_eQe = int_eQe + (e' * Q * e) * dt; 

end

% Plots
lim = rad2deg(Y_max);
beta = rad2deg(X(3,:));  ps = rad2deg(X(4,:));  rs = rad2deg(X(5,:));

figure('Name', 'Adaptive + CBF + tightening');
subplot(2,2,1);
plot(rad2deg(U(1,:)), rad2deg(U(2,:)), 'b--', 'LineWidth', 0.5); hold on; grid on;
rectangle('Position', [-lim(1) -lim(2) 2*lim(1) 2*lim(2)], 'EdgeColor', 'r', 'LineWidth', 1);
xlabel('dela (deg)'); ylabel('delr (deg)'); title('Control Operational Bounds');

subplot(2,2,3);
plot(ps, beta, 'LineWidth', 0.5); hold on; grid on;
rectangle('Position', [-lim(3) -lim(4) 2*lim(3) 2*lim(4)], 'EdgeColor', 'r', 'LineWidth', 1);
xlabel('ps (dps)'); ylabel('beta (deg)'); title('Limited Output Operational Bounds');

subplot(4,2,2); plot(t, beta, 'LineWidth', 1.5); hold on; grid on;
yline([-lim(4) lim(4)], 'k-.'); ylabel('beta, deg');

subplot(4,2,4); h = plot(t, rad2deg(Y_cmd(1,:)), t, ps, t, rad2deg(X_ref(4,:)), '--', 'LineWidth', 1.5); grid on;
yline([-lim(3) lim(3)], 'k-.'); ylabel('ps, dps'); legend(h, 'Cmd', 'ps', 'ps_{ref}');

subplot(4,2,6); plot(t, rs, 'LineWidth', 1.5); grid on; ylabel('rs, dps');

subplot(4,2,8); plot(t, vecnorm(X - X_ref), 'LineWidth', 1.5); grid on;
ylabel('||e||'); xlabel('Time, sec');
                   


% K_hat(t)
figure('Name', 'Adaptive gains');
for i = 1:2
    for j = 1:5
        subplot(2,5,(i-1)*5 + j);
        plot(t, squeeze(K_hist(i,j,:)), 'b', 'LineWidth', 1.2); grid on;
        title(sprintf('K_{hat}(%d,%d)', i, j));
        if i == 2, xlabel('Time, sec'); end
    end
end
sgtitle('Adaptive gains: K_I (cols 1-2), K_P (cols 3-5)');


%
figure('Name', 'Aileron and rudder');

subplot(2,1,1);
plot(t, rad2deg(U(1,:)), 'b', 'LineWidth', 1.2); hold on; grid on;


subplot(2,1,2);
plot(t, rad2deg(U(2,:)), 'b', 'LineWidth', 1.2); hold on; grid on;




