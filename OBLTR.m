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


eig(A_p)
Co = ctrb(A_p, B_p);
rank(Co);

% Plant
A = [zeros(2, 2) C_p_reg;
    zeros(3, 2) A_p];

B = [D_p_reg;
    B_p];

B_ref = [-eye(2);
    zeros(3,2)];


rank(ctrb(A,B));


% LQR
Q_lqr = diag([1.025 1.0289 0 0 1.6021]);

R_lqr = diag([1 0.49129]);

K_lqr = lqr(A, B, Q_lqr, R_lqr);

eig(A)
eig(A - B * K_lqr);

% Closed Loop System
A_cl = A - B * K_lqr;
B_cl = B_ref;
C_cl = [zeros(2,2) C_p_reg] - D_p_reg * K_lqr;
D_cl = zeros(2,2);


sys_cl = ss(A_cl, B_cl, C_cl, D_cl);


sys_cl.InputName = {'ps_{cmd}', 'Ny_{cmd}'};
sys_cl.OutputName = {'ps', 'Ny'};

figure;
step(sys_cl, 8);
grid on


% C - matrix (beta is not directly measured)
C = [1 0 0 0 0;
    0 1 0 0 0;
    0 0 0 1 0;
    0 0 0 0 1];

% B bar
B2 = [-2.036 0.1124;
    -0.1531 -1.452;
    0.1482 1.406;
    -0.001833 -0.01739;
    -0.02397 -0.2274];

B_bar = [B B2];

% ARE weights
v = 0.1; R_0 = 1.4 * eye(4); Q_0 = eye(5);

Q_v = Q_0 + ((v + 1)/v) * B_bar * (B_bar)';

R_v = (v/(v+1)) * R_0;

Pv = icare(A', C', Q_v, R_v);       %(ARE)

Lv = Pv * C' / R_v;                 % Observer Gain


% combined State - Actual and estimated 
A_ob = [A               -B * K_lqr;
       Lv * C                A - B * K_lqr - Lv * C];

B_ob = [B_ref; B_ref];

C_ob = [[zeros(2, 2), C_p_reg],    -D_p_reg * K_lqr];

D_ob = zeros(2, 2);

sys_ob = ss(A_ob, B_ob, C_ob, D_ob);

sys_ob.InputName = {'ps_{cmd}', 'Ny_{cmd}'};
sys_ob.OutputName = {'ps', 'Ny'};

% Compare full state LQR and observer based LQR
figure;
step(sys_cl, 'r', sys_ob, 'g--', 8);
legend('Full State LQR', 'OBLTR');
grid on;


eig(A - Lv*C);           % error tends to zero (e = x - x_hat)


% Adaptive Law Constants
R_0_half = R_0^(-0.5);

M = B_bar' * C' * R_0_half;

[U, Sig, V] = svd(M);

W = (U * V')';

Gamma = 10 * norm(Lv);      % Adaptation Rate

S = [eye(2) zeros(2)];


% Regressor Phi
c_rbf = deg2rad(-4:1:4);            % 9 centers
sig_rbf = deg2rad(4/3);
Phi = @(xh) exp(-(xh(3) - c_rbf').^2 / (2 * (sig_rbf)^2));     % xh = beta


% uncertainty model

Lambda = 0.1 * eye(2);

f_beta = @(b) 1 - exp(-(b/0.0233)^2); 

Kp = K_lqr(:, 3:5);  % proportional gain
%d = @(x) Kp * x(3:5) + f_beta(x(3)) * [1; 1];
d = @(x) f_beta(x(3)) * [1; 1];



% Time loop simulation
dt = 0.003; T = 40;
t = 0:dt:T; nT = numel(t);

%ycmd = [deg2rad(10)*((t>1 & t<11) - (t>21 & t<31));   % Commands: ps pulses +/-10 deg/s, Ny pulse 0.1 g
        %0.1*(t>6 & t<16)];

 ycmd = [deg2rad(10)*((t>1 & t<11)- (t>21 & t<31));
 0.0*(t>6 & t<16)];

unc_on = [0 1 1];           % nominal, uncertain + baseline, uncertain + baseline +adaptive
adapt_on = [0 0 1];

Y = zeros(2, nT, 3);        % row = ps, Ny
Yh = zeros(2, nT, 3);       % observer's estimated [ps; Ny]
U = zeros(2, nT, 3);        % control input [aileron; rudder]

for c = 1:3
    Lam = eye(2) - unc_on(c)* 0.9 * eye(2);

    x = zeros(5,1);
    xh = zeros(5,1);
    Th = zeros(9,2);
     

    for k = 1:nT
        u_bl = -K_lqr * xh;
        ph = Phi(xh);
        u = u_bl - adapt_on(c)* (Th' * ph);
        dk = unc_on(c) * d(x);
        e_y = C * (xh - x);
    
        Y(:,k,c) = C_p_reg*x(3:5) + D_p_reg*Lam*(u + dk);
        Yh(:,k,c) = C_p_reg*xh(3:5) + D_p_reg*u_bl;
        
        x_dot = A * x + B * Lam * (u + dk) + B_ref * ycmd(:,k);
        xh_dot = A * xh + B * u_bl + Lv * (C * x - C * xh) + B_ref * ycmd(:,k);
        Th_dot = -Gamma * ph * e_y' * R_0_half * W * S';
    
        
        x = x + x_dot * dt;
        xh = xh + xh_dot * dt;
        Th = Th + adapt_on(c) * Th_dot * dt;
        U(:,k,c) = u;
    
     end
end


%Plot
figure;
names = {'Nominal', 'Baseline + unc', 'Adaptive + unc'};
%names = {'Nominal', 'Baseline + unc', 'Adaptive + unc', ...
         %'Estimated (baseline + unc)', 'Estimated (adaptive + unc)'};
%names = {'Nominal', 'Baseline + unc', 'Adaptive + unc', 'Estimated (adaptive + unc)'};
for i = 1:2
    subplot(2,1,i); hold on; grid on
    plot(t, ycmd(i,:), 'k');
    plot(t, Y(i,:,1), 'b', t, Y(i,:,2), 'r', t, Y(i,:,3), 'g--', 'LineWidth', 1.5);
    %plot(t, Yh(i,:,2), 'm:', t, Yh(i,:,3), 'k-.', 'LineWidth', 1.5);
    %plot(t, Yh(i,:,3), 'k-.', 'LineWidth', 1.5);
end
subplot(2,1,1); ylabel('p_s (rad/s)'); legend(['Command', names]);
subplot(2,1,2); ylabel('N_y (g)');     xlabel('t (s)');



figure;
ctrl = {'\delta_{ail} (deg)', '\delta_{rud} (deg)'};
for i = 1:2
    subplot(2,1,i); hold on; grid on
    plot(t, rad2deg(U(i,:,1)), 'b', t, rad2deg(U(i,:,2)), 'r', ...
         t, rad2deg(U(i,:,3)), 'g--', 'LineWidth', 1.5);
    ylabel(ctrl{i});
end
subplot(2,1,1); legend(names);
subplot(2,1,2); xlabel('t (s)');






