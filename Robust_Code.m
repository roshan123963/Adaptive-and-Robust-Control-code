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

K_I = K_lqr(:,1:2);
K_P = K_lqr(:,3:5);


% Input and Output Contraints
Y_max = deg2rad([3; 2; 18; 0.5])*1;


Y_min = -Y_max;


C_plim = [0 1 0;
         1 0 0];


C_lim = [-K_I  -K_P;
        zeros(2,2)     C_plim];

B_full = [-B_ref B];


% CBF Design Matrices
alpha = diag([80 8 40 40]); 

H_x = C_lim * A + alpha * C_lim;

% H_u = C_lim * B_full;

H_w = C_plim * A_p * B_p;
H_u = [-K_I -K_I*D_p_reg-K_P*B_p;
         zeros(2,2) H_w];

% Command Signal

dt = 0.0001; T = 40;
t = 0:dt:T; nT = numel(t);

pulse = zeros(1,nT);
pulse(t >= 0 & t < 10) = 1;
pulse(t  >= 20 & t < 30) = -1;


ps_cmd = deg2rad(15) * pulse;

Ny_cmd = 0.0312 * pulse;

Y_cmd = [ps_cmd; Ny_cmd];

% Time loop
case_names = {'Baseline (no saturation)', 'Baseline + hard saturation', 'Baseline + CBF augmentation'};

X = zeros(5,nT,3);
U = zeros(2,nT,3);

for c = 1:3
    x = zeros(5,1);
    for k = 1:nT    
        e_yi = x(1:2);
        x_p = x(3:5);
    
        U_bl = - K_I * e_yi - K_P * x_p;
        V = zeros(2,1);

        if c == 1
            u = U_bl;       % no saturation and no augmentation

        elseif c == 2
            u = U_bl;       % hard saturation

            % Aileron Limit = +/- 3 

            if u(1) > deg2rad(3)
                u(1) = deg2rad(3);
            elseif u(1) < -deg2rad(3)
                u(1) = -deg2rad(3);
            end 

            % Rudder Limit = +/- 2
            if u(2) > deg2rad(2)
                u(2) = deg2rad(2);
            elseif u(2) < -deg2rad(2)
                u(2) = -deg2rad(2);
            end 

        else                % CBF/AW augmentation
            H_min_delta = - H_x * x - H_u * [-Y_cmd(:,k); U_bl] + alpha * Y_min;
            H_max_delta = H_x * x + H_u * [-Y_cmd(:,k); U_bl] - alpha * Y_max;
        
            vw = (H_u) \ (max(0, H_min_delta) - max(0, H_max_delta));
            V = vw(1:2);
            W = vw(3:4);
        
            u = U_bl + W;
        end 
    
        X(:,k,c) = x;
        U(:,k,c) = u;


        e_yi_dot = C_p_reg * x_p -Y_cmd(:,k) + V + D_p_reg * u;
        x_dot_p = A_p * x_p + B_p * u;    
        x = x + [e_yi_dot; x_dot_p] * dt;
    end
end 

lim_deg = rad2deg(Y_max);

for c = 1:3
    beta = rad2deg(X(3,:,c)); ps = rad2deg(X(4,:,c)); rs = rad2deg(X(5,:,c));
    da = rad2deg(U(1,:,c)); dr = rad2deg(U(2,:,c));
    figure('Name', case_names{c});
    subplot(2,2,1);
    plot(da,dr, 'LineWidth', 0.5); 
    hold on; grid on;
    rectangle('Position', [-lim_deg(1) -lim_deg(2) 2*lim_deg(1) 2*lim_deg(2)], 'EdgeColor', 'red', 'LineWidth', 1);
    xlabel('dela(deg)'); ylabel('delr(deg)'); title('Control Operational Bounds');

    subplot(2,2,3); plot(ps, beta, 'LineWidth', 0.5); hold on; grid on;
    rectangle('Position', [-lim_deg(3) -lim_deg(4) 2*lim_deg(3) 2*lim_deg(4)], 'EdgeColor', 'r', 'LineWidth', 1);
    xlabel('ps(dps)'); ylabel('beta(deg)'); title('Limited Output Operational Bounds');

    subplot(4,2,2); plot(t, beta, 'LineWidth', 1.5); hold on; grid on;
    yline([-0.5 0.5], 'k-.'); ylabel('beta, deg'); title('States and Tracking Error Integrators');

    subplot(4,2,4); h = plot(t, rad2deg(ps_cmd), t, ps, 'LineWidth', 1.5); hold on; grid on;
    yline([-18 18], 'k-.'); ylabel('ps, dps'); legend(h, 'Cmd', 'ps');

    subplot(4,2,6); plot(t, rs, 'LineWidth', 1.5); grid on; ylabel('rs, dps');

    subplot(4,2,8); plot(t, X(1,:,c), t, X(2,:,c), 'LineWidth', 1.5); grid on;
    ylabel('eyI'); xlabel('Time, sec'); legend('ePsI', 'eNyI');

    sgtitle(case_names{c});

end 















