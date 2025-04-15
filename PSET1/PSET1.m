%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code does the programming part of the Problem Set from
% Nakamura and Steinsson (2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. CALCULATING a,b,c,d BY GUESS AND VERIFY
% -------------------------------------------------------------------------

% Define parameters
beta  = 0.99;    % Discount factor
kappa = 0.13;    % Phillips Curve slope
rho   = 0.8;     % AR(1) coefficient for M_t
yStar = 0;       % Natural output

% Define the system of equations
func = @(x) [
    -(x(1) - 1) + (-beta * (x(1) - 1) * x(1) + kappa * x(1));                                             % eq1
    (rho - x(2)) + beta * (x(1) - 1) * x(2) - beta * (rho - x(2)) * rho - kappa * x(2);                   % eq2
    (1 - x(3)) - kappa * x(3);                                                                           % eq3
    -x(4) - kappa * (x(4) - 1)                                                                           % eq4
];

% Initial guess: [a, b, c, d]
x0 = [0.5, 0.5, 0.5, 0.5];

% Solve using fsolve
options = optimoptions('fsolve', 'Display', 'iter');
solution = fsolve(func, x0, options);

% Display results
a = solution(1);
b = solution(2);
c = solution(3);
d = solution(4);

fprintf('Solved coefficients:\n');
fprintf('a = %.6f\n', a);
fprintf('b = %.6f\n', b);
fprintf('c = %.6f\n', c);
fprintf('d = %.6f\n', d);


%% 2. TRUE IMPULSE RESPONSE (THEORETICAL) 
% -------------------------------------------------------------------------

% Horizon for IRF
H = 20;  % Number of periods to compute the IRF

% Storage for IRFs
irf_y  = zeros(H+1,1);  % IRF of y
irf_pi = zeros(H+1,1);  % IRF of pi
irf_dm = zeros(H+1,1);  % IRF of Delta M_t = M_t - M_{t-1}

% Initial conditions
Mold = 0;
yold = 0;
deltaMold = 0;  % For t=0

% Apply a one-time monetary shock at t=0
eps_shock = 1;  % Shock to epsilon_t (monetary shock) - assumed
eta_shock = 0;  % No productivity shock for now

% Period t=0
M0 = rho * Mold + eps_shock;
deltaM0 = M0 - Mold;
y0 = a * yold + b * deltaMold + c * eps_shock + d * eta_shock;

irf_y(1)  = y0;
irf_dm(1) = deltaM0;
irf_pi(1) = kappa * (y0 - yStar);

% Forward for t=1,...,H
Mprev = M0;
yprev = y0;
deltaMprev = deltaM0;

for t = 1:H
    % No new shocks after t=0
    eps_now = 0;
    eta_now = 0;
    
    Mnow = rho * Mprev + eps_now;
    
    ynow = a * yprev + b * Mprev + c * eps_now + d * eta_now;
    pinow = kappa * (ynow - yStar);
    
    % Store IRFs
    irf_y(t+1)  = ynow;
    irf_dm(t+1) = Mnow;
    irf_pi(t+1) = pinow;
    
    % Update
    Mprev = Mnow;
    yprev = ynow;
end

%% Plot the IRFs
figure('Name','True IRFs - Your Model','Position',[100 100 1200 400]);

subplot(1,3,1);
plot(0:H, irf_y, '-o','LineWidth',1.5);
title('IRF of y to Monetary Shock');
xlabel('Horizon'); ylabel('Response');

subplot(1,3,2);
plot(0:H, irf_dm, '-o','LineWidth',1.5);
title('IRF of \Delta M to Monetary Shock');
xlabel('Horizon'); ylabel('Response');

subplot(1,3,3);
plot(0:H, irf_pi, '-o','LineWidth',1.5);
title('IRF of \pi to Monetary Shock');
xlabel('Horizon'); ylabel('Response');

% Save figure
saveas(gcf, 'IRFs_Q2_TrueModel.png');


%% 3. SIMULATE DATA AND ESTIMATE MISSPECIFIED MODELS
% --------------------------------------------------

clearvars -except a b c d beta kappa rho yStar irf_y;
rng(12345,'twister');  % for reproducibility - seed

% For IRFs or references
H_irf = 20;   % horizon for IRF

T      = 500;     % Length of each simulated series
Nrep   = 500;     % Number of Monte Carlo replications
Tburn  = 100;     % Burn-in to reduce influence of initial conditions
sigma_e = 0.00066;  % std dev of monetary shock
sigma_n = 0.007;    % std dev of productivity shock

%% 1. SIMULATION SETTINGS
T      = 500;    % Sample size (after burn-in)
Tburn  = 100;    % Burn-in period length
Nrep   = 500;    % Number of Monte Carlo replications

% Preallocate matrices to store IRFs for each model (each column is one replication)
IRF_model1 = zeros(H_irf+1, Nrep);  % For Model 1 (1 lag)
IRF_model2 = zeros(H_irf+1, Nrep);  % For Model 2 (4 lags)
IRF_model3 = zeros(H_irf+1, Nrep);  % For Model 3 (12 lags)

%% 2. MONTE CARLO SIMULATIONS
for irep = 1:Nrep
    %% 2a. Generate Shocks for the entire simulation sample (burn-in+T)
    eps_draw = sigma_e * randn(Tburn + T, 1);  % Shocks for M_t
    eta_draw = sigma_n * randn(Tburn + T, 1);    % Shocks for y_t (eta)
    
    %% 2b. Simulate the True Model
    % Preallocate arrays for M_t and y_t
    M = zeros(Tburn+T, 1);
    y = zeros(Tburn+T, 1);
    
    % Initialize first two values.
    M(1) = 0;
    M(2) = rho * M(1) + eps_draw(2);
    y(1) = 0;
    % For t = 2, assume M(0)=0 so that Delta M(1)=M(1)-0.
    y(2) = a * y(1) + b * (M(1) - 0) + c * eps_draw(2) + d * eta_draw(2);
    
    % Start simulation at t = 3 to avoid indexing below 1.
    for t = 3:(Tburn+T)
        M(t) = rho * M(t-1) + eps_draw(t);
        %dM_tminus1 = M(t-1) - M(t-2); % Compute Delta M_{t-1}
        y(t) = a * y(t-1) + b * M(t-1) + c * eps_draw(t) + d * eta_draw(t);
    end
    
    % Discard burn-in period.
    M = M(Tburn+1:end);
    y = y(Tburn+1:end);
    
    %% 2c. Construct eps_t from simulated M: eps_t = M_t - rho*M_{t-1}.
    eps_data = zeros(T,1);
    for t = 2:T
       eps_data(t) = M(t) - rho * M(t-1);
    end
    % (We will be using observations starting at t=2 whenever we include lagged values.)
    
    %% 2d. Estimate the three misspecified models using OLS with pinv.
    % --- Model 1: y_t = constant + y_{t-1} + eps_t ---
    Y_reg1 = y(2:end);  % Use observations 2:T for y_t.
    X_reg1 = [ones(T-1,1), y(1:end-1), eps_data(2:end)];
    b_hat1 = pinv(X_reg1) * Y_reg1;  % b_hat1 = [alpha; phi1; gamma]
    
    % --- Compute the IRF for Model 1 inline ---
    % For a one-time eps shock: we set eps_1=1 (only in period 1) and 0 otherwise.
    % Let the estimated equation be: y_t = alpha + phi1*y_{t-1} + gamma*eps_t.
    irf_temp = zeros(H_irf+1,1);
    % At t = 1, since lags are zero:
    irf_temp(1) = b_hat1(1) + b_hat1(end);  % = alpha + gamma
    % For t = 2 to H+1, no shock after period 1 (eps=0).
    for t = 2:(H_irf+1)
        lag_sum = 0;
        % p=1 for Model 1.
        if t-1 >= 1
            lag_sum = lag_sum + b_hat1(2)*irf_temp(t-1);
        end
        irf_temp(t) = b_hat1(1) + lag_sum;
    end
    IRF_model1(:,irep) = irf_temp;
    
    % --- Model 2: y_t = constant + 4 lags of y + eps_t ---
    p2 = 4;
    Y_reg2 = y(p2+1:end);
    Xtmp = [];
    for j = 1:p2
        Xtmp = [Xtmp, y(p2+1 - j : end - j)];
    end
    Etmp = eps_data(p2+1:end);
    X_reg2 = [ones(T-p2,1), Xtmp, Etmp];
    b_hat2 = pinv(X_reg2) * Y_reg2;  % b_hat2 = [alpha; phi1; phi2; phi3; phi4; gamma]
    
    % Compute the IRF for Model 2 inline.
    irf_temp = zeros(H_irf+1,1);
    irf_temp(1) = b_hat2(1) + b_hat2(end);  % initial response: alpha + gamma
    for t = 2:(H_irf+1)
        lag_sum = 0;
        % Use up to 4 lags (or fewer if t-1 < 4)
        for j = 1:min(p2, t-1)
            lag_sum = lag_sum + b_hat2(1+j)*irf_temp(t - j);
        end
        irf_temp(t) = b_hat2(1) + lag_sum;
    end
    IRF_model2(:, irep) = irf_temp;
    
    % --- Model 3: y_t = constant + 12 lags of y + eps_t ---
    p3 = 12;
    Y_reg3 = y(p3+1:end);
    Xtmp = [];
    for j = 1:p3
        Xtmp = [Xtmp, y(p3+1 - j : end - j)];
    end
    Etmp = eps_data(p3+1:end);
    X_reg3 = [ones(T-p3,1), Xtmp, Etmp];
    b_hat3 = pinv(X_reg3) * Y_reg3;  % b_hat3 = [alpha; phi1; ...; phi12; gamma]
    
    % Compute the IRF for Model 3 inline.
    irf_temp = zeros(H_irf+1,1);
    irf_temp(1) = b_hat3(1) + b_hat3(end);
    for t = 2:(H_irf+1)
        lag_sum = 0;
        for j = 1:min(p3, t-1)
            lag_sum = lag_sum + b_hat3(1+j)*irf_temp(t - j);
        end
        irf_temp(t) = b_hat3(1) + lag_sum;
    end
    IRF_model3(:, irep) = irf_temp;
    
end % End of Monte Carlo replications

%% 3. Compute Median IRFs Across Replications
medIRF1 = median(IRF_model1, 2);
medIRF2 = median(IRF_model2, 2);
medIRF3 = median(IRF_model3, 2);

%% 4. Plot the Median IRFs
figure('Name','Q3: Median IRFs from Misspecified Models','Position',[100 100 1200 400]);

subplot(1,3,1)
plot(0:H_irf, medIRF1, '-o','LineWidth',1.5);
title('Model 1: 1 lag of y, current \epsilon_t');
xlabel('Horizon'); ylabel('Response');

subplot(1,3,2)
plot(0:H_irf, medIRF2, '-o','LineWidth',1.5);
title('Model 2: 4 lags of y, current \epsilon_t');
xlabel('Horizon'); ylabel('Response');

subplot(1,3,3)
plot(0:H_irf, medIRF3, '-o','LineWidth',1.5);
title('Model 3: 12 lags of y, current \epsilon_t');
xlabel('Horizon'); ylabel('Response');

% Save the figure if desired.
saveas(gcf, 'IRFs_Q3_Misspecified.png');


%% 4. Estimate Model with y(t-1), Delta M(t) to Delta M(t-6)
% --------------------------------------------------

clearvars -except a b c d beta kappa rho yStar irf_y;
rng(12345,'twister');  % for reproducibility - seed

% For IRFs or references
H_irf = 20;   % horizon for IRF


%% 1. SIMULATION SETTINGS
T      = 500;    % Sample size (after burn-in)
Tburn  = 100;    % Burn-in period length
Nrep   = 500;    % Number of Monte Carlo replications

sigma_e = 0.00066;  % std dev of monetary shock
sigma_n = 0.007;    % std dev of productivity shock

scale = 1e4; % MARINIHILATOR

% We will recover eps and include the contemporaneous eps and 6 lags.
p_eps = 6;  % so total eps terms = 7 (lag 0 through lag6)

% Preallocate a matrix for the IRFs from each replication.
IRF_Q4 = zeros(H_irf+1, Nrep);

rng(12345, 'twister');  % For reproducibility

%% 3. MONTE CARLO SIMULATIONS
for irep = 1:Nrep
    %% 3a. Simulate Data from the True Model
    % True model: 
    %   M_t = rho*M_{t-1} + eps_t,
    %   y_t = a*y_{t-1} + b*(M_{t-1} - M_{t-2}) + c*eps_t.
    M = zeros(Tburn+T, 1);
    y = zeros(Tburn+T, 1);
    
    % Generate eps shocks for the full sample:
    eps_draw = sigma_e * randn(Tburn+T, 1);

    y = y * scale; % MARINIHILATOR

    % Initialization (set M(0)=0 implicitly):
    M(1) = 0;
    M(2) = rho * M(1) + eps_draw(2);
    y(1) = 0;
    % At t = 2, assume M(0)=0 so that (M(1)-M(0)) = M(1)
    y(2) = a * y(1) + b * (M(1) - 0) + c * eps_draw(2);
    
    % Simulate from t = 3 to Tburn+T:
    for t = 3:(Tburn+T)
        M(t) = rho * M(t-1) + eps_draw(t);
        %dM_tminus1 = M(t-1) - M(t-2);
        y(t) = a * y(t-1) + b * M(t-1) + c * eps_draw(t);
    end
    
    % Discard burn-in:
    M = M(Tburn+1:end);
    y = y(Tburn+1:end);
    T_data = length(M);  % should equal T
    
    %% 3b. Recover the Monetary Shock Series: eps_t = M_t - rho*M_{t-1}
    eps_data = zeros(T_data, 1);
    for t = 2:T_data
        eps_data(t) = M(t) - rho * M(t-1);
    end
    eps_data = eps_data * scale; % MARINIHILATOR

    % For t = 1, eps_data(1) remains 0.
    
    %% 3c. Construct the Regression Sample for Q4
    % We want to estimate:
    %   y_t = alpha + phi*y_{t-1} + gamma_0*eps_t + ... + gamma_6*eps_{t-6} + error,
    % so we need to start at a time t where all lagged eps values exist.
    t_start = max(2, p_eps+1);   % here, t_start = 7.
    N_sample = T_data - t_start + 1;
    
    % Dependent variable:
    Y_reg = y(t_start:end);
    
    % Build regressor matrix:
    %   Column 1: constant
    %   Column 2: y_{t-1} (one lag of y)
    %   Columns 3 to 9: eps_t, eps_{t-1}, ..., eps_{t-6}
    X_reg = ones(N_sample, 1);        % constant
    X_reg = [X_reg, y(t_start-1:end-1)]; % lagged y
    
    % Append the shock variables:
    for lag = 0:p_eps
        X_reg = [X_reg, eps_data(t_start - lag:end - lag)];
    end
    % X_reg now has 1 + 1 + 7 = 9 columns.
    
    %% 3d. Estimate the Regression by OLS using pinv
    b_hat = pinv(X_reg) * Y_reg;
    % b_hat structure:
    %   b_hat(1) = alpha, b_hat(2) = phi,
    %   b_hat(3) = gamma_0 (contemporaneous),
    %   b_hat(4) = gamma_1, ..., b_hat(9) = gamma_6.
    
    %% 3e. Compute the Impulse Response from the Estimated Regression
    % We now compute the IRF from the estimated dynamic model.
    % We introduce a one-time shock: let the exogenous shock vector e be such that
    %    e(1) = 1 and e(t) = 0 for t > 1.
    
    alpha_hat = b_hat(1);
    phi_hat   = b_hat(2);
    gamma_hat = b_hat(3:3+p_eps);  % gamma_0 to gamma_6
    
    % Preallocate IRF vector:
    irf_est = zeros(H_irf+1, 1);
    e = zeros(H_irf+1, 1);  % shock sequence
    e(1) = 1;  % a one-time positive shock at period 1
    
    % Period 1:
    irf_est(1) = alpha_hat + gamma_hat(1);  
    % For period 1, lagged y = 0 and only the contemporaneous shock contributes.
    
    % For t = 2,...,H_irf+1, compute the response recursively.
    for t = 2:(H_irf+1)
        shock_effect = 0;
        % Sum the contributions from the shock and its lags:
        for j = 0:p_eps
            if (t - j) >= 1
                shock_effect = shock_effect + gamma_hat(j+1) * e(t - j);
            end
        end
        irf_est(t) = alpha_hat + phi_hat * irf_est(t-1) + shock_effect;
    end
    
    % Store this replication's IRF:
    IRF_Q4(:, irep) = irf_est;
end

%% 4. Compute the Median IRF Across Replications and Plot
medIRF_Q4 = median(IRF_Q4, 2);

figure('Name','Q4: IRF using 1 lag of y and (Contemporary + 6 lags of ε)', ...
       'Position',[100 100 600 400]);
plot(0:H_irf, medIRF_Q4, '-o','LineWidth',1.5);
title('IRF: 1 lag of y + (Contemporary + 6 lags of ε)');
xlabel('Horizon'); ylabel('Impulse Response of y');
grid on;
saveas(gcf, 'IRFs_Q4_MonetaryShockLags.png');