%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code does the programming part of the Problem Set from
% Nakamura and Steinsson (2018)
% IMPORTANT:
% For point 1 I used the linear approximation of the output but also for 
% the inflation, thus I assumed some values to get a,b,c,d for the next
% exercises:
% e = 0.3
% f = 0.2
% g = 0.1
% h = 0.05
% (I will double check with Paula if this is correct or I need to properly
% find them, but I still don't know how)
% I should also ask to Paula if I should assume that all the variables in
% SS = 0 (y_t, m_t and pi_t) and also y_star? 
% For Paula: can I assume the Phillips curve to be inflation = kappa(output
% gap) or should I derive the expectation of inflation?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 2. TRUE IMPULSE RESPONSE (THEORETICAL) 
% -------------------------------------------------------------------------

% Horizon for IRF
H = 20;  % Number of periods to compute the IRF

% Parameters from Nakamura and Steinsson (2018)
beta  = 0.99;    % Discount factor
kappa = 0.13;    % Phillips Curve slope
rho   = 0.8;     % AR(1) coefficient for M_t
yStar = 0;       % Natural output

% Inflation coefficients (assumed reasonable values)
e = 0.3;         % Coefficient on y_{t-1}
f = 0.2;         % Coefficient on Delta M_{t-1}
g = 0.1;         % Coefficient on epsilon_t
h = 0.05;        % Coefficient on eta_t

% Compute denominator for a,b,c,d
denom = beta * e + kappa;  % = 0.99 * 0.3 + 0.13 = 0.427

% Compute reduced-form coefficients for y_t
a = e / denom;                         % = 0.7026
b = (f - beta * f * rho) / denom;      % = 0.0974
c = (g - beta * f) / denom;            % = -0.2295
d = (h - kappa) / denom;               % = -0.1874

% Display coefficients
fprintf('Reduced-form coefficients:\n');
fprintf('a = %.4f\n', a);
fprintf('b = %.4f\n', b);
fprintf('c = %.4f\n', c);
fprintf('d = %.4f\n', d);

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
    deltaMnow = Mnow - Mprev;
    
    ynow = a * yprev + b * deltaMprev + c * eps_now + d * eta_now;
    pinow = kappa * (ynow - yStar);
    
    % Store IRFs
    irf_y(t+1)  = ynow;
    irf_dm(t+1) = deltaMnow;
    irf_pi(t+1) = pinow;
    
    % Update
    Mprev = Mnow;
    yprev = ynow;
    deltaMprev = deltaMnow;
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

disp('2) True IRFs computed and plotted using your derived coefficients.');
% Save figure
saveas(gcf, 'IRFs_Q2_TrueModel.png');


%% 3. SIMULATE DATA AND ESTIMATE MISSPECIFIED MODELS
% --------------------------------------------------

clear; clc; close all;

% Set parameters from your model
beta  = 0.99;
kappa = 0.13;
rho   = 0.8;
yStar = 0;

e = 0.3;
f = 0.2;
g = 0.1;
h = 0.05;

denom = beta * e + kappa;

% Coefficients for true model
a = e / denom;
b = (f - beta * f * rho) / denom;
c = (g - beta * f) / denom;
d = (h - kappa) / denom;

fprintf('True Model Coefficients:\n a=%.4f, b=%.4f, c=%.4f, d=%.4f\n', a, b, c, d);

% Simulation settings
T = 500;       % Time series length
reps = 500;    % Number of replications
H_irf = 20;    % IRF horizon

sigEps = 0.00066;   % Std dev of monetary shock
sigEta = 0.007;     % Std dev of productivity shock

rng(123);  % For reproducibility - seed

% Define lag lengths for misspecified models
lagList = [1, 4, 12]; %for Model 1, Model 2 and Model 3

% Storage for IRFs
IRFs_all = zeros(H_irf+1, reps, length(lagList));

% Loop over replications
for rr = 1:reps
    
    % Preallocate
    y_data = zeros(T,1);
    M_data = zeros(T,1);
    
    % Generate shocks
    eps_vec = sigEps * randn(T,1); % monetary
    eta_vec = sigEta * randn(T,1); % productivity
    
    % Initial conditions
    y_data(1) = 0;
    M_data(1) = 0;
    
    % Simulate the true model % from question 1
    for t=2:T
        M_data(t) = rho*M_data(t-1) + eps_vec(t);
        deltaM_t = M_data(t) - M_data(t-1);
        
        y_data(t) = a*y_data(t-1) + b*deltaM_t + c*eps_vec(t) + d*eta_vec(t);
    end
    
    % Estimate misspecified models
    for m = 1:length(lagList)
        
        p_lags = lagList(m);
        
        % Build lag matrix
        Y_lag = lagmatrix(y_data, 1:p_lags);
        Y_dep = y_data(p_lags+1:end);
        X = [ones(T-p_lags,1), Y_lag(p_lags+1:end,:)]; % Add constant
        
        % Estimate AR(p) via OLS
        coeffs = (X' * X) \ (X' * Y_dep);
        
        % Build companion form matrix for IRF
        A = zeros(p_lags, p_lags);
        A(1,:) = coeffs(2:end)';
        A(2:end,1:end-1) = eye(p_lags-1);
        
        % Shock scaling: residual std
        res = Y_dep - X * coeffs;
        sigma_eps = std(res);
        
        % Compute IRF for y_t
        irf = zeros(H_irf+1,1);
        state = zeros(p_lags,1);
        
        % Impact response
        irf(1) = sigma_eps;
        state(1) = sigma_eps;
        
        % Propagate IRF
        for h = 1:H_irf
            state = A * state;
            irf(h+1) = state(1);
        end
        
        IRFs_all(:, rr, m) = irf;
    end
end

figure('Name','Median IRFs - Misspecified Models','Position',[100 100 1000 400]);

for m = 1:length(lagList)
    
    subplot(1,length(lagList),m);
    irf_median = median(IRFs_all(:,:,m),2); % Median across replications
    
    plot(0:H_irf, irf_median, '-o','LineWidth',1.5);
    title(sprintf('Misspecified AR(%d)', lagList(m)));
    xlabel('Horizon'); ylabel('Response of y');
    grid on;
end

disp('3) Median IRFs for misspecified models plotted.');
% Save the figure
saveas(gcf, 'IRFs_Q3_MisspecifiedModels.png');

%% QUESTION 4: Estimate Model with y(t-1), Delta M(t) to Delta M(t-6)
% ------------------------------------------------------------------------

clear; clc; close all;

% -- Suppose you already have these from earlier assignments:
a     = 0.5;     % True model param (from Q1)
b     = 0.0974;     % True model param
c     = -0.2295;    % True model param
d     = -0.1874;    % True model param
rho   = 0.8;      % AR coeff for M_t
sigEps= 0.00066;  % stdev of monetary shock
sigEta= 0.007;    % stdev of productivity shock
T     = 500;      % Sample size
reps  = 500;      % Number of simulations
H_irf = 20;       % IRF horizon
yStar = 0;        % Potential output

rng(123); % Reproducibility

% We'll store IRFs across replications
IRF_money_all = zeros(H_irf+1, reps);

for rr = 1:reps
    
    %% (1) Simulate data from true model
    y_data = zeros(T,1);
    M_data = zeros(T,1);

    eps_vec = sigEps * randn(T,1);  % monetary shocks
    eta_vec = sigEta * randn(T,1);  % productivity shocks

    % Initialize
    y_data(1) = 0;
    M_data(1) = 0;

    for t = 2:T
        M_data(t) = rho * M_data(t-1) + eps_vec(t);
        deltaM_t  = M_data(t) - M_data(t-1);
        
        % True output eq:
        y_data(t) = a * y_data(t-1) + b * deltaM_t + c * eps_vec(t) + d * eta_vec(t);
    end
    
    %% (2) Build the regression: y_t on [y_{t-1}, DeltaM_t, ..., DeltaM_{t-6}]
    deltaM_data = [0; diff(M_data)];  % money growth

    % We'll have 1 lag of y, plus current + 6 lags of DeltaM => 7 DeltaM terms total
    lags_output = 1;  
    lags_money  = 6;  

    Y_lag = lagmatrix(y_data, 1:lags_output);   
    M_lag = lagmatrix(deltaM_data, 0:lags_money); 

    % Drop initial NaN rows
    startIdx = lags_money + 1 + 1;  
    Y_dep    = y_data(startIdx:end);
    X        = [ ones(length(Y_dep),1), ...
                 Y_lag(startIdx:end,:), ...
                 M_lag(startIdx:end,:) ];

    coeffs = (X' * X) \ (X' * Y_dep);
    
    % Residual stdev
    res       = Y_dep - X * coeffs;
    sigma_eps = std(res);

    %% (3) Build companion matrix for IRF
    % State vector: [ y_{t-1}, DeltaM_t, DeltaM_{t-1}, ..., DeltaM_{t-6} ]
    state_size = 8;  
    A = zeros(state_size, state_size);

    % Interpreting 'coeffs':
    %  - coeffs(1) is constant
    %  - coeffs(2) is y_{t-1}'s coeff
    %  - coeffs(3:9) are the 7 money-growth terms (DeltaM_t..DeltaM_{t-6})

    % Row 1 of A has the regression eq (minus constant)
    A(1,1)      = coeffs(2);         % y_{t-1}
    A(1,2:8)    = coeffs(3:9)';      % 7 money-lag coeffs (transpose => row vector)

    % Shift money lags down by 1 each period
    % We'll fill the sub-block (3..8, 2..7) with eye(6)
    A(3:8, 2:7) = eye(6);

    %% (4) Compute IRF: shock is +1 std dev in current DeltaM_t
    irf = zeros(H_irf+1,1);
    
    % Impact state (time 0):
    state = zeros(state_size,1);
    state(2) = sigma_eps;   % put shock in current DeltaM slot

    % Immediate impact on y
    irf(1) = A(1,:) * state; 

    % Evolve
    for h = 1:H_irf
        state = A * state;
        irf(h+1) = state(1);
    end
    
    IRF_money_all(:, rr) = irf;
end

%% (5) Plot median IRF
figure('Name','Q4: IRF of y to Money Growth Shock','Position',[100 100 700 400]);
irf_median = median(IRF_money_all, 2);

plot(0:H_irf, irf_median, '-o','LineWidth',1.5);
title('Median IRF of y to Monetary Shock (Q4 Model)');
xlabel('Horizon'); ylabel('Response of y');
grid on;

% Save figure if desired
saveas(gcf, 'IRF_Q4_Misspecified.png');

%% QUESTION 5: Jordà Local Projection IRF
% ------------------------------------------------------------------------

clearvars -except a b c d rho sigEps sigEta reps T H_irf yStar % Keep parameters

rng(123); % For reproducibility

H_lp = 20;  % Horizon for local projection

IRF_lp_all = zeros(H_lp+1, reps);

for rr = 1:reps
    
    % (1) Simulate data from true model
    y_data = zeros(T,1);
    M_data = zeros(T,1);

    eps_vec = sigEps * randn(T,1); % Monetary shocks (this is your Shock_t)
    eta_vec = sigEta * randn(T,1); % Productivity shocks

    y_data(1) = 0;
    M_data(1) = 0;

    for t = 2:T
        M_data(t) = rho * M_data(t-1) + eps_vec(t);
        deltaM_t  = M_data(t) - M_data(t-1);
        
        y_data(t) = a * y_data(t-1) + b * deltaM_t + c * eps_vec(t) + d * eta_vec(t);
    end
    
    % (2) For each horizon h, run the local projection regression
    beta_h = zeros(H_lp+1, 1); % Store IRF for this replication
    
    for h = 0:H_lp
        idx = 1:(T-h);
        Y_h = y_data(idx + h);
        Shock = eps_vec(idx);  % This is your "identified" shock
        
        % Regress y_{t+h} on shock_t + constant
        X = [ones(length(Shock),1), Shock];
        coeff = (X' * X) \ (X' * Y_h);
        
        beta_h(h+1) = coeff(2); % Coefficient on shock
    end
    
    IRF_lp_all(:, rr) = beta_h;
end

%% Plot median IRF

figure('Name','Local Projection IRF','Position',[100 100 700 400]);
irf_median_lp = median(IRF_lp_all, 2);

plot(0:H_lp, irf_median_lp, '-o','LineWidth',1.5);
title('Median IRF of y to Monetary Shock (Local Projection)');
xlabel('Horizon'); ylabel('Response of y');
grid on;

saveas(gcf, 'IRF_LocalProjection.png');
