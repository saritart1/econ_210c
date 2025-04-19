%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code does the programming part of the Problem Set from
% Nakamura and Steinsson (2018)
% Sara Restrepo Tamayo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;  close all;  clc;

%-----------------
% Parameters
%-----------------

%Define the seed:
rng(18042025,"twister");

%Calibration
beta = 0.99;
kappa = 0.13;
rho = 0.8;
sigma_M = 0.00066;
sigma_eta = 0.007;

%Coefficients: y = a*y(-1) + b*ΔM(-1) + c*eps + d*eta
a = 0;
b = (1-beta*rho)*rho/(1+kappa-beta*rho);
c = (1-beta*rho)/(1+kappa-beta*rho);
d = kappa/(1+kappa);

%------------------------
%% QUESTION 1.3:
%------------------------

% Simulation inputs
T = 600; %Burn first 100 simulations
N_simu = 500;

 %IRF storage matrices: 
 t_irf=21;
 IRF_1lag = zeros(t_irf, N_simu); % y ~ y(-1) + eps
 IRF_4lag = zeros(t_irf, N_simu); % y ~ y(-1) + y(-2) + y(-3) + y(-4) + eps
 IRF_12lag = zeros(t_irf, N_simu); % y ~ y(-1) ... + y(-12) + eps


%------------------------
% Loop over simulations
%-----------------------
for s = 1:N_simu

    % 1. Simulate model:
    %-----------------------
    % Generate shocks
    eta = sigma_eta * randn(T+1, 1); % productivity shocks
    eps = sigma_M * randn(T+1, 1); % monetary shocks

    % Simulate ΔM, y
    %Initiate matrices:
    deltaM = zeros(T+1,1);
    y = zeros(T+1,1);

    for t = 2:T+1
        deltaM(t) = rho * deltaM(t-1) + eps(t);
        y(t) = a * y(t-1) + b * deltaM(t-1) + c * eps(t) + d * eta(t);
    end

    % Save simulated data: burn 100 obs
    Y = y(102:end);
    E = eps(102:end); 

    % 2. Run misspecifications:
    %--------------------------
    % IRF: simulate y after ε_0 = 1 std dev
    eps_irf = zeros(t_irf,1); 
    eps_irf(1) = 1*sigma_M;

    %% Model 1: y ~ y(-1) + eps
    %Estimate AR(1)
    X1 = [ones(T-100,1), y(101:end-1), E];
    b1 = pinv(X1)*Y;
    
    %IRF AR(1)
    y1 = zeros(t_irf,1);
    y1(1) = b1(1) + b1(end)*eps_irf(1); % response in t=1
    for t = 2:t_irf
        y1(t) = b1(1) + b1(2) * y1(t-1) + b1(3) * eps_irf(t);
    end
    IRF_1lag(:, s) = y1;

    %% Model 2: y ~ y(-1) + y(-2) + y(-3) + y(-4) + eps
    %Estimate AR(4)
    Y4lags = lagmatrix(y(102:end), 1:4);
    X4 = [ones(T-104,1), Y4lags(5:end,:), E(5:end)];
    Y4 = Y(5:end);
    b4 = pinv(X4)*Y4;

    %IRF AR(4)
    y4 = zeros(t_irf,1);
    y4(1) = b4(1) + b4(end)*eps_irf(1); % response in t=1
    for t = 2:t_irf
        lag_sum = 0;
        % Use up to 4 lags (or fewer if t-1 < 4)
        for j = 1:min(4, t-1)
            lag_sum = lag_sum + b4(1+j)*y4(t-j);
        end
        y4(t) = b4(1) + lag_sum + b4(end) * eps_irf(t);
    end
    IRF_4lag(:, s) = y4;

    %% Model 3:  y ~ y(-1) ... + y(-12) + eps
    %Estimate AR(12)
    Y12lags = lagmatrix(y(102:end), 1:12);
    X12 = [ones(T-112,1), Y12lags(13:end,:), E(13:end)];
    Y12 = Y(13:end);
    b12 = pinv(X12) * Y12;

    %IRF AR(12)
    y12 = zeros(t_irf,1);
    y12(1) = b12(1) + b12(end)*eps_irf(1); % response in t=1
    for t = 2:t_irf
        lag_sum = 0;
        % Use up to 12 lags (or fewer if t-1 < 12)
        for j = 1:min(12, t-1)
            lag_sum = lag_sum + b12(1+j)*y12(t-j);
        end
        y12(t) = b12(1) + lag_sum + b12(end) * eps_irf(t);
    end

    IRF_12lag(:, s) = y12;
end

%------------------------
% Save median IRF:
%-----------------------
median_1 = median(IRF_1lag, 2);
median_4 = median(IRF_4lag, 2);
median_12 = median(IRF_12lag, 2);


%------------------------
% True IRF:
%-----------------------
true_irf = zeros(t_irf,1);
deltaM = zeros(t_irf,1);
eps_irf = zeros(t_irf,1); eps_irf(1) = 1*sigma_M;

for t = 1:t_irf
    if t == 1
        deltaM(t) = eps_irf(t);
        true_irf(t) = c * eps_irf(t);
    else
        deltaM(t) = rho * deltaM(t-1) + eps_irf(t);
        true_irf(t) = a * true_irf(t-1) + b * deltaM(t-1) + c * eps_irf(t);
    end
end

f = figure('Visible','off'); hold on;
plot(0:t_irf-1, true_irf, 'k--', 'LineWidth', 2);
plot(0:t_irf-1, median_1,'Color', [0.7, 0, 0.2], 'LineWidth', 2);
plot(0:t_irf-1, median_4, 'Color',[0.7, 0.5, 0.85], 'LineWidth', 2);
plot(0:t_irf-1, median_12, 'Color', [0, 0.5, 0.5], 'LineWidth', 2);
legend('True IRF','1 lag', '4 lags', '12 lags','Location','north');
xlabel('Time','FontSize', 12);
ylabel('$\mathbf{IRF} \; (y_{t+h})$', 'Interpreter', 'latex','FontSize', 12);
title('Median IRF from Misspecified Models','FontSize',14);

% Save figure
saveas(f, 'ps1_3_IRF_misspec_y.png');
close(f);
%------------------------
%% Question 1.4:
%------------------------

 %IRF storage matrices: 
 t_irf=21;
 IRF_1_6lag = zeros(t_irf, N_simu); % y ~ y(-1) + eps + eps(-1)+...+eps(-6)

%------------------------
% Loop over simulations
%-----------------------
for s = 1:N_simu

    % 1. Simulate model:
    %-----------------------
    % Generate shocks
    eta = sigma_eta * randn(T+7, 1); % productivity shocks
    eps = sigma_M * randn(T+7, 1); % monetary shocks

    % Simulate ΔM, y
    %Initiate matrices:
    deltaM = zeros(T+7,1);
    y = zeros(T+7,1);

    for t = 2:T+7
        deltaM(t) = rho * deltaM(t-1) + eps(t);
        y(t) = a * y(t-1) + b * deltaM(t-1) + c * eps(t) + d * eta(t);
    end

    % Save simulated data: burn 100 obs
    Y = y(108:end);
    E = eps(108:end); 

    % 2. Run misspecifications:
    %--------------------------
    % IRF: simulate y after ε_0 = 1 std dev
    eps_irf = zeros(t_irf,1); 
    eps_irf(1) = 1*sigma_M;

    % ε_t to ε_{t-6}
    E_lags = zeros(T-100, 7);
    for j = 0:6
        E_lags(:, j+1) = eps(108-j:end-j);  % eps(t-j)
    end

    %% Model: y(-1) + eps  + eps(-1)+...+eps(-6)
    %Estimate coefficients
    X1 = [ones(T-100,1), y(107:end-1),E_lags];
    b1 = pinv(X1)*Y;
    
    %IRF 
    y1 = zeros(t_irf,1);
    y1(1) = b1(1) + b1(3)*eps_irf(1); % response in t=1

    for t = 2:t_irf
    y1(t) = b1(1) + b1(2) * y1(t-1); 
        for j = 0:min(6, t-1)
            y1(t) = y1(t) + b1(3+j) * eps_irf(t-j);
        end
    end
    IRF_1_6lag(:, s) = y1;
end
%------------------------
% Save median IRF:
%-----------------------
median_1_6 = median(IRF_1_6lag, 2);

f = figure('Visible','off'); hold on;
plot(0:t_irf-1, true_irf, 'k--', 'LineWidth', 2);
plot(0:t_irf-1, median_1_6,'Color', [0.29, 0.33, 0.13], 'LineWidth', 2);
legend({'True IRF', '$y_{t-1} + \sum_{s = t - 6}^{t} \Delta M_s+ \varepsilon_t$'}, ...
       'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast');
xlabel('Time','FontSize', 12);
ylabel('$\mathbf{IRF} \; (y_{t+h})$', 'Interpreter', 'latex', 'FontSize', 12);
title('$\mathrm{Misspecified\ Model:}\ y_{t-1} + \sum_{s = t - 6}^{t} \Delta M_s+ \varepsilon_t$', ...
      'Interpreter', 'latex', 'FontSize', 14);
% Save figure
saveas(f, 'ps1_4_IRF_misspec_y.png');
close(f);

%------------------------
%% Question 1.5:
%------------------------
 %IRF storage matrices: 
 H=20;
 IRF_jorda = zeros(H+1, N_simu); 
 IRF_jorda_control = zeros(H+1, N_simu); %control: output lag

for s = 1:N_simu
    % 1. Simulate model:
    %-----------------------
    % Generate shocks
    eta = sigma_eta * randn(T + H + 1, 1);
    eps = sigma_M * randn(T + H + 1, 1);
   
    % Simulate ΔM, y
    %Initiate matrices:
    deltaM = zeros(T + H + 1, 1);
    y = zeros(T + H + 1, 1);

    for t = 2:T + H + 1
        deltaM(t) = rho * deltaM(t-1) + eps(t);
        y(t) = a * y(t-1) + b * deltaM(t-1) + c * eps(t) + d * eta(t);
    end

    % 2. Local Projection: estimation at each horizon h
    %-------------------------------------------------
    %WITH CONTROL
    for h = 0:(H+1)
        Y_h = y((1+h):(T+h));                % y_{t+h}
        eps_t = eps(2:T+1);                  % ε_t
        y_lag = y(1:T);                      % y_{t-1}
        X = [ones(T,1), eps_t, y_lag];       % Add constant + controls
        b_lp = pinv(X)*Y_h;
        IRF_jorda_control(h+1, s) = b_lp(2)*sigma_M;    % β_h on ε_t
    end

    %WITHOUT CONTROL
    for h = 0:(H+1)
        Y_h = y((1+h):(T+h));                % y_{t+h}
        eps_t = eps(2:T+1);                  % ε_t
        y_lag = y(1:T);                      % y_{t-1}
        X = [ones(T,1), eps_t];       % Add constant
        b_lp = pinv(X)*Y_h;
        IRF_jorda(h+1, s) = b_lp(2)*sigma_M;    % β_h on ε_t
    end
end

% Median IRF
median_LP_control = median(IRF_jorda_control, 2);
median_LP = median(IRF_jorda, 2);

% Plot comparison with true IRF
f = figure('Visible','off'); hold on;
plot(0:H, true_irf(1:H+1), 'k--', 'LineWidth', 3);
plot(0:H, median_LP_control(2:end), 'Color', [0.7, 0, 0.2], 'LineWidth', 2);
plot(0:H, median_LP(2:end), 'Color', [0, 0.5, 0.5], 'LineWidth', 2);
legend({'True IRF', 'LP with control:$y_{t-1}$', 'LP w/o control'}, 'Interpreter', 'latex','Location', 'NorthEast');
xlabel('Time','FontSize', 12);
ylabel('$\mathbf{IRF} \; (y_{t+h})$', 'Interpreter', 'latex', 'FontSize', 12);
title('Jordà Local Projection vs. True IRF','FontSize', 14);

% Save figure
saveas(f, 'ps1_5_Jorda_y.png');
close(f);