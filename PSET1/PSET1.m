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