%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MINIMAL TEST OF POSITIVE IRF
% (No productivity shock, just monetary shock, large enough to see an effect)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. CHOOSE SIMPLE, CLEAR PARAMETERS
T        = 5000;       % Big sample so the regression can detect the shock
Tburn    = 100;        % Burn-in
Nrep     = 1;          % Just one replication for clarity
rho      = 0.8;        % AR(1) for money
a        = 0.9;        % y(t) depends heavily on y(t-1)
b        = 0.1;        % moderate effect of delta M_{t-1}
c        = 0.5;        % CLEARLY POSITIVE effect of eps_t
sigma_e  = 0.1;        % a decent-sized shock so it's clearly visible
H_irf    = 20;         % horizon for IRF
p_eps    = 6;          % number of shock lags in the regression

rng(123,'twister');

%% 2. SIMULATE THE TRUE MODEL (NO PRODUCTIVITY SHOCK)
%   M_t = rho*M_{t-1} + eps_t
%   y_t = a*y_{t-1} + b*(M_{t-1}-M_{t-2}) + c*eps_t
M = zeros(Tburn+T,1);
y = zeros(Tburn+T,1);

eps_draw = sigma_e*randn(Tburn+T,1);

% Initialize
M(1) = 0; 
M(2) = rho*M(1) + eps_draw(2);
y(1) = 0;
% t=2
y(2) = a*y(1) + b*(M(1)-0) + c*eps_draw(2);

for t = 3:(Tburn+T)
    M(t) = rho*M(t-1) + eps_draw(t);
    %dM_tminus1 = M(t-1) - M(t-2);
    y(t) = a*y(t-1) + b*M(t-1) + c*eps_draw(t);
end

% Discard burn-in
M = M(Tburn+1:end);
y = y(Tburn+1:end);
T_data = length(M);

%% 3. RECOVER THE MONETARY SHOCK: eps_t = M_t - rho*M_{t-1}
eps_data = zeros(T_data,1);
for t=2:T_data
    eps_data(t) = M(t) - rho*M(t-1);
end
% eps_data(1)=0 by default

%% 4. REGRESSION: y_t = alpha + phi*y_{t-1} + gamma0*eps_t + ... + gamma6*eps_{t-6}
t_start = p_eps+1;   % =7
N_sample = T_data - t_start + 1;

Y_reg = y(t_start:end);

X_reg = ones(N_sample,1);                     % constant
X_reg = [X_reg, y(t_start-1:end-1)];          % 1 lag of y
for lag=0:p_eps
    X_reg = [X_reg, eps_data(t_start-lag:end-lag)];
end
% So X_reg has 1 + 1 + 7 = 9 columns.

b_hat = pinv(X_reg)*Y_reg; 
% b_hat(1)= alpha, b_hat(2)= phi, b_hat(3)= gamma_0, b_hat(4)=gamma_1,...,b_hat(9)= gamma_6

%% 5. COMPUTE THE IRF
%   IRF means a one-time shock e(1)=+1, e(t)=0 for t>1.
alpha_hat = b_hat(1);
phi_hat   = b_hat(2);
gamma_hat = b_hat(3:end);  % gamma_0..gamma_6

irf_est = zeros(H_irf+1,1);
e       = zeros(H_irf+1,1);
e(1)    = 1;   % one-time positive shock

% t=1
irf_est(1) = alpha_hat + gamma_hat(1);

% t=2..(H_irf+1)
for t=2:(H_irf+1)
   shock_effect = 0;
   for j=0:p_eps
       if (t-j)>=1
           shock_effect = shock_effect + gamma_hat(j+1)*e(t-j);
       end
   end
   irf_est(t) = alpha_hat + phi_hat*irf_est(t-1) + shock_effect;
end

%% 6. PLOT THE IRF
figure('Name','Minimal Test: IRF Should be Clearly Positive','Position',[100 100 600 400])
plot(0:H_irf, irf_est,'-o','LineWidth',1.5);
title('IRF: a=0.9, b=0.1, c=0.5, sigma_e=0.1');
xlabel('Horizon'); ylabel('y response');
grid on;