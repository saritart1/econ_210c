%clear all 

%define the endogeneous variables of the model (no shocks)
var pi, ygap, i, rn, g; %var from variables and you put the variables, inflation, output gap, natural rate of interest and government spending

%define the shocks (exogeneous shocks)
varexo u_g; %varexo from variable exogeneous (variable spending shock

%parameters of our model
parameters gamma, phi, lambda, beta, sg, rho_g, phi_pi, sigma_g, omega, phi_ygap;


%-------------------------------------------------------
% 2. Calibration %tell Dynare what values they should take to our model
%-------------------------------------------------------

load paramfile; %call other file that has the values of the parameters but you can also just put the value (this is iin the main document Paula sent us)
%if we dont want to have a separate file we can put it directly here for example; gamma = .1; etc

set_param_value('gamma',gamma);
set_param_value('phi',phi);
set_param_value('lambda',lambda);
set_param_value('beta',beta);
set_param_value('sg',sg);
set_param_value('rho_g',rho_g);
set_param_value('phi_pi',phi_pi);
set_param_value('sigma_g',sigma_g);
set_param_value('omega',omega);
set_param_value('phi_ygap',phi_ygap);

%-------------------------------------------------------
% 3. Model
%-------------------------------------------------------
model(linear); 
%this is a program so you use end at the end. This is using the already log linearized version of the model, but Dynare usually log linearizes it, so to stop this from happening, we say model(linear) to tell the model that is is already linearlized
    
    % (1) Phillips Curve
    pi = lambda * omega * ygap + beta * pi(+1) ;

    % (2) Dynamic IS
    ygap = ygap(+1) - ( (1-sg) / gamma ) * (i - pi(+1) - rn) ;

    % (3) Natural interest rate
    rn = - (gamma * sg * phi * ( rho_g - 1) * g) / ( phi * (1-sg) + gamma ) ;
    
    % (4) Taylor Rule
    %i = phi_pi * pi ;
    i = phi_pi * pi + phi_ygap * ygap;

    % (5) Exogenous process for government spending
    g = rho_g * g(-1) + u_g;

  
end;

%-------------------------------------------------------
% 4. Steady State file
%------------------------------------------------------

%steady_state_model; 



%end;

%-------------------------------------------------------
% 5. Shocks
%-------------------------------------------------------

shocks;

    
    var u_g = sigma_g^2; 
%what is the process for the exogeneous variable

end;

%-------------------------------------------------------
% 6. Computation
%-------------------------------------------------------

steady;

check;

stoch_simul(irf=20); %solve and simulate the model for 20 horizons
save nkm_3eq.mat


