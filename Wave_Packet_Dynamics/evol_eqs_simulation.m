% Evolution equations in terms of a arbitrary potential surface. The system
% described by the following equations:
% q_first_deriv = (1/mass) * p
% p_first_deriv = -del_q * V
% gamma_deriv = (1/mass) * eta
% eta_deriv = (9 * reduced_planck_const^2) / (4 * mass^2 * gamma^3) - partial_V/partial_gamma
%
% NOTE: In our case, the reduced_planck_const = 1
% NOTE: V is the potential function
%
% The point:
% The evolution of the parameters can be modeled by this simple system and
% then plugged into the Guassian equation which being used as a good
% approximation/approach to the Wave Function


% environment/material configuration
mass = 1;
reduced_planck_constant = 1;

% initial values
pos_0 = 0;

dt = 0.001; % integration time

% simulation parameters
simulation_steps = 100000; % the number of integration steps to take

% set up the system of ODES
% dqdt = (1/mass) * INTEGRATE(dpdt);
% dpdt = -1 * gradient(INTEGRATE(dqdt)) * V
% dGammadt = (1/mass) * INTEGRATE(dEtadt)
% dEtadt = 9 * reduced_planck_constant * reduced_planck_constant *
% (1/9*mass*mass*gamma*gamma*gamma) - partialVpartialGamma

for step = 1:simulation_steps
    
end


function dqdt = diff_eq_DQDT(p)
    dqdt = (1/mass) * p;
end

function dpdt = diff_eq_DPDT(del_q,V)
    dpdt = -1 * del_q * V;
end

function dGammadt = diff_eq_DGAMMADT(eta)
    dGammadt = (1/mass) * eta;
end

function dEtadt = diff_eq_DETADT(gamma, pVpGamma)
    dEtadt = 9 * reduced_planck_constant * reduced_planck_constant * ...
        (1 / (9 * mass * mass * gamma * gamma * gamma)) - pVpGamma; 
end

