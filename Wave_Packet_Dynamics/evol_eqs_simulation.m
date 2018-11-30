% Wave Packet Dynamics
%
% The evolution of a system of equations (used to approximate the
% Schrodinger equation) have been mapped onto an effective classical system
%
%
% Evolution equations in terms of a arbitrary potential surface. The system
% described by the following equations:
% dqdt = q_first_deriv = (1/mass) * p
% dpdt = p_first_deriv = -del_q * V
% dGammadt = gamma_first_deriv = (1/mass) * eta
% dEtadt = eta_first_deriv = (9 * reduced_planck_const^2) /...
%                               (4 * mass^2 * gamma^3) - partial_V/partial_gamma
%
% NOTE: In our case, the reduced_planck_const = 1
% NOTE: V is the potential function
%
% How do we relate this system of ODEs to position, velocity, acceleration
% and force so that we can use the Velocity Verlet algorithm to solve the
% system in simulation?
% q(t) = position
% p(t) = velocity
%
% dqdt = velocity
% dpdt = acceleration
% V = force
% 
% gamma(t) = wave packet width 
% eta(t) = wave packet width momentum
%
% The point:
% The evolution of the parameters can be modeled by this simple system and
% then plugged into the Guassian equation which being used as a good
% approximation/approach to the Wave Function
%
% Resources Used:
% 1. http://www.cchem.berkeley.edu/chem195/_n_v_e___verlet_8m.html
% 2. https://people.sc.fsu.edu/~jburkardt/m_src/md/md.m


% environment/material configuration
mass = 1;
reduced_planck_constant = 1;
num_particles = 1;
num_dimensions = 1;

% initial values
pos_0 = 0;
pos_OLD = -0.001; % the position in the previous time step

dt = 0.001; % integration time

% simulation parameters
simulation_steps = 100000; % the number of integration steps to take


for step = 0:simulation_steps
    
    % STEP 1 - Solve for Q (position) and P (velocity) using Velocity
    % Verley
    
    % initialize the positions, velocities and accelerations if starting
    % the simulation
    if (step == 0)
        [positions, velocities, accelerations] = initialize(num_particles, num_dimensions);
    else
        % update the positions, velocities and acceleration using the
        % velocity verlet algorithm
        [positions, velocities, accelerations] = velocity_verlet(num_particles, num_dimensions,...
            positions, velocities, forces, accelerations, mass, dt);
    end
    
    % compute the forces and energies at every step
    [force, potential_energy, kinetic_enery] = compute_forces(num_particles, num_dimensions,...
        positions, velocities, mass); % this call may be updated to reflect the ability to specific different potential operators
    
    
    % STEP 2 - Solve for GAMMA and ETA
    
    % STEP 3 - Solve GAUSSIAN TEST FUNCTION
    
    % LOG RESULTS
    
    % TODO - add some data reporting here on a set step frequency basis -
    % see resource linked above for example
    
end


% NOT SURE HOW/IF THESE FUNCTIONS FACTOR IN HERE

% The velocity equation for the Velocity Verlet algorithm
function dqdt = diff_eq_DQDT(p)
    dqdt = (1/mass) * p;
end

% The acceleration equation for the Velocity Verlet algorithm
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

