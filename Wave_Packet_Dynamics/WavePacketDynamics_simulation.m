% Wave Packet Dynamics
%
% The evolution of a system of equations (used to approximate the
% Schrodinger equation) have been mapped onto an effective classical system
% in 8 equations
% dq/dt = (1/m)*p --- 3 components
% dp/dt = -del(q)V --- 3 components
% dGamma/dt = (1/m)*eta
% dEta/dt = (9 * reduced_planck_const^2) /(4 * mass^2 * gamma^3) - ...
%                                               partial_V/partial_gamma   
%                                     
% NOTE: V is the potential function
%
% We use the velocity verlet integrator to solve the system in a simulation
% 
% gamma(t) = wave packet width 
% eta(t) = wave packet width momentum
%
% Resources Used:
% 1. http://www.cchem.berkeley.edu/chem195/_n_v_e___verlet_8m.html
% 2. https://people.sc.fsu.edu/~jburkardt/m_src/md/md.m


% environment/material configuration
mass = 1;
reduced_planck_constant = 1;
num_particles = 1;
num_dimensions = 3;
epsilon = 1;
delta_t = 0.001; % integration time

% simulation parameters
simulation_steps = 100000; % the number of integration steps to take
potential_operator_idx = 1; % we'll use this to toggle different equations of V
log_frequency = 10000;


for simulation_step = 1:simulation_steps
    
    % STEP 1 - Solve for Q (position) and P (velocity) using Velocity
    % Verlet
    
    % initialize the positions, velocities and accelerations if starting
    % the simulation
    if (simulation_step == 1)
        % the starting position for the (single) particle is (1,1,1)
        [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = initialize(num_particles, num_dimensions, simulation_steps);
        
        % compute the first forces pPRIME and etaPRIME
        if potential_operator_idx == 1
            pPRIME_acc = compute_force_pPRIME_1(pPRIME_acc, num_dimensions, num_particles, q_pos, simulation_step, epsilon, mass);
            etaPRIME_acc = compute_force_etaPRIME_1(etaPRIME_acc, num_particles, gamma_packet_width, simulation_step, epsilon, mass, reduced_planck_constant);
        end
    else
        % update the positions, velocities and acceleration using the
        % velocity verlet algorithm
        [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = ...
            velocity_verlet(q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc, ...
            mass, delta_t, epsilon, reduced_planck_constant, simulation_step,...
            num_dimensions, num_particles, ...
            potential_operator_idx);
    end
    
    
    % LOG RESULTS
    if mod(simulation_step, log_frequency) == 0
       disp("Updating plot at step: " + simulation_step);
       
       % we'll plot the positions as a continuous line in 3D with
       % annotations on the line at these intervals so that we can
       % discern the trajectory over time
       plot3(reshape(q_pos(1, 1, 1:simulation_step), [1, simulation_step]),...
           reshape(q_pos(2, 1, 1:simulation_step), [1, simulation_step]),...
           reshape(q_pos(3, 1, 1:simulation_step), [1, simulation_step]),...
           '-*', 'MarkerIndices', 1:1000:simulation_step);
    end
    
    % TODO - add some data report   ing here on a set step frequency basis -
    % see resource linked above for example
    
end