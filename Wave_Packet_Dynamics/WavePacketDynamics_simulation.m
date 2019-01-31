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
omega = sqrt((2 * epsilon) / mass);
gamma_0 = 1;
A_0 = (9 * reduced_planck_constant * reduced_planck_constant) / ...
    (8 * mass * mass * gamma_0 * gamma_0 * epsilon);

% simulation parameters
simulation_steps = 10000; % the number of integration steps to take
delta_t = 0.001; % integration time
t_total = simulation_steps * delta_t;
potential_operator_idx = 1; % we'll use this to toggle different equations of V

t = linspace(0, t_total, simulation_steps); % t represents all the integration times in the interval 0 through 10

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
            pPRIME_acc = compute_force_pPRIME_1(pPRIME_acc, num_dimensions, num_particles, q_pos, simulation_step, epsilon);
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
end % end simulation loop


% plot the components of q separately, against the know solutions:
%
% q_1(t) = cos(omega*t), q_2(t) = 0, q_3(t) = 0
% eta(t) = root(A_0 * sin^2(omega * t) + eta_0^2 * cos^2(omega * t))
%
% omega = root(2*epsilon / mass)
% A_0 = 9*reduced_planck_constant^2 / 8 * mass^2 * eta_0^2 * epsilon
%

% log the results
figure(1);
plot(t, reshape(q_pos(1, 1, :), 1, simulation_steps)); % need to reduce the dimensionality
hold on
plot(t, cos(omega*t));
hold off
xlabel('time');
ylabel('q1');
title('q1 - position, component 1');
drawnow

figure(2);
plot(t, reshape(q_pos(2, 1, :), 1, simulation_steps));
xlabel('time');
ylabel('q2');
title('q2 - position, component 2');
drawnow

figure(3);
plot(t, reshape(q_pos(3, 1, :), 1, simulation_steps));
xlabel('time');
ylabel('q3');
title('q3 - position, component 3');
drawnow


% plot gamma
figure(4);
plot(t, gamma_packet_width);
hold on
plot(t, sqrt(...
    (A_0 * sin(omega * (t)).^2) + ...
    (gamma_0 * gamma_0 * cos(omega * t).^2 )));
hold off
xlabel('time');
ylabel('gamma');
title('gamma - packet width');
drawnow