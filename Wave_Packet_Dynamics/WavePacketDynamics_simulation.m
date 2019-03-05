%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Wave Packet Dynamics
%
% Simulate the dynamics of a single electron orbiting a nucleus of
% arbitrary atomic number Z
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
% NOTE: Atomic units are assumed and used throughout the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% environment/material configuration
mass = 1;
reduced_planck_constant = 1;
num_particles = 1;
num_dimensions = 3;
epsilon = 1;
% -----------------------------------

% simulation parameters
simulation_steps = 100; % the number of integration steps to take
delta_t = .1; % integration time change
t_total = simulation_steps * delta_t;
t = linspace(0, t_total, simulation_steps); % t represents all the integration times in the interval 0 through t_total
% ---------------------------------------------

% values used for Coulomb interaction operator
Z = 1; % atomic number
e = 1;
A = 1;
x_0 = .5;
gamma_0 = 1;
%A = (9 * reduced_planck_constant * reduced_planck_constant) / ...
 %   (4 * mass * mass * Z * e * e); % the transformation introduced in eqs (28), (29) allow us to represent all
                                    % tunable parameters in the simulation
                                    % with this new variable A
% toggles for which plots we want to see (1-yes, 0-no)
% 1. q over time
% 2. gamma over time
% 3. q vs. gamma
% 4. isocontours of V
% 5. surface of V
% 6. pseudocolor plot of V
visualization_toggles = [0, 0, 0, 1, 1, 1];
% ---------------------------------------------


% values used for quadratic well operator
gamma_0_well = 1;
A_0 = (9 * reduced_planck_constant * reduced_planck_constant) / ...
    (8 * mass * mass * gamma_0_well * gamma_0_well * epsilon);
omega = sqrt((2 * epsilon) / mass);
% ----------------------------------------

% The following index designates which potential operator the simulation
% will apply to the electron:
% Quadratic Well = 1
% Coulomb Interaction = 2
potential_operator_idx = 2;

for simulation_step = 1:simulation_steps
    % STEP 1 - Solve for Q (position) and P (velocity) using Velocity
    % Verlet
    
    % initialize the positions, velocities and accelerations if starting
    % the simulation
    if (simulation_step == 1)
        if potential_operator_idx == 1
            % initialize the first values according to the starting values
            % designated for the quadratic well operator
            [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = initialize_quadraticWell(num_particles, num_dimensions, simulation_steps);
        
            % compute the first forces pPRIME and etaPRIME given the
            % quadratic well operator
            pPRIME_acc = compute_force_pPRIME_quadraticWell(pPRIME_acc, num_dimensions, num_particles, q_pos, simulation_step, epsilon);
            etaPRIME_acc = compute_force_etaPRIME_quadraticWell(etaPRIME_acc, num_particles, gamma_packet_width, simulation_step, epsilon, mass, reduced_planck_constant);
        elseif potential_operator_idx == 2
            [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = initialize_coulombInteraction(num_particles, simulation_steps, A, x_0, gamma_0);
            
            % compute the first forces pPRIME and etaPRIME given the
            % Coulomb interaction operator
            pPRIME_acc = compute_force_pPRIME_coulombInteraction(pPRIME_acc, num_particles, Z, q_pos, gamma_packet_width, simulation_step);
            etaPRIME_acc = compute_force_etaPRIME_coulombInteraction(etaPRIME_acc, num_particles, A, Z, reduced_planck_constant, mass, gamma_packet_width, q_pos, simulation_step);
        else
            disp('No potential operator index given, or incorrectly set. Returning');
            return;
        end
    else
        % update the positions, velocities and acceleration using the
        % velocity verlet algorithm
        [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = ...
            velocity_verlet(q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc, ...
            mass, delta_t, epsilon, reduced_planck_constant, A, Z, simulation_step,...
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
if potential_operator_idx == 1
    figure(1);
    plot(t, reshape(q_pos(1, 1, :), 1, simulation_steps)); % need to reduce the dimensionality
    hold on
    plot(t, cos(omega*t));
    hold off
    xlabel('time');
    ylabel('q1');
    title('Quadratic Well Operator: q1 - position, component 1');
    drawnow

    figure(2);
    plot(t, reshape(q_pos(2, 1, :), 1, simulation_steps));
    xlabel('time');
    ylabel('q2');
    title('Quadratic Well Operator: q2 - position, component 2');
    drawnow

    figure(3);
    plot(t, reshape(q_pos(3, 1, :), 1, simulation_steps));
    xlabel('time');
    ylabel('q3');
    title('Quadratic Well Operator: q3 - position, component 3');
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
    title('Quadratic Well Operator: gamma - packet width');
    drawnow
elseif potential_operator_idx == 2
    visualize_coulomb(t, q_pos, gamma_packet_width, Z, A, e, visualization_toggles);
end
