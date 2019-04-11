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
simulation_steps = 1000; % the number of integration steps to take
delta_t = .1; % integration time change
t_total = simulation_steps * delta_t;
t = linspace(0, t_total, simulation_steps); % t represents all the integration times in the interval 0 through t_total
% ---------------------------------------------

% values used for quadratic well operator
gamma_0_well = 1;
A_0_well = (9 * reduced_planck_constant * reduced_planck_constant) / ...
    (8 * mass * mass * gamma_0_well * gamma_0_well * epsilon);
omega_well = sqrt((2 * epsilon) / mass);
% ----------------------------------------

% values used for Coulomb interaction operator
Z_CI = 1; % atomic number
e_CI = 1;
A_CI = 1;
x_0_CI = 1;
gamma_0_CI = 1;
%A = (9 * reduced_planck_constant * reduced_planck_constant) / ...
 %   (4 * mass * mass * Z * e * e); % the transformation introduced in eqs (28), (29) allow us to represent all
                                    % tunable parameters in the simulation
                                    % with this new variable A
                                    

% values used for screened Coulomb interaction
Z_SCI = 1;
e_SCI = 1;
q_0_SCI = zeros(num_dimensions, num_particles);
q_0_SCI(1,1) = 1;
q_0_SCI(2, 1) = .3;
q_0_SCI(3, 1) = .1;

gamma_0_SCI = 1.3;
lambda_SCI = 10;
% ---------------------------------------------

% values used for square well potential
q_0_SW = -.2;
gamma_0_SW = .5;
V0_SW = 1;
a_SW = 3;
% ---------------------------------------------

% values used for Morse potential (MP). 1D - from ChemPhysLetter (CPL)
q_0_MP_CPL = 5;
gamma_0_MP_CPL = 1;
D_MP_CPL = 50; % well depth
a_MP_CPL = .15; % width of the potential
b_MP_CPL = .05;
% ---------------------------------------------
                                    
% toggles for which plots we want to see (1-yes, 0-no)
% 1. q over time
% 2. gamma over time
% 3. q vs. gamma
% 4. isocontours of V
% 5. surface of V
% 6. pseudocolor plot of V
visualization_toggles = [1, 1, 1, 0, 1, 0];
% ---------------------------------------------

% The following index designates which potential operator the simulation
% will apply to the electron:
% Quadratic Well = 1
% Coulomb Interaction = 2
% Screened Coulomb Interaction = 3
% Square Well = 4  ----- NOTE: uses ODEs defined in ChemPhysLetters (CPL) paper
% Morse Potential (from CPL) = 5
potential_operator_idx = 2;

close all; % close any open figures
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
            [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = initialize_coulombInteraction(num_particles, simulation_steps, A_CI, x_0_CI, gamma_0_CI);
            
            % compute the first forces pPRIME and etaPRIME given the
            % Coulomb interaction operator
            pPRIME_acc = compute_force_pPRIME_coulombInteraction(pPRIME_acc, num_particles, Z_CI, q_pos, gamma_packet_width, simulation_step);
            etaPRIME_acc = compute_force_etaPRIME_coulombInteraction(etaPRIME_acc, num_particles, A_CI, Z_CI, reduced_planck_constant, mass, gamma_packet_width, q_pos, simulation_step);
        elseif potential_operator_idx == 3
            [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = initialize_screenedCoulombInteraction(num_particles, num_dimensions, simulation_steps, q_0_SCI, gamma_0_SCI);
            
            pPRIME_acc = compute_force_pPRIME_screenedCoulombInteraction(pPRIME_acc, num_dimensions, num_particles, q_pos, simulation_step, lambda_SCI, Z_SCI, e_SCI, gamma_packet_width);
            etaPRIME_acc = compute_force_etaPRIME_screenedCoulombInteraction(etaPRIME_acc, num_particles, num_dimensions, gamma_packet_width, simulation_step, mass, reduced_planck_constant, Z_SCI, e_SCI, q_pos, lambda_SCI);
        elseif potential_operator_idx == 4
            [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = initialize_squareWell(simulation_steps, q_0_SW, gamma_0_SW);
            
            pPRIME_acc = compute_force_pPRIME_squareWell(pPRIME_acc, q_pos, gamma_packet_width, V0_SW, a_SW, simulation_step);
            etaPRIME_acc = compute_force_etaPRIME_squareWell(etaPRIME_acc, gamma_packet_width, q_pos, mass, V0_SW, a_SW, simulation_step);
        
        elseif potential_operator_idx == 5
            [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = initialize_morse_CPL(simulation_steps, q_0_MP_CPL, gamma_0_MP_CPL);
            
            pPRIME_acc = compute_force_pPRIME_morse_CPL(pPRIME_acc, q_pos, gamma_packet_width, D_MP_CPL, a_MP_CPL, b_MP_CPL, simulation_step);
            etaPRIME_acc = compute_force_etaPRIME_morse_CPL(etaPRIME_acc, gamma_packet_width, q_pos, mass, D_MP_CPL, a_MP_CPL, b_MP_CPL, simulation_step);
        else
            disp('No potential operator index given, or incorrectly set. Returning');
            return;
        end
    else
        % update the positions, velocities and acceleration using the
        % velocity verlet algorithm
        [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = ...
            velocity_verlet(q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc, ...
            mass, delta_t, epsilon, reduced_planck_constant, A_CI, Z_CI, Z_SCI, e_CI, e_SCI, lambda_SCI, V0_SW, a_SW, D_MP_CPL, a_MP_CPL, b_MP_CPL,...
            simulation_step,...
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
if potential_operator_idx == 1 % quadratic well
    visualize_quadraticWell(t, q_pos, gamma_packet_width, simulation_steps, A_0_well, omega_well, gamma_0_well, epsilon);
elseif potential_operator_idx == 2 % coulomb interaction
    visualize_coulomb(t, q_pos, gamma_packet_width, Z_CI, A_CI, e_CI, visualization_toggles);
elseif potential_operator_idx == 3 % screened coulomb interaction
    visualize_screenedCoulomb(t, q_pos, gamma_packet_width, q_0_SCI, gamma_0_SCI, Z_SCI, e_SCI, lambda_SCI, simulation_steps, visualization_toggles);
elseif potential_operator_idx == 4 % square well
    visualize_squareWell(t, q_pos, gamma_packet_width, V0_SW, a_SW, mass, visualization_toggles);
elseif potential_operator_idx == 5 % morse potential (CPL)
    visualize_morse_CPL(t, q_pos, gamma_packet_width, D_MP_CPL, a_MP_CPL, b_MP_CPL, mass, visualization_toggles);
end
