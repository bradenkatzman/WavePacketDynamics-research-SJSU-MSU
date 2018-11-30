% Wave Packet Dynamics

function [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = ...
    initialize(num_particles, num_dimensions, simulation_steps)
    
    % POSITION - start at 1,1,1
    q_pos = zeros(num_dimensions, num_particles, simulation_steps);
    % make the starting point (1,1,1)
    q_pos(1:num_dimensions, 1:num_particles, 1) = 1;
    
    % VELOCITY - start at 0,0,0
    p_vel = zeros(num_dimensions, num_particles, simulation_steps);
    
    % ACCELERATION - start at 0,0,0
    pPRIME_acc = zeros(num_dimensions, num_particles, simulation_steps);
    
    % PACKET WIDTH - start at 1
    gamma_packet_width = zeros(num_particles, simulation_steps);
    gamma_packet_width(1:num_particles, 1) = 1;
    
    % PACKET MOMEMTUM
    eta_packet_momentum = zeros(num_particles, simulation_steps);
    etaPRIME_acc = zeros(num_particles, simulation_steps);
end