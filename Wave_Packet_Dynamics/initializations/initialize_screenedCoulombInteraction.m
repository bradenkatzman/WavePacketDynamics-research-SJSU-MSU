% Wave Packet Dynamics

function [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = ...
    initialize_screenedCoulombInteraction(num_particles, num_dimensions, simulation_steps,...
    q_0, gamma_0)
    
    % POSITION
    q_pos = zeros(num_dimensions, num_particles, simulation_steps);
    q_pos(1:num_dimensions, 1:num_particles, 1) = q_0;
    
    % VELOCITY - start at 0,0,0
    p_vel = zeros(num_dimensions, num_particles, simulation_steps);
    
    % ACCELERATION - start at 0,0,0
    pPRIME_acc = zeros(num_dimensions, num_particles, simulation_steps);
    
    % PACKET WIDTH - start at gamma_0
    gamma_packet_width = zeros(num_particles, simulation_steps);
    gamma_packet_width(1:num_particles, 1) = gamma_0;
    
    % PACKET MOMEMTUM - start at 0
    eta_packet_momentum = zeros(num_particles, simulation_steps);
    etaPRIME_acc = zeros(num_particles, simulation_steps);
end