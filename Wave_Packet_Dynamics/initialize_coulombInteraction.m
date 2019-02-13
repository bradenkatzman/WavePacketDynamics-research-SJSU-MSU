% Wave Packet Dynamics

function [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = ...
    initialize_coulombInteraction(num_particles, simulation_steps, A)
    
    % POSITION - start at 1,0,0 --> really just 1 because the 2nd and 3rd
    % components will always be 0 so we're not going to keep track of that
    % so that we don't have dimensionality issues
    q_pos = zeros(num_particles, simulation_steps);
    % make the starting point (1,0,0)
    q_pos(1, 1) = 0.001;
    
    % VELOCITY - start at 0,0,0
    p_vel = zeros(num_particles, simulation_steps);
    
    % ACCELERATION - start at 0,0,0
    pPRIME_acc = zeros(num_particles, simulation_steps);
    
    % PACKET WIDTH - start at 
    gamma_packet_width = zeros(num_particles, simulation_steps);
    gamma_packet_width(1, 1) = sqrt(pi/6) * A;
    
    % PACKET MOMEMTUM - start at 0
    eta_packet_momentum = zeros(num_particles, simulation_steps);
    etaPRIME_acc = zeros(num_particles, simulation_steps);
end