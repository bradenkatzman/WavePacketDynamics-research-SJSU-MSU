% Wave Packet Dynamics

function [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = ...
    initialize_squareWell(simulation_steps, q_0, gamma_0)
    
    % POSITION - start at 1,0,0 --> really just 1 because the 2nd and 3rd
    % components will always be 0 so we're not going to keep track of that
    % so that we don't have dimensionality issues
    q_pos = zeros(1, simulation_steps);
    % make the starting point (1,0,0)
    q_pos(1, 1) = q_0;
    
    % VELOCITY - start at 0,0,0
    p_vel = zeros(1, simulation_steps);
    
    % ACCELERATION - start at 0,0,0
    pPRIME_acc = zeros(1, simulation_steps);
    
    % PACKET WIDTH - start at sqrt(pi/6) * A
    gamma_packet_width = zeros(1, simulation_steps);
    gamma_packet_width(1, 1) = gamma_0;
    
    % PACKET MOMEMTUM - start at 0
    eta_packet_momentum = zeros(1, simulation_steps);
    etaPRIME_acc = zeros(1, simulation_steps);
end