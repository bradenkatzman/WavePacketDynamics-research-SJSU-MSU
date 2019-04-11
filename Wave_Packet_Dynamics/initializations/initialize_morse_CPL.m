function [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = ...
    initialize_morse_CPL(simulation_steps, q_0, gamma_0)

    % position
    q_pos = zeros(1, simulation_steps);
    q_pos(1, 1) = q_0;
    
    % velocity
    p_vel = zeros(1, simulation_steps);
    
    % acceleration
    pPRIME_acc = zeros(1, simulation_steps);
    
    % packet width
    gamma_packet_width = zeros(1, simulation_steps);
    gamma_packet_width(1, 1) = gamma_0;
    
    % packet momentum
    eta_packet_momentum = zeros(1, simulation_steps);
    
    % packet acceleration
    etaPRIME_acc = zeros(1, simulation_steps);
end