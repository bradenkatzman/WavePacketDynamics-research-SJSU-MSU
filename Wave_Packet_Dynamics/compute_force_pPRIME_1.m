% Wave Packet Dynamics
% Calculate pPRIME for the potential operator:
% dp/dt = -del(q)V
% V = epsilon(q^2 + gamma^2)
% --> del(q)V = 2*epsilon*q

function pPRIME_acc = compute_force_pPRIME_1(pPRIME_acc, num_dimensions, num_particles,...
    q_pos, simulation_step,...
    epsilon)
    
    for particle_itr = 1:num_particles
       % dp/dt at t = -2*epsilon / mass * q(t) 
       pPRIME_acc(1:num_dimensions, particle_itr, simulation_step) = ...
           (-2*epsilon) * q_pos(1:num_dimensions, particle_itr, simulation_step);
    end
end