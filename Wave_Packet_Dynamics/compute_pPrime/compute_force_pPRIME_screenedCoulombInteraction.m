% Note: q = ||q_vec||

function pPRIME_acc = compute_force_pPRIME_screenedCoulombInteraction(pPRIME_acc, num_dimensions, num_particles,...
    q_pos, simulation_step,...
    lambda, Z, e, gamma_packet_width)

    for particle_itr = 1:num_particles
        % compute the norm now to cut down on the size of the equation
        norm_q = norm(q_pos(1:num_dimensions, particle_itr, simulation_step));
        
        gamma = gamma_packet_width(particle_itr, simulation_step);
        q_bold = q_pos(1:num_dimensions, particle_itr, simulation_step);
        
        pPRIME_acc(1:num_dimensions, particle_itr, simulation_step) = -1 * ...
            (norm_q + lambda) .* ((Z*e^2)./(2*lambda.*norm_q.^3)) .* erfc((gamma)./(sqrt(6)*lambda) - (sqrt(3/2).*(norm_q./gamma))) .* e^(-norm_q./lambda + ((gamma.^2)./(6.*lambda^2))) .* q_bold - ...
            (norm_q - lambda) .* ((Z*e^2)./(2*lambda.*norm_q.^3)) .* erfc((gamma)./(sqrt(6)*lambda) + (sqrt(3/2).*(norm_q./gamma))) .* e^(norm_q./lambda + ((gamma.^2)./(6*lambda^2))) .* q_bold + ...
            sqrt(6/pi) * ((Z*e^2)/(gamma.*norm_q.^2)) * e^((-3 .* norm_q.^2)./(2 .* gamma.^2)) .* q_bold;
    end

end