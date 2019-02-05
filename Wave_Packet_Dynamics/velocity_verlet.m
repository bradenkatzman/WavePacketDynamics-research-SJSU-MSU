% Update the positions, velocities, accelerations, wave packet widths, and
% wave packet momentums of the particles using the Velocity Verlet algorithm
%
% Parameters:
%
% Velocity Verlet Algorithm
%    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
%    a(t+dt) =  
%    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt

function [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = ...
velocity_verlet(q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc,...
    mass, delta_t, epsilon, reduced_planck_constant, Z, simulation_step, ...
    num_dimensions, num_particles, ...
    potential_operator_idx)

    for particle_itr = 1:num_particles

        % compute the new position q -- q(t+1)
        q_pos(1:num_dimensions, particle_itr, simulation_step) = ...
            q_pos(1:num_dimensions, particle_itr, simulation_step-1) + ...
                (delta_t / mass) * p_vel(1:num_dimensions, particle_itr, simulation_step-1) + ...
                    0.5 * (delta_t * delta_t) * pPRIME_acc(1:num_dimensions, particle_itr, simulation_step-1);

        % compute the new acceleration pPRIME -- pPRIME(t+1)
        if potential_operator_idx == 1
           pPRIME_acc = compute_force_pPRIME_quadraticWell(pPRIME_acc, num_dimensions, num_particles,...
                                                                                q_pos, simulation_step, epsilon);
        elseif potential_operator_idx == 2
            pPRIME_acc = compute_force_pPRIME_coulombInteraction(pPRIME_acc, num_dimensions, num_particles, Z, q_pos, gamma_packet_width, simulation_step);
        else
            disp('Potential operator idx not correctly set. Returning');
            return;
        end

        % compute the new velocity p -- p(t+1)
        p_vel(1:num_dimensions, particle_itr, simulation_step) = ...
            p_vel(1:num_dimensions, particle_itr, simulation_step-1) + ...
            0.5 * delta_t * (pPRIME_acc(1:num_dimensions, particle_itr, simulation_step-1) + pPRIME_acc(1:num_dimensions, particle_itr, simulation_step));


        % compute the new packet width gamma -- gamma(t+1)
        gamma_packet_width(particle_itr, simulation_step) = gamma_packet_width(particle_itr, simulation_step-1) + ...
            delta_t * eta_packet_momentum(particle_itr, simulation_step-1) + ...
            0.5 * (delta_t * delta_t) * etaPRIME_acc(particle_itr, simulation_step-1);

        % compute the new acceleration etaPRIME -- etaPRIME(t+1)
        if potential_operator_idx == 1
            etaPRIME_acc = compute_force_etaPRIME_quadraticWell(etaPRIME_acc, num_particles, ...
                                                                    gamma_packet_width, simulation_step, epsilon, mass, reduced_planck_constant);
        elseif potential_operator_idx == 2
            etaPRIME_acc = compute_force_ePRIME_coulombInteraction(etaPRIME_acc, num_particles, Z, reduced_planck_constant, mass, ...
                                                                     gamma_packet_width, q_pos, simulation_step);
        else
            disp('Potential operator idx not correctly set. Returning');
            return;
        end

        % compute the new velocity eta -- eta_packet_momentum(t+1)
        eta_packet_momentum(particle_itr, simulation_step) = eta_packet_momentum(particle_itr, simulation_step-1) + ...
            0.5 * delta_t * (etaPRIME_acc(particle_itr, simulation_step-1) + etaPRIME_acc(particle_itr, simulation_step));
    end

end