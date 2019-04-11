% Update the positions, velocities, accelerations, wave packet widths, and
% wave packet momentums of the particles using the Velocity Verlet algorithm
%
% Parameters:
%
% Velocity Verlet Algorithm
%    1. x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
%    2. a(t+dt) = however the forces are computed given the potential surface supplied 
%    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt

function [q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc] = ...
velocity_verlet(q_pos, p_vel, pPRIME_acc, gamma_packet_width, eta_packet_momentum, etaPRIME_acc,...
    mass, delta_t, epsilon, reduced_planck_constant, A_CI, Z_CI, Z_SCI, e_CI, e_SCI, lambda_SCI, V0_SW, a_SW, D_MP_CPL, a_MP_CPL, b_MP_CPL,...
    simulation_step, ...
    num_dimensions, num_particles, ...
    potential_operator_idx)

    for particle_itr = 1:num_particles
        % ---- VELOCITY VERLET STEP 1 --------
        % compute the new position q -- q(t+1)
        if potential_operator_idx == 1 || potential_operator_idx == 3 % position is in three dimensions
            q_pos(1:num_dimensions, particle_itr, simulation_step) = ...
                q_pos(1:num_dimensions, particle_itr, simulation_step-1) + ...
                    (delta_t / mass) * p_vel(1:num_dimensions, particle_itr, simulation_step-1) + ...
                        0.5 * (delta_t * delta_t) * pPRIME_acc(1:num_dimensions, particle_itr, simulation_step-1);
        elseif potential_operator_idx == 2 || potential_operator_idx == 4 || potential_operator_idx == 5 % position is in one dimension
             q_pos(particle_itr, simulation_step) = ...
                q_pos(particle_itr, simulation_step-1) + ...
                    (delta_t / mass) * p_vel(particle_itr, simulation_step-1) + ...
                        0.5 * (delta_t * delta_t) * pPRIME_acc(particle_itr, simulation_step-1);
        end
        
        % compute the new packet width gamma -- gamma(t+1)
            gamma_packet_width(particle_itr, simulation_step) = gamma_packet_width(particle_itr, simulation_step-1) + ...
                delta_t * eta_packet_momentum(particle_itr, simulation_step-1) + ...
                0.5 * (delta_t * delta_t) * etaPRIME_acc(particle_itr, simulation_step-1);
        % --------------------- END VELOCITY VERLET STEP 1 ------------------------------

        % ------------ VELOCITY VERLET STEP 2 ------------------------
        % compute the new acceleration pPRIME -- pPRIME(t+1)
        if potential_operator_idx == 1
           pPRIME_acc = compute_force_pPRIME_quadraticWell(pPRIME_acc, num_dimensions, num_particles,...
                                                                                q_pos, simulation_step, epsilon);
        elseif potential_operator_idx == 2
            pPRIME_acc = compute_force_pPRIME_coulombInteraction(pPRIME_acc, num_particles, Z_CI, q_pos, gamma_packet_width, simulation_step);
        elseif potential_operator_idx == 3
            pPRIME_acc = compute_force_pPRIME_screenedCoulombInteraction(pPRIME_acc, num_dimensions, num_particles, ...
                                                                                 q_pos, simulation_step, ...
                                                                                 lambda_SCI, Z_SCI, e_SCI, gamma_packet_width);
        elseif potential_operator_idx == 4
            pPRIME_acc = compute_force_pPRIME_squareWell(pPRIME_acc, q_pos, gamma_packet_width, V0_SW, a_SW, simulation_step);
        elseif potential_operator_idx == 5
            pPRIME_acc = compute_force_pPRIME_morse_CPL(pPRIME_acc, q_pos, gamma_packet_width, D_MP_CPL, a_MP_CPL, b_MP_CPL, simulation_step);
        else
            disp('Potential operator idx not correctly set. Returning');
            return;
        end
        
        % compute the new acceleration etaPRIME -- etaPRIME(t+1)
        if potential_operator_idx == 1
            etaPRIME_acc = compute_force_etaPRIME_quadraticWell(etaPRIME_acc, num_particles, ...
                                                                    gamma_packet_width, simulation_step, epsilon, mass, reduced_planck_constant);
        elseif potential_operator_idx == 2
            etaPRIME_acc = compute_force_etaPRIME_coulombInteraction(etaPRIME_acc, num_particles, A_CI, Z_CI, reduced_planck_constant, mass, ...
                                                                     gamma_packet_width, q_pos, simulation_step);
        elseif potential_operator_idx == 3
            etaPRIME_acc = compute_force_etaPRIME_screenedCoulombInteraction(etaPRIME_acc, num_particles, num_dimensions,...
                                                                    gamma_packet_width, simulation_step, mass, reduced_planck_constant, Z_SCI, e_SCI, q_pos, lambda_SCI);
        elseif potential_operator_idx == 4
            etaPRIME_acc = compute_force_etaPRIME_squareWell(etaPRIME_acc, gamma_packet_width, q_pos, mass, V0_SW, a_SW, simulation_step);
        elseif potential_operator_idx == 5
            etaPRIME_acc = compute_force_etaPRIME_morse_CPL(etaPRIME_acc, gamma_packet_width, q_pos, mass, D_MP_CPL, a_MP_CPL, b_MP_CPL, simulation_step);
        else
            disp('Potential operator idx not correctly set. Returning');
            return;
        end
        % ------------------------ END VELOCITY VERLET STEP 2 -------------

        % --------------- VELOCITY VERLET STEP 3 --------------------
        if potential_operator_idx == 1 || potential_operator_idx == 3 % we have to check the operator here because the quad well is a 3D problem and Coulomb is 1D
            % compute the new velocity p -- p(t+1)
            p_vel(1:num_dimensions, particle_itr, simulation_step) = ...
                p_vel(1:num_dimensions, particle_itr, simulation_step-1) + ...
                0.5 * delta_t * (pPRIME_acc(1:num_dimensions, particle_itr, simulation_step-1) + pPRIME_acc(1:num_dimensions, particle_itr, simulation_step));
        elseif potential_operator_idx == 2 || potential_operator_idx == 4 || potential_operator_idx == 5
            % compute the new velocity p -- p(t+1)
            p_vel(particle_itr, simulation_step) = ...
                p_vel(particle_itr, simulation_step-1) + ...
                0.5 * delta_t * (pPRIME_acc(particle_itr, simulation_step-1) + pPRIME_acc(particle_itr, simulation_step));
        end

        % compute the new velocity eta -- eta_packet_momentum(t+1)
        eta_packet_momentum(particle_itr, simulation_step) = eta_packet_momentum(particle_itr, simulation_step-1) + ...
            0.5 * delta_t * (etaPRIME_acc(particle_itr, simulation_step-1) + etaPRIME_acc(particle_itr, simulation_step));
        % ----------------- END VELOCITY VERLET STEP 3 --------------
    end
end