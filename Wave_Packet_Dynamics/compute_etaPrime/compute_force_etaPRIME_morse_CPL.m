%
% Compute eta prime for the Morse potential
% NOTE: this is computed according to the system of ODEs
% presented in the ChemPhysLetters paper
%
% NOTE: this paper solves all problems in one spatial dimension

function etaPRIME_acc = compute_force_etaPRIME_morse_CPL(etaPRIME_acc,...
    gamma_packet_width, q_pos, ...
    m, D, a, b, ...
    simulation_step)
    
    q = q_pos(1, simulation_step);
    gamma = gamma_packet_width(1, simulation_step);
    
    etaPRIME_acc(1, simulation_step) = -1 * ...
        ((-1/(4*m*(gamma^3))) + ...
        (2*D*(b^2)*gamma)*(2*exp((-2*b*q) + (2*b*a) + (2*(b^2)*(gamma^2))) - exp((-b*q) + (b*a) + (.5*(b^2)*(gamma^2)))));

end