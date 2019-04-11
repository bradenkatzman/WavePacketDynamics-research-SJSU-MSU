%
% Compute p prime for the Morse potential
% NOTE: this is computed according to the system of ODEs
% presented in the ChemPhysLetters paper
%
% NOTE: this paper solves all problems in one spatial dimension

function pPRIME_acc = compute_force_pPRIME_morse_CPL(pPRIME_acc,...
    q_pos, gamma_packet_width, ...
    D, a, b, ...
    simulation_step)

    q = q_pos(1, simulation_step);
    gamma = gamma_packet_width(1, simulation_step);
    
    pPRIME_acc(1, simulation_step) = -1 * ...
        ((2*D*b) * (exp((-b*q) + (b*a) + (.5*(b^2)*(gamma^2))) - exp((-2*b*q) + (2*b*a) + (2*(b^2)*(gamma^2)))));

end