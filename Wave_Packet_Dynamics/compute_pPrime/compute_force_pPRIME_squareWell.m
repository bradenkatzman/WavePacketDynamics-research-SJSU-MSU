%
% Compute p prime for the square well operator
% NOTE: This is computed according to the system of ODEs
% presented in the ChemPhysLetter paper (in the case of p',
% the equations are the same between their paper in ours.
% However, the equation for eta' for the square well operator
% differs so we note that here for consistency, since the
% two functions are coupled
%
% NOTE: This paper solves all problems in one spatial dimension

function pPRIME_acc = compute_force_pPRIME_squareWell(pPRIME_acc,...
    q_pos, gamma_packet_width, ...
    V0, a, ...
    simulation_step)


    q = q_pos(1, simulation_step);
    gamma = gamma_packet_width(1, simulation_step);

    pPRIME_acc(1, simulation_step) = -1 * ...
        ((V0/(sqrt(2*pi) * gamma)) * (exp(-((a - q)/(sqrt(2) * gamma))^2) - (exp(-((a + q)/(sqrt(2) * gamma))^2))));
    
end