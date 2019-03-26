%
% Compute eta prime for the square well operator
% NOTE: This is computed according to the system of ODEs
% presented in the ChemPhysLetter paper
%
% NOTE: This paper solves all problems in one spatial dimension

function etaPRIME_acc = compute_force_etaPRIME_squareWell(etaPRIME_acc, ...
    gamma_packet_width, q_pos, ...
    m, V0, a, ...
    simulation_step)

    q = q_pos(1, simulation_step);
    gamma = gamma_packet_width(1, simulation_step);
    
    etaPRIME_acc(1, simulation_step) = -1 * ...
        ((-1/(4*m*gamma^3)) + ((V0/sqrt(2*pi)*gamma^2) * (((a + q)*exp((-1)*((a + q)/(sqrt(2)*gamma)))^2) + (exp((-1)*((a - q)/(sqrt(2)*gamma)))^2))));
end