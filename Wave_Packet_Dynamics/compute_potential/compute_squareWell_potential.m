%
% Eqn (30)
%
function V = compute_squareWell_potential(q, gamma, ...
    V0, a, m)

    V = (1./(8*m.*gamma.^2)) - (0.5*V0).*(erf(((a + q)./sqrt(2)) .* gamma) + erf(((a - q)./sqrt(2)) .* gamma));
end