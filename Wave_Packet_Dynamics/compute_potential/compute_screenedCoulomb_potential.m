%
% Eqn (30)
%
function V = compute_screenedCoulomb_potential(q, gamma, ...
    Z, e, lambda)

    V = -1 .* ((Z*e^2) ./ (2.*q)) * e.^((-q./lambda) + ((gamma.^2)./(6*lambda^2))) .* erfc((gamma./(sqrt(6).*lambda)) - (sqrt(3/2) .* (q ./ gamma))) + ...
        ((Z*e^2) ./ (2.*q)) .* e.^((q./lambda) + ((gamma.^2)/(6.*lambda.^2))) .* erfc((gamma./(sqrt(6).*lambda)) + (sqrt(3/2) .* (q ./ gamma)));

end