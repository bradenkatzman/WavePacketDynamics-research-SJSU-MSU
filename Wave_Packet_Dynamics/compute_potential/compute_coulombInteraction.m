function V = compute_coulombInteraction(q, gamma, ...
    Z, e)

    V = -1 * ((Z*e)./q) * erf(sqrt(3/2) .* q./gamma); 
end