function V = compute_Morse_potential_CPL(q, gamma, ...
    D, a, b, m)

    V = (1./(8*m.*gamma.^2)) + D.*(exp(-2*b.*(q - a - gamma.^2 .* b)) - ...
        2 .* exp(-b .* (q - a - (1/2).* gamma.^2 .* b)));

end