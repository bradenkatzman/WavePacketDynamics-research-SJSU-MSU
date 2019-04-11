function [] = visualize_morse_CPL(t, q_pos, gamma_packet_width, ...
    D, a, b, m, ...
    visualization_toggles)

    [X, Y] = meshgrid(-5:.1:5, .1:.1:5);
    V = compute_Morse_potential_CPL(X, Y, D, a, b, m);
    
    if visualization_toggles(1) == 1
        % plot q over time
        figure(1);
        plot(t, q_pos);
        xlabel('time');
        ylabel('q_1');
        title(['Morse Potential (CPL): q with D=', num2str(D), ' a=', num2str(a), ' b=', num2str(b)]);
    end
    
    if visualization_toggles(2) == 1
       % plot gamma over time
       figure(2);
       plot(t, gamma_packet_width);
       xlabel('time');
       ylabel('\gamma');
       title(['Morse Potential (CPL): \gamma with D=', num2str(D), ' a=', num2str(a), ' b=', num2str(b)]);
    end
    
    if visualization_toggles(3) == 1
       % plot q against gamma
       figure(3);
       plot(q_pos, gamma_packet_width);
       hold on
       plot(q_pos(1,1), gamma_packet_width(1,1), 'r*');
       text(q_pos(1,1), gamma_packet_width(1,1), '(q_0, \gamma_0)', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
       xlabel('q', 'FontSize', 25);
       ylabel('\gamma', 'FontSize', 25);
       title(['Morse Potential (CPL): q vs. gamma with D=', num2str(D), ' a=', num2str(a), ' b=', num2str(b)]);
    end
    
    if visualization_toggles(4) == 1
       % plot the isocontours of V(q, gamma)
       figure(4);
       contourf(X, Y, V, 20);
       xlabel('q', 'FontSize', 15);
       ylabel('\gamma', 'FontSize', 25);
       title(['Morse Potential (CPL): Isocontours of V(q, \gamma) with D=', num2str(D), ' a=', num2str(a), ' b=', num2str(b)]);
    end

    if visualization_toggles(5) == 1
        % plot the surface of V
        figure(5);
        surf(X,Y,V);
        shading interp;
        title(['Morse Potential (CPL): Surface of V(q, \gamma) with D=', num2str(D), ' a=', num2str(a), ' b=', num2str(b)]);
    end
    
    if visualization_toggles(6) == 1
        % plot the pseudocolor plot
        figure(6);
        pcolor(X, Y, V);
        shading interp;
        title(['Morse Potential (CPL): Pseudocolor plot of V(q, \gamma) with D=', num2str(D), ' a=', num2str(a), ' b=', num2str(b)]);
end