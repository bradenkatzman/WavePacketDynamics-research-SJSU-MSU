% Visualization and plotting
function [] = visualize_coulomb(t, q_pos, gamma_packet_width, ...
    Z, A, e, ...
    visualization_toggles)

        [X, Y] = meshgrid(-2:.1:2, .1:.1:2);
        V = ((-Z*e) ./ X) .* erf((sqrt(3/2)) * (X ./ Y));

    if visualization_toggles(1) == 1
        % plot q over time
        figure(1);
        plot(t, q_pos);
        xlabel('time');
        ylabel('q1');
        title(['Coulomb: x (position), Z=', num2str(Z)]);
        grid on;
    end
    
    if visualization_toggles(2) == 1
        % plot gamma over time
        figure(2);
        plot(t, gamma_packet_width);
        xlabel('time');
        ylabel('gamma');
        title(['Coulomb: gamma - packet width, Z=', num2str(Z)]);
        grid on;
    end
    
    if visualization_toggles(3) == 1
        % plot q against gamma
        figure(3);
        plot(q_pos, gamma_packet_width);
        hold on
        plot(q_pos(1,1), gamma_packet_width(1,1), 'r*');
        text(q_pos(1,1),gamma_packet_width(1,1),'(x_0, \gamma_0)','VerticalAlignment','top','HorizontalAlignment','left')
        xlabel('x', 'FontSize', 15);
        ylabel('\gamma', 'FontSize', 25);
        title(['Coulomb: x vs. gamma, A=', num2str(A)]);
        grid on;
    end
    
    if visualization_toggles(4) == 1
        % plot the isocontours of V(q, gamma)
       
        figure(4);
        contourf(X, Y, V, 20);
        xlabel('q', 'FontSize', 15);
        ylabel('\gamma', 'FontSize', 25);
        title(['Coulomb: Isocontours of V(q, \gamma) with A=', num2str(A)]);
    end
    
    if visualization_toggles(5) == 1
        % plot the surface of V
        figure(5);
        surf(X, Y, V);
        shading interp;
        title(['Coulomb: Surface of V(q, \gamma) with A=', num2str(A)]);
    end
    
    if visualization_toggles(6) == 1
        % plot the pseudocolor plot
        figure(6);
        pcolor(X, Y, V);
        shading interp;
        title(['Coulomb: Pseudocolor plot (with shading interp) of V(q, \gamma) with A=', num2str(A)]);
    end
end