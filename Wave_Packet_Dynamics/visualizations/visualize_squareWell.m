% Visualization and plotting
function [] = visualize_squareWell(t, q_pos, gamma_packet_width, ...
    V0, a, m, ...
    visualization_toggles)

    [X, Y] = meshgrid(q_pos, gamma_packet_width);
    V = compute_squareWell_potential(X, Y, V0, a, m);

    if visualization_toggles(1) == 1
        % plot q over time
        figure(1);
        plot(t, q_pos);
        xlabel('time');
        ylabel('q1');
        title(['Square Well: q (position), V0=', num2str(V0), ', a=', num2str(a)]);
        grid on;
    end
    
    if visualization_toggles(2) == 1
        % plot gamma over time
        figure(2);
        plot(t, gamma_packet_width);
        xlabel('time');
        ylabel('gamma');
        title(['Square Well: gamma - packet width, V0=', num2str(V0), ', a=', num2str(a)]);
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
        title(['Square Well: x vs. gamma, , V0=', num2str(V0), ', a=', num2str(a)]);
        grid on;
    end
    
    if visualization_toggles(4) == 1
        % plot the isocontours of V(q, gamma)
        figure(4);
        contourf(X, Y, V, 20);
        xlabel('q', 'FontSize', 15);
        ylabel('\gamma', 'FontSize', 25);
        title(['Square Well: Isocontours of V(q, \gamma) with V0=', num2str(V0), ', a=', num2str(a)]);
    end
    
    if visualization_toggles(5) == 1
        % plot the surface of V
        figure(5);
        surf(X, Y, V);
        shading interp;
        title(['Square Well: Surface of V(q, \gamma) with , V0=', num2str(V0), ', a=', num2str(a)]);
    end
    
    if visualization_toggles(6) == 1
        % plot the pseudocolor plot
        figure(6);
        pcolor(X, Y, V);
        shading interp;
        title(['Square Well: Pseudocolor plot of V(q, \gamma) with , V0=', num2str(V0), ', a=', num2str(a)]);
    end
end