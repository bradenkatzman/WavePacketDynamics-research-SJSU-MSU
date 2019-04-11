function [] = visualize_screenedCoulomb(t, q_pos, gamma_packet_width, ...
    q_0, gamma_0, Z, e, lambda,...
    simulation_steps, visualization_toggles)

    % make a grid of the magnitude of q and gamma and compute the potential
    %q_norm = sqrt(sum(q_pos.^2, 2));
    %[qX, gammaY] = meshgrid(reshape(q_norm(1, 1, :), 1, simulation_steps), gamma_packet_width);
    [X, Y] = meshgrid(.01:.05:10.01, .01:.05:10.01);
    V = compute_screenedCoulomb_potential(X, Y, Z, e, lambda);

    if visualization_toggles(1) == 1
        % plot q over time
        figure(1);
        plot(t, reshape(q_pos(1, 1, :), 1, simulation_steps));
        xlabel('time');
        ylabel('q1');
        title(['Screened Coulomb: q_01=', num2str(q_0(1,1)), ', Z=', num2str(Z), ', e=', num2str(e), ', \lambda', '=', num2str(lambda)]);
        grid on;
        
         figure(2);
         plot(t, reshape(q_pos(2, 1, :), 1, simulation_steps));
         xlabel('time');
         ylabel('q2');
         title(['Screened Coulomb: q_02=', num2str(q_0(2,1)), ', Z=', num2str(Z), ', e=', num2str(e), ', \lambda', '=', num2str(lambda)]);
         grid on;
        
         figure(3);
         plot(t, reshape(q_pos(3, 1, :), 1, simulation_steps));
         xlabel('time');
         ylabel('q3');
         title(['Screened Coulomb: q_03=', num2str(q_0(3,1)), ', Z=', num2str(Z), ', e=', num2str(e), ', \lambda', '=', num2str(lambda)]);
         grid on;
    end
    
    if visualization_toggles(2) == 1
        % plot gamma over time
        figure(4);
        plot(t, gamma_packet_width);
        xlabel('time');
        ylabel('\gamma');
        title(['Screened Coulomb: \gamma_0= ', num2str(gamma_0), ', Z=', num2str(Z), ', e=', num2str(e), ', \lambda', '=', num2str(lambda)]);
        grid on;
    end
    
    if visualization_toggles(3) == 1
        % plot q (component 1) against gamma
        figure(5);
        plot(reshape(q_pos(1, 1, :), 1, simulation_steps), gamma_packet_width);
        hold on
        plot(q_0(1,1), gamma_0, 'r*');
        text(q_0(1,1),gamma_0,'(q_01, \gamma_0)','VerticalAlignment','top','HorizontalAlignment','left')
        xlabel('x', 'FontSize', 15);
        ylabel('\gamma', 'FontSize', 25);
        title(['Screened Coulomb: q_01 vs. \gamma',', Z=', num2str(Z), ', e=', num2str(e), ', \lambda', '=', num2str(lambda)]);
        grid on;
        
        figure(6);
        plot(reshape(q_pos(2, 1, :), 1, simulation_steps), gamma_packet_width);
        hold on
        plot(q_0(2,1), gamma_0, 'r*');
        text(q_0(2,1),gamma_0,'(q_03, \gamma_0)','VerticalAlignment','top','HorizontalAlignment','left')
        xlabel('x', 'FontSize', 15);
        ylabel('\gamma', 'FontSize', 25);
        title(['Screened Coulomb: q_02 vs. \gamma',', Z=', num2str(Z), ', e=', num2str(e), ', \lambda', '=', num2str(lambda)]);
        grid on;
        
        figure(7);
        plot(reshape(q_pos(3, 1, :), 1, simulation_steps), gamma_packet_width);
        hold on
        plot(q_0(3,1), gamma_0, 'r*');
        text(q_0(3,1),gamma_0,'(q_03, \gamma_0)','VerticalAlignment','top','HorizontalAlignment','left')
        xlabel('x', 'FontSize', 15);
        ylabel('\gamma', 'FontSize', 25);
        title(['Screened Coulomb: q_03 vs. \gamma',', Z=', num2str(Z), ', e=', num2str(e), ', \lambda', '=', num2str(lambda)]);
        grid on;
    end
    
    if visualization_toggles(4) == 1
        % plot the isocontours of V(q_norm, gamma)        
        figure(8);
        contourf(X, Y, V, 20);
        shading interp;
        xlabel('q (vector-magnitude)', 'FontSize', 15);
        ylabel('\gamma', 'FontSize', 25);
        title(['Screened Coulomb: Isocontours of V(q1, \gamma) with Z=', num2str(Z), ', e=', num2str(e), ', \lambda', '=', num2str(lambda)]);
    end
    
    if visualization_toggles(5) == 1
        % plot the surface of V
        figure(9);
        surf(X, Y, V);
        shading interp;
        title(['Coulomb: Surface of V(q, \gamma) with Z=', num2str(Z), ', e=', num2str(e), ', \lambda', '=', num2str(lambda)]);
    end
    
    if visualization_toggles(6) == 1
        % plot the pseudocolor plot
        figure(10);
        pcolor(X, Y, V);
        shading interp;
        title(['Coulomb: Pseudocolor plot (with shading interp) of V(q, \gamma) with Z=', num2str(Z), ', e=', num2str(e), ', \lambda', '=', num2str(lambda)]);
    end
end