
function [] = visualize_quadraticWell(t, q_pos, gamma_packet_width, simulation_steps, ...
    A_0, omega, gamma_0, eps)

    [X, Y] = meshgrid(-1:.1:1, 1:.01:1.2);
    V = eps .* (X.^2 + Y.^2);

    figure(1);
    plot(t, reshape(q_pos(1, 1, :), 1, simulation_steps)); % need to reduce the dimensionality
    hold on
    plot(t, cos(omega*t));
    hold off
    xlabel('time');
    ylabel('q1');
    title('Quadratic Well Operator: q1 - position, component 1');
    drawnow

    figure(2);
    plot(t, reshape(q_pos(2, 1, :), 1, simulation_steps));
    xlabel('time');
    ylabel('q2');
    title('Quadratic Well Operator: q2 - position, component 2');
    drawnow

    figure(3);
    plot(t, reshape(q_pos(3, 1, :), 1, simulation_steps));
    xlabel('time');
    ylabel('q3');
    title('Quadratic Well Operator: q3 - position, component 3');
    drawnow


    % plot gamma
    figure(4);
    plot(t, gamma_packet_width);
    hold on
    plot(t, sqrt(...
        (A_0 * sin(omega * (t)).^2) + ...
        (gamma_0 * gamma_0 * cos(omega * t).^2 )));
    hold off
    xlabel('time');
    ylabel('gamma');
    title('Quadratic Well Operator: gamma - packet width');
    drawnow
    
    figure(5);
    surf(X, Y, V);
    shading interp;
    title('Quadratic Well: Surface of V');