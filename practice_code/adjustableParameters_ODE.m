% dydt = -ky + Fsin(omega*t)
%   We might want to solve this equation with different values of k, F, and
%   omega so we can explore how the system behaves
%

function adjustableParameters_ODE

% solve the ODE above from t=start_time through t=end_time, y(0) =
% initial_y

% define the parameter values
k = 10;
F = 1;
omega = 1;
time_range = [0, 20];
y0 = 0;

% NOTE: we have to specify that y is differentiated with respect to t with
% @(t,y)
[t_values, sol_values] = ode45(@(t,y) diff_eq(t,y,k,F,omega), time_range, y0);

plot(t_values, sol_values);

end

function dydt = diff_eq(t, y, k, F, omega)
    dydt = -k*y + F*sin(omega*t);
end