% Simple first order differential equation
% dydt = -10y + sin(t)

% evaluate the simple ode at (0,0)
z1 = simple_ode(0, 0); % this equation computes the time derivative

% using ode45(@(variables)function_name(variables), 
%             [start_time end_time],
%             intial_value_y)
%
% Return value:
% each row in the solution array y corresponds to a value return in column
% vector t
%
% *** Note: The 'function handle', the @, is used to tell MATLAB the name
% of the function to be integrated
[t_values, y_values] = ode45(@(t,y) simple_ode(t,y), ...
                            [0, 20], ...
                            0);
                            % param 1 - the ODE, and the variables it
                            % requires. And the variable 'y' and what is
                            % differentiated with respect to, 't'
                            % 
                            % param 2 - the range of values for the
                            % variable that the derivative is wrt.
                            % Here, wrt time. I.e. the limits of
                            % integration
                            %
                            % param 3 - the starting value of y
 
% plot the values. Graphy of y as a function of t
plot(t_values, y_values);
                            
