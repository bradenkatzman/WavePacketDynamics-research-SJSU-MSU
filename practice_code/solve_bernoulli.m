% solve the bernoulli differential equation dydt = (1/6)(t*(y^4)+ 2*y)
%   with the initial condition y(0) = -2

function solve_bernoulli
    initial_y = -2;
    time_range = [0,20];
    
    [t_values, solution_values] = ode45(@(t,y) diff_eq(t,y), ...
                                        time_range,...
                                        initial_y);
    plot(t_values, solution_values);                                
end

function dydt = diff_eq(t,y)
    dydt = (t*(y^4) + 2*y)/6;
end

% Notes:
%   The function handle @(t,y) diff_eq(t,y) tells MATLAB two things:
%       1. Which function to use to compute the value of dydt. It will use
%       this to compute dydt between t=0 and t=20
%       2. The @(t,y) tells MATLAB that the variable y is differentiated
%       with respect to t in the equation. These are dummy variables, they 
%       only need to match the symbols used in the diffential equation


% STOPPED AT 11.3
