% Update the positions, velocities and accelerations of the particles using
% the Velocity Verlet algorithm
%
% Parameters:
% num_particles: int
% num_dimensions: int
% positions: matrix with rows=num_dimensions, cols=num_particles
% velocities: matrix with rows=num_dimensions, cols=num_particles
% forces: matrix with rows=num_dimensions, cols=num_particles
% accelerations: matrix with rows=num_dimensions,cols=num_particles
% mass: int (the mass of each particles - assumed to be uniform)
% delta_t: int
%
% The Velocity Verlet Algorithm
%    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
%    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
%    a(t+dt) = f(t) / m
%
% Resource Used: https://people.sc.fsu.edu/~jburkardt/m_src/md/md.m

function [positions, velocities, accelerations] = ...
velocity_verlet(num_particles, num_dimensions, ...
    positions, velocities, forces, accelerations,...
    mass, delta_t)

r_mass = 1.0 / mass;

% compute the new positions
positions = posititions + ...
    velocities * delta_t + ...
    0.5 * accelerations * dt * dt;

% compute the new velocities
velocities = velocities + ...
    0.5 * (forces * r_mass * accelerations) * dt;

% compute the new accelerations
accelerations = forces * r_mass;

end