% Wave Packet Dynamics
%
% Resources Used: https://people.sc.fsu.edu/~jburkardt/m_src/md/md.m
function [positions, velocities, accelerations] = initialize(num_particles, num_dimensions)

    % THE RESOURCE LINKED ABOVE SEEDS AND GENERATES RANDOM POSITIONS -
    % should I implement that? We will need that when we jump from 1 to
    % more than 1 particle I assume? FOR NOW JUST SET IT AT A REASONABLE
    % PLACE SINCE JUST STARTING WITH 1 ELECTRON
    positions = zeros(num_particles_num_dimensions);
    
    velocities = zeros(num_particles, num_dimensions);
    accelerations = zeros(num_particles, num_dimensions);
end