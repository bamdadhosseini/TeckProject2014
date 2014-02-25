%% Script for generation of dt for various 
%% wind speed

domainsize = [2000, 2000];
meshsize = [100, 100];

% assume wind is stored in the structure 
% wind.vel, wind.dir, wind.time

dx = domainsize(1)/meshsize(1);
dy = domainsize(2)/meshsize(2);

% compute windspeed in x, y directions

velx = abs(wind.vel .* cos(wind.dir));
vely = abs(wind.vel .* sin(wind.dir));

nu = 0.9; % desirable courant number

numb_delta_t = 5; % number of delta_t intervals
max_delta_t = 200;
delta_t = min(min( nu*dx./velx, nu*dy./vely), max_delta_t );

delta_t_bnd = ...
    [0 10 30 50 100 200 Inf];
%[0, linspace(min(delta_t), max(delta_t), numb_delta_t), Inf];
wind.dt = 0;

for i = 2:numb_delta_t+1
   
    wind.dt = wind.dt + (delta_t > delta_t_bnd(i-1)).*...
        (delta_t <= delta_t_bnd(i)).*delta_t_bnd(i);
    
end

wind.dt = wind.dt + (delta_t > delta_t_bnd(end-1)).*...
        (delta_t < delta_t_bnd(end)).*delta_t_bnd(end-1);