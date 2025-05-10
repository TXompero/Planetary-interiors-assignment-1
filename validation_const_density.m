clearvars; clc; close all

%% 2 layer planet

rho_core = 7245;        % in [kg/m^3]
rho_mantle = 3182;      % in [kg/m^3]

R_planet = 2439.7e3;    % in [m]
alpha = 0.8209;         % adim [-]

% Generate layers compatible with the function
core = struct;
mantle = struct;

% Define core
core.const_density = 1;

core.R1 = 0;
core.R2 = R_planet*alpha;
core.n = 5e3;
core.rho_initial_guess = rho_core;

% Define mantle
mantle.const_density = 1;
mantle.rho_initial_guess = rho_mantle;

mantle.R1 = R_planet*alpha;
mantle.R2 = R_planet;
mantle.n = 5e3;

% Define crust
crust.thermal_env.is_convective = 0;

crust.const_density = 1;
crust.rho_initial_guess = 2900;

crust.R1 = 2340e3;
crust.R2 = 2440e3;
crust.n = 1e2;

%% Compute profiles

planet = {core; mantle};

[rho_planet,C_planet,r_vec,rho_vec,m_vec,g_vec,p_vec] = solve_planet(planet);

% Analytical solution
rho_analytical = rho_core*alpha^3 + rho_mantle*(1-alpha^3);
C_analytical = 2*(rho_core*alpha^5/rho_analytical + rho_mantle*(1-alpha^5)/rho_analytical)/5;

% Plot solutions
fprintf('2 layer planet\n')
fprintf('Relative error in the bulk modulus:   %e \n',abs(rho_planet - rho_analytical)/rho_analytical)
fprintf('Relative error in the inertia moment: %e \n',abs(C_planet - C_analytical)/C_analytical)

%% 3 layer planet

rho_core = 7245;        % in [kg/m^3]
rho_mantle = 3182;      % in [kg/m^3]
rho_crust = 2700;       % in [kg/m^3]

R_planet = 2439.7e3;    % in [m]
alpha = 0.8209;         % adim [-]
beta = 0.993;           % adim [-]

% Generate layers compatible with the function
core = struct;
mantle = struct;
crust = struct;

% Define core
core.const_density = 1;

core.R1 = 0;
core.R2 = R_planet*alpha;
core.n = 5e3;
core.rho_initial_guess = rho_core;

% Define mantle
mantle.const_density = 1;
mantle.rho_initial_guess = rho_mantle;

mantle.R1 = R_planet*alpha;
mantle.R2 = R_planet*beta;
mantle.n = 5e3;

% Define crust
crust.const_density = 1;
crust.rho_initial_guess = rho_crust;

crust.R1 = R_planet*beta;
crust.R2 = R_planet;
crust.n = 1e2;

% Compute profiles
planet = {core; mantle; crust};

[rho_planet_3l,C_planet_3l,r_vec_3l,rho_vec_3l,m_vec_3l,g_vec_3l,p_vec_3l] = solve_planet(planet);

% Analytical solution
rho_analytical_3l = rho_core*alpha^3 + rho_mantle*(beta^3-alpha^3) + rho_crust*(1-beta^3);
C_analytical_3l = 2*(rho_core*alpha^5/rho_analytical_3l + rho_mantle*(beta^5-alpha^5)/rho_analytical_3l + rho_crust*(1-beta^5)/rho_analytical_3l)/5;

% Plot solutions
fprintf('\n3 layer planet\n')
fprintf('Relative error in the bulk modulus:   %e \n',abs(rho_planet_3l - rho_analytical_3l)/rho_analytical_3l)
fprintf('Relative error in the inertia moment: %e \n\n',abs(C_planet_3l - C_analytical_3l)/C_analytical_3l)