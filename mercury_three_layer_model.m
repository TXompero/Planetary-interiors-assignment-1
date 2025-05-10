clearvars; clc; close all

%% Constant density model

R_planet = 2440e3;

rho_core = 7034.32;        
rho_mantle = 3343.35;     
rho_crust = 2903.03;

alpha = 0.8295;
beta = 0.984;

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
crust.rho_initial_guess = 2900;

crust.R1 = R_planet*beta;
crust.R2 = R_planet;
crust.n = 1e2;

% Compute profiles
planet = {core; mantle; crust};
[rho_planet,C_planet,r_vec,rho_vec,m_vec,g_vec,p_vec,T_vec] = solve_planet(planet);

figure(1)
subplot(1,3,1)
grid on
hold on
plot(m_vec,r_vec/1e3,'k--','LineWidth',1.5)
xlabel('M(r) [kg]')
ylabel('r [km]')
subplot(1,3,2)
grid on
hold on
plot(g_vec,r_vec/1e3,'k--','LineWidth',1.5)
xlabel('g(r) [m/s^2]')
ylabel('r [km]')
subplot(1,3,3)
grid on
hold on
plot(p_vec/1e9,r_vec/1e3,'k--','LineWidth',1.5)
xlabel('p(r) [GPa]')
ylabel('r [km]')

fprintf('bulk density: %f\nAdimensional moment of inertia: %f\n\n',rho_planet,C_planet)