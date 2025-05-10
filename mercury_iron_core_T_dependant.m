%% Test Adiabatic profile

clearvars; clc; close all

%% Define material used in the core

% Solid gamma-iron
Fe = struct;
Fe.ref_density = 8201.84;
Fe.ref_T = 2500;
Fe.ref_p = 20e9;
Fe.thermal_exp = 5.7953e-5;
Fe.K = 129.02;
Fe.c_p = 850;


%% Define layers

% Set i = 1 for min. temperature profile, 2 for average and 3 for maximum
% i = 0 or any other values to select own profile, see else statement below
% rho_core is treated as initial guess to obtain the temperature and pressure
% dependant density profile
i = 1;

if i == 1

    T_inner_core = 2200;

    rho_core = 7034;        
    rho_mantle = 4064.5;
    rho_crust = 3300;
    
    alpha = 0.6715;
    beta = 0.984;

elseif i == 2

    T_inner_core = 2375;

    rho_core = 7034;        
    rho_mantle = 4067;
    rho_crust = 3300;
    
    alpha = 0.6755;
    beta = 0.984;

elseif i == 3

    T_inner_core = 2550;

    rho_core = 7034;        
    rho_mantle = 4066;
    rho_crust = 3300;
    
    alpha = 0.68;
    beta = 0.984;

else

    %-------------------------------
    % Define here own model to test
    %-------------------------------

    T_inner_core = 2550;

    rho_core = 7034;        
    rho_mantle = 4066;
    rho_crust = 3300;
    
    alpha = 0.68;
    beta = 0.984;

end

R_planet = 2440e3;

core = struct;
mantle = struct;
crust = struct;

% Define core
core.material = Fe;

core.const_density = 0;
core.thermal_env.is_convective = 1;
core.thermal_env.T_lower = T_inner_core;

core.R1 = 0;
core.R2 = alpha*R_planet;
core.n = 2e3;
core.rho_initial_guess = 7500;

% Define mantle
mantle.thermal_env.is_convective = 0;

mantle.const_density = 1;
mantle.rho_initial_guess = rho_mantle;

mantle.R1 = alpha*R_planet;
mantle.R2 = beta*R_planet;
mantle.n = 1e3;

% Define crust
crust.thermal_env.is_convective = 0;

crust.const_density = 1;
crust.rho_initial_guess = rho_crust;

crust.R1 = beta*R_planet;
crust.R2 = R_planet;
crust.n = 1e2;

%% Compute profiles

planet = {core; mantle; crust};

[rho_planet,C_planet,r_vec,rho_vec,m_vec,g_vec,p_vec,T_vec] = solve_planet(planet);

%% Plot results

figure(1)
subplot(1,3,1)
grid on
hold on
plot(m_vec,r_vec/1e3,'k','LineWidth',1.5)
xlabel('M(r) [kg]')
ylabel('r [km]')
subplot(1,3,2)
grid on
hold on
plot(g_vec,r_vec/1e3,'k','LineWidth',1.5)
xlabel('g(r) [m/s^2]')
ylabel('r [km]')
subplot(1,3,3)
grid on
hold on
plot(p_vec/1e9,r_vec/1e3,'k','LineWidth',1.5)
xlabel('p(r) [GPa]')
ylabel('r [km]')

T_Fe = @(p) 1811 + 28.95*p - 0.32*p.^2 + 0.0013*p.^3;

figure(2)
plot(r_vec/1e3,T_vec,'k','LineWidth',1.5)
hold on
grid on
plot(r_vec/1e3,T_Fe(p_vec./1e9),'r','LineWidth',1.5)
legend('T','solidus line')
xlabel('r [km]')
ylabel('T(r) [K]')
xlim([0 R_planet*alpha/1e3])

figure(3)
plot(r_vec/1e3,rho_vec,'k','LineWidth',1.5)
xlabel('r [km]')
ylabel('\rho(r) [kg/m^3]')

fprintf('Inner core temperature: %f\nbulk density: %f\nAdimensional moment of inertia: %f\n\n',T_inner_core,rho_planet,C_planet)