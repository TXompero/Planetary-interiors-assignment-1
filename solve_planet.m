function [rho_planet,C_planet,r_vec,rho_vec,m_vec,g_vec,p_vec,T_vec] = solve_planet(layers)

G = 6.6743e-11;

n_idx_vec = ones(1,length(layers)+1);
n_tot = 1;

const_condition = 1;

r_vec = [];
rho_vec_new = [];
for i = 1:length(layers)

    r_vec = [r_vec, linspace(layers{i}.R1,layers{i}.R2,layers{i}.n)];
    rho_vec_new = [rho_vec_new, layers{i}.rho_initial_guess*ones(1,layers{i}.n)];

    n_tot = n_tot + layers{i}.n;
    n_idx_vec(i+1) = n_tot;

    if ~layers{i}.const_density
        const_condition = 0;
    end

end

it = 1;
err = 1e5;

while err > 100 && it < 100

    rho_vec = rho_vec_new;
    m_vec = zeros(1,length(r_vec));

    % Upward integration of mass
    for i = 1:length(r_vec)-1

        dMdR = 4*pi*rho_vec(i)*r_vec(i)^2;
        m_vec(i+1) = m_vec(i) + dMdR*(r_vec(i+1)-r_vec(i));

    end

    g_vec = G .* m_vec ./ (r_vec.^2);
    p_vec = zeros(1,length(r_vec));
    T_vec = zeros(1,length(r_vec));
    
    % Downward integration of pressure
    for i = length(r_vec):-1:2
    
        dpdR = rho_vec(i)*g_vec(i);
        p_vec(i-1) = p_vec(i) + dpdR*(r_vec(i)-r_vec(i-1));
    
    end

    % Break cycle if all layers have constant densities
    if const_condition
        break
    end
        
    for i = 1:length(layers)
    
        if ~layers{i}.thermal_env.is_convective
            continue
        end

        idx_active = n_idx_vec(i):n_idx_vec(i+1)-1;
    
        T_vec(n_idx_vec(i)) = layers{i}.thermal_env.T_lower;
        
        for j = n_idx_vec(i):n_idx_vec(i+1)-2
        
            dTdz = layers{i}.material.thermal_exp * g_vec(j) * T_vec(j) / layers{i}.material.c_p;

            % Avoid g = NaN, at center the core is ~isothermal
            if j == 1
                dTdz = 0;
            end

            T_vec(j+1) = T_vec(j) + dTdz * (r_vec(j) - r_vec(j+1));
        
        end
           
        rho_vec_new(idx_active) = layers{i}.material.ref_density * (1 - layers{i}.material.thermal_exp * ...
            (T_vec(idx_active) - layers{i}.material.ref_T) + (p_vec(idx_active) - layers{i}.material.ref_p)*1e-9/layers{i}.material.K);
    
    end

    err = max(abs(rho_vec - rho_vec_new));
    it = it + 1;

end

fprintf('Process converged after %d iterations\n',it)

rho_planet = trapz(r_vec,rho_vec.*r_vec.^2);
rho_planet = rho_planet*3/r_vec(end)^3;

C_planet = trapz(r_vec,rho_vec.*r_vec.^4);
C_planet = C_planet * 2 / (rho_planet*r_vec(end)^5);