clc;
clear;
close all;

%% Define parameters
%% membrane properties
kappa = 25;  %% bending rigidity, Unit: k_BT
sigma = 0.012;
l_m = 4;
A_p = 100;    %% Unit: nm^2
mu = 5;   %% membrane-protein interaction, Unit: k_BT/nm^2
Chi_pp = -1*0.36;  %% protein-protein interaction, Unit: k_BT
%% Polymer properties
a = 10;      %% monomer size,  Unit: nm
Chi = 0;   %% Flory parameter
v = (1-2*Chi)*a^3;     %% Excluded volume parameter, Unit: nm^3
nu = 3/5;

% N_span = 20;     %% Total number of monomers
N_span = 0:1:30;
i = 0;
for N = N_span
    N
    i = i+1;
    N_m = N;
    c_p = sqrt(pi/6)*a*N^(nu)/4/kappa/A_p + 1/2/l_m;
    phi_span = 0.0000:0.0000001:1;
    j = 0;
    for phi_v = phi_span
        j = j+1;
        y_v(j) = -mu + A_p*kappa*c_p^2*phi_v + (1+N_m)*Chi_pp*phi_v + log(phi_v/(1-phi_v));
    end
    [yv,index_minv] = min(abs(y_v));
    phiv_min(i) = phi_span(index_minv);
    sigma_eff = sigma - kappa/2*c_p^2*(phi_span(index_minv))^2 + 1/A_p*(log(1-phi_span(index_minv)) - (1+N_m)/2*Chi_pp*(phi_span(index_minv))^2);
    gamma = kappa*c_p^2 + 1/(A_p*phi_span(index_minv)*(1-phi_span(index_minv))) + (1+N_m)*Chi_pp/A_p;
    kappa_eff = kappa*(1-kappa^2*c_p^2/gamma^2*(Chi_pp*(1+N_m)/kappa/A_p+c_p^2));
    R_t_an_eff(i) = sqrt(kappa_eff/(2*sigma_eff));
    Delta_phi_an_eff = kappa*c_p/(kappa*c_p^2+1/(A_p*phi_span(index_minv)*(1-phi_span(index_minv)))+(1+N_m)*Chi_pp/A_p)/R_t_an_eff(i);
    phi_t_an_eff(i) = phi_span(index_minv) + Delta_phi_an_eff;
    S_an_eff(i) = 1 + Delta_phi_an_eff/phi_span(index_minv);
end

plot(N_span,S_an_eff,'b','linewidth',2)
xlabel('\itN');
ylabel('\itS');