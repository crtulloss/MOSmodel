close all;
clear;
clc;

%Determine how long the model takes to generate the current from voltages
%input to it using the extracted parameters

W = 25e-4;
L = 25e-4;
temp = 300;
VGS = 2;
VDS = 1;
VSB = 0.5;

load('params.mat');
IDS = current(W, L, gamma, gamma_a, gamma_b, gamma_c, gamma_d,...
    VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c, a_theta_d,...
    eta_a, eta_b, eta_c, eta_d, delta_L, Ec, la, VE, temp, VGS, VDS, VSB);

%current plotting function
function IDS = current(W, L, gamma, gamma_a, gamma_b, gamma_c, gamma_d,...
    VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c, a_theta_d,...
    eta_a, eta_b, eta_c, eta_d, delta_L, Ec, la, VE, temp,VGS, VDS, VSB)

gamma = gamma*(gamma_a.*sqrt(VSB) + gamma_b/L^2 + gamma_c.*sqrt(VSB)...
    /L^2 + gamma_d);
a_theta = a_theta_a/sqrt(L) + a_theta_b/L + a_theta_c/L^2 + a_theta_d;
eta_E = eta_a/sqrt(L) + eta_b/L + eta_c/L^2 + eta_d;

L = L - delta_L;

% only matters for appendix K WI VDS-dependence test,
% where temperature is allowed to change
% otherwise phit = 0.026V
phit = constants.k * temp / constants.q;
phiF = phiF*phit/constants.phit;

VGB = VGS + VSB;
VDB = VDS + VSB;

% calculate drain and source surface potentials

% note that this works whether VSB is scalar or vector
% same for VDB-VSB below
func_psi_s0 = @(psi_s0_val) VGB - VFB -...
    gamma.*sqrt(psi_s0_val + phit*exp(...
    (psi_s0_val-2*phiF-VSB)/phit)) - psi_s0_val;
psi_s0 = fsolve(func_psi_s0, ones(size(VGB)).*VSB);

% the difference psi_sL - psi_s0 can be small, so instead
% of making that the difference between two solutions, we solve
% for it directly and us it to calculate psi_sL
func_delta_psi_s = @(delta_psi_s_val) -delta_psi_s_val...
    -gamma.*sqrt(delta_psi_s_val + psi_s0 +...
    phit*exp((delta_psi_s_val+psi_s0-...
    2*phiF-VDB)/phit)) +...
    gamma.*sqrt(psi_s0 + phit*exp(...
    (psi_s0-2*phiF-VSB)/phit));
delta_psi_s = fsolve(func_delta_psi_s, ones(size(VGB)).*(VDB-VSB));

% now calculate psi_sL
psi_sL = psi_s0 + delta_psi_s;

% alpha calculation - note that this is a DIFFERENT definition
% from the way alpha is defined in the book
alpha = 1 + gamma./(sqrt(psi_sL) + sqrt(psi_s0));

QB0 = -gamma.*parameters.Cox.*sqrt(psi_s0);
QBL = -gamma.*parameters.Cox.*sqrt(psi_sL);

QB_avg = (QB0 + QBL)/2;
QI_avg = -parameters.Cox.*(VGB - VFB - (psi_s0+psi_sL)/2) - QB_avg;

mu = mu0./(1-a_theta./constants.eps.*(QB_avg + eta_E*QI_avg));

% for calculating VDS_eff for vel sat and CLM
VW = (-gamma/2 + sqrt(gamma.^2 / 4 + (VGS + abs(VSB))...
    - VFB)).^2 - 2*phiF;
VDS_prime = VW - VSB;

A = 8;

VDSeff = VDS./((1 + (abs(VDS)./VDS_prime).^A).^(1/A));
VDSeff_CLM = VDS_prime.*((1 + (abs(VDS)./VDS_prime).^A).^(1/A));

vel_sat_factor = (1./(0.5 * (1 + sqrt(1 + 2*(VDSeff/(L*Ec)).^2))));
lp = la*log(1 + (VDSeff_CLM-VDS_prime)./VE);
CLM_factor = (1./(1 - lp./L));

% calculate drain current - components "due to" drift and diffusion
IDS1 = W/L.*mu.*parameters.Cox .* (VGB - VFB...
    - psi_s0 - gamma.*sqrt(psi_s0) -...
    alpha.*delta_psi_s/2) .* delta_psi_s;
IDS2 = W/L.*mu.*parameters.Cox.*phit.*alpha.*delta_psi_s;
IDS = CLM_factor.*vel_sat_factor.*(IDS1 + IDS2);
%IDS = IDS1 + IDS2;

end