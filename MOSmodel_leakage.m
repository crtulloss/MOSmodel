% Caleb Rees Tulloss
% Jagannaath Shiva Letchumanan
% ELEN 6302 MOS
% Project: Simplified All-Region MOSFET Model

% Main project code

clear;
clc;

%% Notes

% Parameters and physical constants are each in their own file

%% Plots
close all

% name format: data_G/D_W_L
% file columns:  VDS	VGS     VSB     IDS

% long-channel data
data_G_25_25 = dlmread('W25000_L25000_idvg.txt');
data_D_25_25 = dlmread('W25000_L25000_idvd.txt');
this_W = 25e-4;
this_L = 25e-4;

leakage_lim = 1.4e-11;

% plot long-channel measured and modeled IDS vs VGS
num=73;
figure
hold on

num_data_sets = 7;
% for rms error calculation, assuming all weights are 1
rms_error_vgs = zeros(num_data_sets, 1);
for i = 1:num_data_sets
    this_VGS = data_G_25_25(num*(i-1)+1:num*i, 2);
    this_IDS = data_G_25_25(num*(i-1)+1:num*i, 4);
    
    this_VDS = data_G_25_25(num*i, 1);
    this_VSB = data_G_25_25(num*i, 3);
    
    modeled_IDS = current_from_VGS(this_W, this_L, parameters.gamma,...
        parameters.VFB, parameters.phiF, parameters.u0, 0, 0,...
        constants.roomTemp, this_VGS, this_VDS, this_VSB);
    
    plot(this_VGS, this_IDS*1e6);
    plot(this_VGS, modeled_IDS*1e6,'*');
    
    modeled_IDS_notbelow_leakage =...
        ((modeled_IDS < leakage_lim) == 0);
    
    modeled_IDS = modeled_IDS(modeled_IDS_notbelow_leakage);
    this_IDS = this_IDS(modeled_IDS_notbelow_leakage);
    
    num_notbelow_leakage = length(this_IDS);
    
    % calculate rms error
    difference = this_IDS-modeled_IDS;
    normalized_difference = difference./this_IDS;
    sum_sq = sum(normalized_difference.^2);
    this_rms_error = sqrt(sum_sq/num_notbelow_leakage);
    
    rms_error_vgs(i) = this_rms_error;
end

title('I_{DS} vs. V_{GS}');
xlabel('V_{GS} (V)');
ylabel('I_{DS} (\muA)');

% plot long-channel measured and modeled IDS vs VDS
num=37;
figure
hold on

num_data_sets = 5;
% for rms error calculation, assuming all weights are 1
rms_error_vds = zeros(num_data_sets, 1);
for i = 1:num_data_sets
    this_VDS = data_D_25_25(num*(i-1)+1:num*i, 1);
    this_IDS = data_D_25_25(num*(i-1)+1:num*i, 4);
    
    this_VGS = data_D_25_25(num*i, 2);
    this_VSB = data_D_25_25(num*i, 3);
    
    modeled_IDS = current_from_VDS(this_W, this_L, parameters.gamma,...
        parameters.VFB, parameters.phiF, parameters.u0, 0, 0,...
        constants.roomTemp, this_VGS, this_VDS, this_VSB);
    
    plot(this_VDS, this_IDS*1e6);
    plot(this_VDS, modeled_IDS*1e6,'*');
    
    modeled_IDS_notbelow_leakage =...
        ((modeled_IDS < leakage_lim) == 0);
    
    modeled_IDS = modeled_IDS(modeled_IDS_notbelow_leakage);
    this_IDS = this_IDS(modeled_IDS_notbelow_leakage);
    
    num_notbelow_leakage = length(this_IDS);
    
    % calculate rms error
    difference = this_IDS-modeled_IDS;
    normalized_difference = difference./this_IDS;
    sum_sq = sum(normalized_difference.^2);
    this_rms_error = sqrt(sum_sq/num_notbelow_leakage);
    
    rms_error_vds(i) = this_rms_error;
end

title('I_{DS} vs. V_{DS}');
xlabel('V_{DS} (V)');
ylabel('I_{DS} (\muA)');

% run parameter extraction
num_data_sets = 7;
num = 73;
[gamma, NA, VFB] = extract_gamma(data_G_25_25, num_data_sets, num);
phiF = constants.phit*log(NA/1e10);

% re-plot using extracted parameters

% plot long-channel measured and modeled IDS vs VGS
num=73;
figure
hold on

num_data_sets = 7;
% for rms error calculation, assuming all weights are 1
rms_error_vgs2 = zeros(num_data_sets, 1);
for i = 1:num_data_sets
    this_VGS = data_G_25_25(num*(i-1)+1:num*i, 2);
    this_IDS = data_G_25_25(num*(i-1)+1:num*i, 4);
    
    this_VDS = data_G_25_25(num*i, 1);
    this_VSB = data_G_25_25(num*i, 3);
    
    modeled_IDS = current_from_VGS(this_W, this_L, gamma,...
        VFB, phiF, parameters.u0, 0, 0, constants.roomTemp,...
        this_VGS, this_VDS, this_VSB);
    
    plot(this_VGS, this_IDS*1e6);
    plot(this_VGS, modeled_IDS*1e6,'*');
    
    modeled_IDS_notbelow_leakage =...
        ((modeled_IDS < leakage_lim) == 0);
    
    modeled_IDS = modeled_IDS(modeled_IDS_notbelow_leakage);
    this_IDS = this_IDS(modeled_IDS_notbelow_leakage);
    
    num_notbelow_leakage = length(this_IDS);
    
    % calculate rms error
    difference = this_IDS-modeled_IDS;
    normalized_difference = difference./this_IDS;
    sum_sq = sum(normalized_difference.^2);
    this_rms_error = sqrt(sum_sq/num_notbelow_leakage);
    
    rms_error_vgs2(i) = this_rms_error;
end

title('I_{DS} vs. V_{GS}');
xlabel('V_{GS} (V)');
ylabel('I_{DS} (\muA)');

% plot long-channel measured and modeled IDS vs VDS
num=37;
figure
hold on

num_data_sets = 5;
% for rms error calculation, assuming all weights are 1
rms_error_vds2 = zeros(num_data_sets, 1);
for i = 1:num_data_sets
    this_VDS = data_D_25_25(num*(i-1)+1:num*i, 1);
    this_IDS = data_D_25_25(num*(i-1)+1:num*i, 4);
    
    this_VGS = data_D_25_25(num*i, 2);
    this_VSB = data_D_25_25(num*i, 3);
    
    modeled_IDS = current_from_VDS(this_W, this_L, gamma,...
        VFB, phiF, parameters.u0, 0, 0, constants.roomTemp,...
        this_VGS, this_VDS, this_VSB);
    
    plot(this_VDS, this_IDS*1e6);
    plot(this_VDS, modeled_IDS*1e6,'*');
    
    modeled_IDS_notbelow_leakage =...
        ((modeled_IDS < leakage_lim) == 0);
    
    modeled_IDS = modeled_IDS(modeled_IDS_notbelow_leakage);
    this_IDS = this_IDS(modeled_IDS_notbelow_leakage);
    
    num_notbelow_leakage = length(this_IDS);
    
    % calculate rms error
    difference = this_IDS-modeled_IDS;
    normalized_difference = difference./this_IDS;
    sum_sq = sum(normalized_difference.^2);
    this_rms_error = sqrt(sum_sq/num_notbelow_leakage);
    
    rms_error_vds2(i) = this_rms_error;
end

title('I_{DS} vs. V_{DS}');
xlabel('V_{DS} (V)');
ylabel('I_{DS} (\muA)');

%Extract the required mobility parameters
num=73;
num_data_sets = 7;

[mu0, a_theta, eta_E] = extract_mu(this_W, this_L, gamma, VFB, phiF,... 
    constants.roomTemp, data_G_25_25, num_data_sets, num);

figure

% for rms error calculation, assuming all weights are 1
rms_error_vgs3 = zeros(num_data_sets, 1);
for i = 1:num_data_sets
    this_VGS = data_G_25_25(num*(i-1)+1:num*i, 2);
    this_IDS = data_G_25_25(num*(i-1)+1:num*i, 4);
    
    this_VDS = data_G_25_25(num*i, 1);
    this_VSB = data_G_25_25(num*i, 3);
    
    modeled_IDS = current_from_VGS(this_W, this_L, gamma,...
        VFB, phiF, mu0, a_theta, eta_E, constants.roomTemp,...
        this_VGS, this_VDS, this_VSB);
    
    semilogy(this_VGS, this_IDS*1e6);
    hold on
    semilogy(this_VGS, modeled_IDS*1e6,'*');
    
    modeled_IDS_notbelow_leakage =...
        ((modeled_IDS < leakage_lim) == 0);
    
    modeled_IDS = modeled_IDS(modeled_IDS_notbelow_leakage);
    this_IDS = this_IDS(modeled_IDS_notbelow_leakage);
    
    num_notbelow_leakage = length(this_IDS);
    
    % calculate rms error
    difference = this_IDS-modeled_IDS;
    normalized_difference = difference./this_IDS;
    sum_sq = sum(normalized_difference.^2);
    this_rms_error = sqrt(sum_sq/num_notbelow_leakage);
    
    rms_error_vgs3(i) = this_rms_error;
end

title('I_{DS} vs. V_{GS}');
xlabel('V_{GS} (V)');
ylabel('I_{DS} (\muA)');

% plot long-channel measured and modeled IDS vs VDS
num=37;
figure
hold on

num_data_sets = 5;
% for rms error calculation, assuming all weights are 1
rms_error_vds3 = zeros(num_data_sets, 1);
for i = 1:num_data_sets
    this_VDS = data_D_25_25(num*(i-1)+1:num*i, 1);
    this_IDS = data_D_25_25(num*(i-1)+1:num*i, 4);
    
    this_VGS = data_D_25_25(num*i, 2);
    this_VSB = data_D_25_25(num*i, 3);
    
    modeled_IDS = current_from_VDS(this_W, this_L, gamma,...
        VFB, phiF, mu0, a_theta, eta_E, constants.roomTemp,...
        this_VGS, this_VDS, this_VSB);
    
    plot(this_VDS, this_IDS*1e6);
    plot(this_VDS, modeled_IDS*1e6,'*');
    
    modeled_IDS_notbelow_leakage =...
        ((modeled_IDS < leakage_lim) == 0);
    
    modeled_IDS = modeled_IDS(modeled_IDS_notbelow_leakage);
    this_IDS = this_IDS(modeled_IDS_notbelow_leakage);
    
    num_notbelow_leakage = length(this_IDS);
    
    % calculate rms error
    difference = this_IDS-modeled_IDS;
    normalized_difference = difference./this_IDS;
    sum_sq = sum(normalized_difference.^2);
    this_rms_error = sqrt(sum_sq/num_notbelow_leakage);
    
    rms_error_vds3(i) = this_rms_error;
end

title('I_{DS} vs. V_{DS}');
xlabel('V_{DS} (V)');
ylabel('I_{DS} (\muA)');

%% Validation: Appendix K and more

% all these tests for W=L=25um
this_W = 25e-4;
this_L = 25e-4;

% rms error - remove leakage points bc we aren't doing leakage?

% remember to add log plots back!

% % K1: DC tests
% 
% % continuity and smooth behavior - densely spaced plots
% % IDS vs. VDS for fixed VGS
% cont_VGS = 1;
% cont_VSB = 0;
% cont_VDS = linspace(-0.1, 4, 1000)';
% cont_IDS = current_from_VDS(this_W, this_L, gamma,...
%     VFB, phiF, mu0, a_theta, eta_E, constants.roomTemp,...
%     cont_VGS, cont_VDS, cont_VSB);
% figure
% plot(cont_VDS, cont_IDS*1e6);
% title('Continuity Test: I_{DS} vs. V_{DS}');
% xlabel('V_{DS} (V)');
% ylabel('I_{DS} (\muA)');
% % logID vs VGS for fixed VDS, and inspect plots for discont/kinks
% cont_VDS = 1;
% cont_VSB = 0;
% cont_VGS = linspace(0, 4, 1000)';
% cont_log_IDS = log(current_from_VGS(this_W, this_L, gamma,...
%     VFB, phiF, mu0, a_theta, eta_E, constants.roomTemp,...
%     cont_VGS, cont_VDS, cont_VSB));
% figure
% plot(cont_VGS, cont_log_IDS);
% title('Continuity Test: ln(I_{DS}) vs. V_{GS}');
% xlabel('V_{GS} (V)');
% ylabel('ln(I_{DS})');
% % inspect plots for discontinuities (problem with I),
% % and kinks (problem with derivative of I).
% % especially between regions of inversion (VGS plot),
% % around VDS=0, and between sat/nonsat (VDS plots)
% % when model is done, zoom in on these regions to check
% 
% % behavior at zero bias
% IDS_zeroBias_vgs = current_from_VGS(this_W, this_L, gamma,...
%     VFB, phiF, mu0, a_theta, eta_E, constants.roomTemp,...
%     0, 0, 0);
% IDS_zeroBias_vds = current_from_VDS(this_W, this_L, gamma,...
%     VFB, phiF, mu0, a_theta, eta_E, constants.roomTemp,...
%     0, 0, 0);
% % both should be 0!
% 
% % weak-inversion behavior
% % SR = (IDS/VDS)/(dID/dVDS) as a function of VGS
% % for small VDS and several temps
% T_sweep = 273+[-40 27 125]';
% WI_VGS = linspace(0, 4, 1000);
% WI_VSB = 0;
% delta_VDS = 0.00001;
% figure
% hold on
% for i = 1:3
%     T = T_sweep(i);
%     this_phit = constants.k * T / constants.q;
%     % fixed VDS, but need to find slope dID/dVDS,
%     % so pick two closely spaced VDS points
%     WI_VDS_1 = this_phit*0.5;
%     WI_VDS_2 = WI_VDS_1 + delta_VDS;
%     % calculate current at each of those point,
%     % as a function of VGS
%     WI_IDS_1 = current_from_VGS(this_W, this_L, gamma,...
%         VFB, phiF, mu0, a_theta, eta_E, T,...
%         WI_VGS, WI_VDS_1, WI_VSB);
%     WI_IDS_2 = current_from_VGS(this_W, this_L, gamma,...
%         VFB, phiF, mu0, a_theta, eta_E, T,...
%         WI_VGS, WI_VDS_2, WI_VSB);
%     % calculate dID/dVDS as a function of VGS
%     derivID = (WI_IDS_2 - WI_IDS_1)./delta_VDS;
%     % calculate SR
%     WI_SR = WI_IDS_1 ./ WI_VDS_1 ./ derivID;
%     
%     plot(WI_VGS, WI_SR);
% end
% plot(WI_VGS, 1.297*ones(size(WI_VGS)), '--');
% title('WI V_{DS} Dependence Test: S_R vs. V_{GS}');
% xlabel('V_{GS} (V)');
% ylabel('S_R');
% legend('T=-40\circC', 'T=27\circC', 'T=125\circC', '2(\surd\ite - 1)');
% % should be 2(sqrt(e)-1) = 1.297 in weak inversion, 
% % decrease toward 1 in strong inversion
% % check against figure K.1
% 
% % symmetry
% % voltages defined as in Fig. K.2,
% % except VB, which has the opposite polarity for consistency
% % with the rest of this model
% sym_VG = 1;
% sym_VB = -1;
% sym_VX = linspace(-0.2, 0.2, 1000)';
% sym_VD = sym_VX;
% sym_VS = -sym_VX;
% sym_VGS = sym_VG - sym_VS;
% sym_VDS = sym_VD - sym_VS;
% sym_VSB = sym_VS - sym_VB;
% sym_IDS = current_from_VDS(this_W, this_L, gamma,...
%     VFB, phiF, mu0, a_theta, eta_E, constants.roomTemp,...
%     sym_VGS, sym_VDS, sym_VSB);
% diff_VX = diff(sym_VX);
% sym_IDS_deriv1 = diff(sym_IDS)./diff_VX;
% sym_IDS_deriv2 = diff(sym_IDS_deriv1)./diff_VX(1:end-1);
% figure
% plot(sym_VX, sym_IDS);
% title('Symmetry Test: I_{X} vs. V_{X} (Figure K.2)');
% xlabel('V_{X} (V)');
% ylabel('I_{X} (A)');
% figure
% plot(sym_VX(1:end-1), sym_IDS_deriv1);
% title('Symmetry Test: dI_{X}/dV_{X} vs. V_{X} (Figure K.2)');
% xlabel('V_{X} (V)');
% ylabel('dI_{X}/dV_{X} (S)');
% figure
% plot(sym_VX(1:end-2), sym_IDS_deriv2);
% title('Symmetry Test: d^2I_{X}/dV_{X}^2 vs. V_{X} (Figure K.2)');
% xlabel('V_{X} (V)');
% ylabel('d^2I_{X}/dV_{X}^2 (S/V)');
% % check that ID is an odd function of Vx: VD or (-VS)
% % if not, how to measure (ID-IS)/2?
% % plot Ix, first and second derivatives as a function of Vx
% % this may be a problem for velocity saturation!
% 
% % K2: Conductance Tests
% 
% % weak and moderate inversion behavior
% % not for very low VGS because we are not modeling leakage,
% % and dividing by a very small IDS leads to unphysically high
% % transconductance efficiency
% gm_VGS = linspace(0.01, 4, 1000)';
% num_VSB = 7;
% gm_VSB = linspace(0, 3, num_VSB)';
% gm_VDS = 0.1;
% diff_gm_VGS = diff(gm_VGS);
% % plot gm/ID (transconductance efficiency) vs ID - semilogx
% % for several VSB, finely spaced VGS
% % from weak to strong inversion
% figure
% for i = 1:num_VSB
%     this_VSB = gm_VSB(i);
%     gm_IDS = current_from_VGS(this_W, this_L, gamma,...
%         VFB, phiF, mu0, a_theta, eta_E, constants.roomTemp,...
%         gm_VGS, gm_VDS, this_VSB);
%     
%     diff_gm_IDS = diff(gm_IDS);
%     gm = diff_gm_IDS./diff_gm_VGS;
%     gmeffic = gm./(gm_IDS(1:end-1));
%     
%     semilogx(gm_IDS(1:end-1), gmeffic);
%     hold on
% end
% semilogx(logspace(-11, -4, 10)', ones(10,1)/constants.phit, '--');
% title('Transconductance Test: g_m/I_{DS} vs. I_{DS}');
% xlabel('I_{DS} (A)');
% ylabel('g_m/I_{DS} (V^{-1})');
% axis([1e-11, 1e-4, 0, 40]);
% % curves should vary smoothly with VGS, and in weak inversion should
% % approach 1/phit as VSB increases
% % curves should peak and then decrease as ID decreases - not constant!
% 
% % output conductance
% num_VGS = 5;
% g0_VGS = linspace(2, 4, num_VGS)';
% g0_VSB = 0;
% g0_VDS = linspace(0, 4, 1000);
% diff_g0_VDS = diff(g0_VDS);
% figure
% hold on
% for i = 1:num_VGS
%     this_VGS = g0_VGS(i);
%     g0_IDS = current_from_VDS(this_W, this_L, gamma,...
%         VFB, phiF, mu0, a_theta, eta_E, constants.roomTemp,...
%         this_VGS, g0_VDS, g0_VSB);
%     
%     diff_g0_IDS = diff(g0_IDS);
%     g0 = diff_g0_IDS./diff_g0_VDS;
%     
%     plot(g0_VDS(1:end-1), g0*1e3);
% end
% title('Output Conductance Test: g_0 vs. V_{DS}');
% xlabel('V_{DS} (V)');
% ylabel('g_0 (mS)');
% % zoom in around sat transition
% % should be smooth - see Fig. K.5
% % may be problematic if we use any interpolation e.g. VDSeff


%% Functions: Long-Channel Model

function IDS = current_from_VGS(W, L, gamma, VFB, phiF, mu0, a_theta,...
    eta_E, temp, VGS, VDS, VSB)

% only matters for appendix K WI VDS-dependence test,
% where temperature is allowed to change
% otherwise phit = 0.026V
phit = constants.k * temp / constants.q;
phiF = phiF*phit/constants.phit;

VGB = VGS + VSB;
VDB = VDS + VSB;

% calculate drain and source surface potentials
func_psi_s0 = @(psi_s0_val) VGB - VFB -...
    gamma*sqrt(psi_s0_val + phit*exp(...
    (psi_s0_val-2*phiF-VSB)/phit)) - psi_s0_val;
psi_s0 = fsolve(func_psi_s0, ones(size(VGB))*VSB);

% the difference psi_sL - psi_s0 can be small, so instead
% of making that the difference between two solutions, we solve
% for it directly and us it to calculate psi_sL
func_delta_psi_s = @(delta_psi_s_val) -delta_psi_s_val...
    -gamma*sqrt(delta_psi_s_val + psi_s0 +...
    phit*exp((delta_psi_s_val+psi_s0-...
    2*phiF-VDB)/phit)) +...
    gamma*sqrt(psi_s0 + phit*exp(...
    (psi_s0-2*phiF-VSB)/phit));
delta_psi_s = fsolve(func_delta_psi_s, ones(size(VGB))*(VDB-VSB));

% now calculate psi_sL
psi_sL = psi_s0 + delta_psi_s;

% alpha calculation - note that this is a DIFFERENT definition
% from the way alpha is defined in the book
alpha = 1 + gamma./(sqrt(psi_sL) + sqrt(psi_s0));

QB0 = -gamma*parameters.Cox.*sqrt(psi_s0);
QBL = -gamma*parameters.Cox.*sqrt(psi_sL);
QB_avg = (QB0 + QBL)/2;

QI_avg = -parameters.Cox.*(VGB - VFB - (psi_s0+psi_sL)/2) - QB_avg;

mu = mu0./(1-a_theta/constants.eps.*(QB_avg + eta_E*QI_avg));

% calculate drain current
IDS1 = W/L*mu*parameters.Cox .* (VGB - VFB...
    - psi_s0 - gamma*sqrt(psi_s0) - alpha.*delta_psi_s/2)...
    .*delta_psi_s;
IDS2 = W/L*mu*parameters.Cox*phit.*alpha.*delta_psi_s;
IDS = IDS1+IDS2;

end

function IDS = current_from_VDS(W, L, gamma, VFB, phiF, mu0, a_theta, ...
    eta_E, temp, VGS, VDS, VSB)

% only matters for appendix K WI VDS-dependence test,
% where temperature is allowed to change
% otherwise phit = 0.026V
phit = constants.k * temp / constants.q;
phiF = phiF*phit/constants.phit;

VGB = VGS + VSB;
VDB = VDS + VSB;

% calculate drain and source surface potentials

% VSB can be a vector, so psi_s0 is still a vector
func_psi_s0 = @(psi_s0_val) VGB - VFB -...
    gamma*sqrt(psi_s0_val + phit*exp(...
    (psi_s0_val-2*phiF-VSB)/phit)) - psi_s0_val;
psi_s0 = fsolve(func_psi_s0, VSB);

% the difference psi_sL - psi_s0 can be small, so instead
% of making that the difference between two solutions, we solve
% for it directly and us it to calculate psi_sL
func_delta_psi_s = @(delta_psi_s_val) -delta_psi_s_val...
    -gamma*sqrt(delta_psi_s_val + psi_s0 +...
    phit*exp((delta_psi_s_val+psi_s0-...
    2*phiF-VDB)/phit)) +...
    gamma*sqrt(psi_s0 + phit*exp(...
    (psi_s0-2*phiF-VSB)/phit));
delta_psi_s = fsolve(func_delta_psi_s, VDB-VSB);

% now calculate psi_sL
psi_sL = psi_s0 + delta_psi_s;

% alpha calculation - note that this is a DIFFERENT definition
% from the way alpha is defined in the book
alpha = 1 + gamma./(sqrt(psi_sL) + sqrt(psi_s0));

QB0 = -gamma*parameters.Cox.*sqrt(psi_s0);
QBL = -gamma*parameters.Cox.*sqrt(psi_sL);
QB_avg = (QB0 + QBL)/2;

QI_avg = -parameters.Cox.*(VGB - VFB - (psi_s0+psi_sL)/2) - QB_avg;

mu = mu0./(1-a_theta/constants.eps.*(QB_avg + eta_E*QI_avg));

% calculate drain current - components "due to" drift and diffusion
IDS1 = W/L.*mu*parameters.Cox .* (VGB - VFB...
    - psi_s0 - gamma*sqrt(psi_s0) -...
    alpha.*delta_psi_s/2) .* delta_psi_s;
IDS2 = W/L.*mu*parameters.Cox*phit.*alpha.*delta_psi_s;
IDS = IDS1 + IDS2;

end

%Extracting gamma from plots of Vt against sqrt(Vsb + phi0)

function [gamma, Na, VFB] = extract_gamma(data_G, num_data_sets, num) 
    
index_max_gm = zeros(num_data_sets,1);
Vt = zeros(num_data_sets,1);
Vsb = zeros(num_data_sets,1);

for i=1:num_data_sets
    this_VGS = data_G(num*(i-1)+1:num*i, 2);
    this_IDS = data_G(num*(i-1)+1:num*i, 4);

    Vsb(i) = data_G(num*i, 3);

    %Using maximum transconductance to linearly extrapolate the threshold
    %voltage
    gm = diff(this_IDS)./diff(this_VGS);
    [gm_max , index_max_gm(i)] = max(gm);
    %Assuming alpha = 1
    Vt(i) = this_VGS(index_max_gm(i)) - this_IDS(index_max_gm(i))/gm_max - 0.05;
end

phi0_guess = 2*parameters.phiF + 5*constants.phit;
phi_0_value = phi0_guess-0.5:0.01:phi0_guess+0.5;
error_gamma = zeros(1,numel(phi_0_value));
    
%The linear least squares algorithm is used to obtain the best line fit
for i=1:numel(phi_0_value)
    A = [sqrt(phi_0_value(i)+Vsb) ones(num_data_sets,1)];
    P = A/(A'*A)*A';
    B = Vt;
    e = B - P*B;
    error_gamma(i) = norm(e);
end
    
%Linear fit with lowest error to a straight line gives the value of phi0
[~, index_phi] = min(error_gamma);
phi0 = phi_0_value(index_phi);

A = [sqrt(phi0+Vsb) ones(num_data_sets,1)];

%The result vector will contain two elements - the slope and the
%y-intercept of the best line fit
lms_result = ((A'*A)\A')*Vt;
gamma = lms_result(1);

Na = (gamma*parameters.Cox/sqrt(2*constants.q*constants.eps))^2;
VFB = Vt(1) - phi0 - gamma*sqrt(phi0);

end

function [mu0_avg, a_theta_avg, eta_E] = extract_mu(W, L, gamma, VFB, phiF,...
    temp, data_G, num_data_sets, num)

eta_E = 0.5;
mu0 = zeros(num_data_sets,1);
a_theta = zeros(num_data_sets,1);

% only matters for appendix K WI VDS-dependence test,
% where temperature is allowed to change
% otherwise phit = 0.026V
phit = constants.k * temp / constants.q;
phiF = phiF*phit/constants.phit;

for i = 1:num_data_sets
    VDS = data_G(num*i, 1);
    VSB = data_G(num*i, 3);
    
    Vm = VFB + 2*phiF + gamma*sqrt(2*phiF + VSB);
    index = find(data_G(:,2)>Vm+2,1);
    
    VGS = data_G(num*(i-1)+index:num*i, 2);
    IDS = data_G(num*(i-1)+index:num*i, 4);
    
    VGB = VGS + VSB;
    VDB = VDS + VSB;

    % calculate drain and source surface potentials
    func_psi_s0 = @(psi_s0_val) VGB - VFB -...
        gamma*sqrt(psi_s0_val + phit*exp(...
        (psi_s0_val-2*phiF-VSB)/phit)) - psi_s0_val;
    psi_s0 = fsolve(func_psi_s0, ones(size(VGB))*VSB);

    % the difference psi_sL - psi_s0 can be small, so instead
    % of making that the difference between two solutions, we solve
    % for it directly and us it to calculate psi_sL
    func_delta_psi_s = @(delta_psi_s_val) -delta_psi_s_val...
        -gamma*sqrt(delta_psi_s_val + psi_s0 +...
        phit*exp((delta_psi_s_val+psi_s0-...
        2*phiF-VDB)/phit)) +...
        gamma*sqrt(psi_s0 + phit*exp(...
        (psi_s0-2*phiF-VSB)/phit));
    delta_psi_s = fsolve(func_delta_psi_s, ones(size(VGB))*(VDB-VSB));

    % now calculate psi_sL
    psi_sL = psi_s0 + delta_psi_s;

    % alpha calculation - note that this is a DIFFERENT definition
    % from the way alpha is defined in the book
    alpha = 1 + gamma./(sqrt(psi_sL) + sqrt(psi_s0));
    
    QB0 = -gamma*parameters.Cox.*sqrt(psi_s0);
    QBL = -gamma*parameters.Cox.*sqrt(psi_sL);
    QB_avg = (QB0 + QBL)/2;

    QI_avg = -parameters.Cox.*(VGB - VFB - (psi_s0+psi_sL)/2) - QB_avg;

    % calculate drain current - components "due to" drift and diffusion
    % calculate drain current
    IDS1 = W/L*parameters.Cox * (VGB - VFB - psi_s0 - gamma*sqrt(psi_s0)...
        - alpha.*delta_psi_s/2).*delta_psi_s;
    IDS2 = W/L*parameters.Cox*phit*alpha.*delta_psi_s;
    IDS_model = IDS1+IDS2;

    B = IDS_model./IDS;
    A = [ones(numel(VGS),1) -QB_avg-eta_E*QI_avg];
    lms_result = (A'*A)\A'*B;
    mu0(i) = 1/lms_result(1);
    a_theta(i) = lms_result(2)*mu0(i)*constants.eps;
end

mu0_avg = mean(mu0);
a_theta_avg = mean(a_theta);

end

%% Short-Channel Effects

% not sure if these will be separate functions or multipliers
% that are turned on/off above,
% so this is just a placeholder