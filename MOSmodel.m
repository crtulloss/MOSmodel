% Caleb Rees Tulloss
% Jagannaath Shiva Letchumanan
% ELEN 6302 MOS
% Project: Simplified All-Region MOSFET Model

% Main project code

%% Notes

% note to Shiva:
% for now I have just used some sample parameters from homework problems,
% to make sure my long-channel model works

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
        parameters.VFB, parameters.phiF, this_VGS, this_VDS, this_VSB);
    
    plot(this_VGS, this_IDS*1e6);
    plot(this_VGS, modeled_IDS*1e6,'*');
    
    % calculate rms error
    difference = this_IDS-modeled_IDS;
    normalized_difference = difference./this_IDS;
    sum_sq = sum(normalized_difference.^2);
    this_rms_error = sqrt(sum_sq/num);
    
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
        parameters.VFB, parameters.phiF, this_VGS, this_VDS, this_VSB);
    
    plot(this_VDS, this_IDS*1e6);
    plot(this_VDS, modeled_IDS*1e6,'*');
    
    % calculate rms error
    difference = this_IDS-modeled_IDS;
    normalized_difference = difference./this_IDS;
    sum_sq = sum(normalized_difference.^2);
    this_rms_error = sqrt(sum_sq/num);
    
    rms_error_vds(i) = this_rms_error;
end

title('I_{DS} vs. V_{DS}');
xlabel('V_{DS} (V)');
ylabel('I_{DS} (\muA)');

num_data_sets = 7;
num = 73;
[gamma, NA, VFB] = extract_gamma(data_G_25_25, num_data_sets, num);
phiF = constants.phit*log(NA/1e10);

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
        VFB, phiF, this_VGS, this_VDS, this_VSB);
    
    plot(this_VGS, this_IDS*1e6);
    plot(this_VGS, modeled_IDS*1e6,'*');
    
    % calculate rms error
    difference = this_IDS-modeled_IDS;
    normalized_difference = difference./this_IDS;
    sum_sq = sum(normalized_difference.^2);
    this_rms_error = sqrt(sum_sq/num);
    
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
        VFB, phiF, this_VGS, this_VDS, this_VSB);
    
    plot(this_VDS, this_IDS*1e6);
    plot(this_VDS, modeled_IDS*1e6,'*');
    
    % calculate rms error
    difference = this_IDS-modeled_IDS;
    normalized_difference = difference./this_IDS;
    sum_sq = sum(normalized_difference.^2);
    this_rms_error = sqrt(sum_sq/num);
    
    rms_error_vds2(i) = this_rms_error;
end

title('I_{DS} vs. V_{DS}');
xlabel('V_{DS} (V)');
ylabel('I_{DS} (\muA)');

%% Functions: Long-Channel Model

function IDS = current_from_VGS(W, L, gamma, VFB, phiF, VGS, VDS, VSB)

VGB = VGS + VSB;
VDB = VDS + VSB;

% calculate drain and source surface potentials
func_psi_s0 = @(psi_s0_val) VGB - VFB -...
    gamma*sqrt(psi_s0_val + constants.phit*exp(...
    (psi_s0_val-2*phiF-VSB)/constants.phit)) - psi_s0_val;
psi_s0 = fsolve(func_psi_s0, ones(size(VGB))*VSB);

% the difference psi_sL - psi_s0 can be small, so instead
% of making that the difference between two solutions, we solve
% for it directly and us it to calculate psi_sL
func_delta_psi_s = @(delta_psi_s_val) -delta_psi_s_val...
    -gamma*sqrt(delta_psi_s_val + psi_s0 +...
    constants.phit*exp((delta_psi_s_val+psi_s0-...
    2*phiF-VDB)/constants.phit)) +...
    gamma*sqrt(psi_s0 + constants.phit*exp(...
    (psi_s0-2*phiF-VSB)/constants.phit));
delta_psi_s = fsolve(func_delta_psi_s, ones(size(VGB))*(VDB-VSB));

% now calculate psi_sL
psi_sL = psi_s0 + delta_psi_s;

% alpha calculation - note that this is a DIFFERENT definition
% from the way alpha is defined in the book
alpha = 1 + gamma./(sqrt(psi_sL) + sqrt(psi_s0));

% calculate drain current
IDS1 = W/L*parameters.u0*parameters.Cox * (VGB - VFB...
    - psi_s0 - gamma*sqrt(psi_s0) - alpha.*delta_psi_s/2)...
    .*delta_psi_s;
IDS2 = W/L*parameters.u0*parameters.Cox*constants.phit*alpha.*delta_psi_s;
IDS = IDS1+IDS2;

end

function IDS = current_from_VDS(W, L, gamma, VFB, phiF, VGS, VDS, VSB)

VGB = VGS + VSB;
VDB = VDS + VSB;

% calculate drain and source surface potentials

% note that in this case psi_s0 is a constant since
% VGB and VSB are both fixed
func_psi_s0 = @(psi_s0_val) VGB - VFB -...
    gamma*sqrt(psi_s0_val + constants.phit*exp(...
    (psi_s0_val-2*phiF-VSB)/constants.phit)) - psi_s0_val;
psi_s0 = fsolve(func_psi_s0, VSB);

% the difference psi_sL - psi_s0 can be small, so instead
% of making that the difference between two solutions, we solve
% for it directly and us it to calculate psi_sL
func_delta_psi_s = @(delta_psi_s_val) -delta_psi_s_val...
    -gamma*sqrt(delta_psi_s_val + psi_s0 +...
    constants.phit*exp((delta_psi_s_val+psi_s0-...
    2*phiF-VDB)/constants.phit)) +...
    gamma*sqrt(psi_s0 + constants.phit*exp(...
    (psi_s0-2*phiF-VSB)/constants.phit));
delta_psi_s = fsolve(func_delta_psi_s, VDB-VSB);

% now calculate psi_sL
psi_sL = psi_s0 + delta_psi_s;

% alpha calculation - note that this is a DIFFERENT definition
% from the way alpha is defined in the book
alpha = 1 + gamma./(sqrt(psi_sL) + sqrt(psi_s0));

% calculate drain current - components "due to" drift and diffusion
IDS1 = W/L*parameters.u0*parameters.Cox * (VGB - VFB...
    - psi_s0 - gamma*sqrt(psi_s0) -...
    alpha.*delta_psi_s/2) .* delta_psi_s;
IDS2 = W/L*parameters.u0*parameters.Cox*constants.phit*alpha.*delta_psi_s;
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

Na = (gamma*parameters.Cox/constants.sqrt2qeps)^2;
VFB = Vt(1) - phi0 - gamma*sqrt(phi0);

end

%% Short-Channel Effects

% not sure if these will be separate functions or multipliers
% that are turned on/off above,
% so this is just a placeholder