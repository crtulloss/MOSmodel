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
    
    modeled_IDS = current_from_VGS(this_W, this_L, this_VGS,...
        this_VDS, this_VSB);
    
    plot(this_VGS, this_IDS*1e6);
    plot(this_VGS, modeled_IDS*1e6);
    
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
    
    modeled_IDS = current_from_VDS(this_W, this_L, this_VGS,...
        this_VDS, this_VSB);
    
    plot(this_VDS, this_IDS*1e6);
    plot(this_VDS, modeled_IDS*1e6);
    
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

%% Functions: Long-Channel Model

function IDS = current_from_VGS(W, L, VGS, VDS, VSB)

VGB = VGS + VSB;
VDB = VDS + VSB;

% calculate drain and source surface potentials
func_psi_s0 = @(psi_s0_val) VGB - parameters.VFB -...
    parameters.gamma*sqrt(psi_s0_val + constants.phit*exp(...
    (psi_s0_val-2*parameters.phiF-VSB)/constants.phit)) - psi_s0_val;
psi_s0 = fsolve(func_psi_s0, ones(size(VGB))*VSB);

% the difference psi_sL - psi_s0 can be small, so instead
% of making that the difference between two solutions, we solve
% for it directly and us it to calculate psi_sL
func_delta_psi_s = @(delta_psi_s_val) -delta_psi_s_val...
    -parameters.gamma*sqrt(delta_psi_s_val + psi_s0 +...
    constants.phit*exp((delta_psi_s_val+psi_s0-...
    2*parameters.phiF-VDB)/constants.phit)) +...
    parameters.gamma*sqrt(psi_s0 + constants.phit*exp(...
    (psi_s0-2*parameters.phiF-VSB)/constants.phit));
delta_psi_s = fsolve(func_delta_psi_s, ones(size(VGB))*(VDB-VSB));

% now calculate psi_sL
psi_sL = psi_s0 + delta_psi_s;

% calculate drain current
IDS1 = W/L*parameters.u*parameters.Cox * (VGB - parameters.VFB...
    - (psi_sL+psi_s0)/2 - parameters.gamma/2*(sqrt(psi_sL)+sqrt(psi_s0)))...
    .*delta_psi_s;
IDS2 = W/L*parameters.u*parameters.Cox*constants.phit*(...
    delta_psi_s + parameters.gamma*(sqrt(psi_sL) - sqrt(psi_s0)));
IDS = IDS1+IDS2;
end

function IDS = current_from_VDS(W, L, VGS, VDS, VSB)

VGB = VGS + VSB;
VDB = VDS + VSB;

% calculate drain and source surface potentials

% note that in this case psi_s0 is a constant since
% VGB and VSB are both fixed
func_psi_s0 = @(psi_s0_val) VGB - parameters.VFB -...
    parameters.gamma*sqrt(psi_s0_val + constants.phit*exp(...
    (psi_s0_val-2*parameters.phiF-VSB)/constants.phit)) - psi_s0_val;
psi_s0 = fsolve(func_psi_s0, VSB);

% the difference psi_sL - psi_s0 can be small, so instead
% of making that the difference between two solutions, we solve
% for it directly and us it to calculate psi_sL
func_delta_psi_s = @(delta_psi_s_val) -delta_psi_s_val...
    -parameters.gamma*sqrt(delta_psi_s_val + psi_s0 +...
    constants.phit*exp((delta_psi_s_val+psi_s0-...
    2*parameters.phiF-VDB)/constants.phit)) +...
    parameters.gamma*sqrt(psi_s0 + constants.phit*exp(...
    (psi_s0-2*parameters.phiF-VSB)/constants.phit));
delta_psi_s = fsolve(func_delta_psi_s, VDB-VSB);

% now calculate psi_sL
psi_sL = psi_s0 + delta_psi_s;

% calculate drain current
IDS1 = W/L*parameters.u*parameters.Cox * (VGB - parameters.VFB...
    - (psi_sL+psi_s0)/2 - parameters.gamma/2*(sqrt(psi_sL)+sqrt(psi_s0)))...
    .*delta_psi_s;
IDS2 = W/L*parameters.u*parameters.Cox*constants.phit*(...
    delta_psi_s + parameters.gamma*(sqrt(psi_sL) - sqrt(psi_s0)));
IDS = IDS1+IDS2;

end

%% Short-Channel Effects

% not sure if these will be separate functions or multipliers
% that are turned on/off above,
% so this is just a placeholder