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
for i = 1:7
    this_VGS = data_G_25_25(num*(i-1)+1:num*i, 2);
    this_IDS = data_G_25_25(num*(i-1)+1:num*i, 4);
    
    this_VDS = data_G_25_25(num*i, 1);
    this_VSB = data_G_25_25(num*i, 3);
    
    modeled_IDS = current_from_VGS(this_W, this_L, this_VGS,...
        this_VDS, this_VSB);
    
    plot(this_VGS, this_IDS*1e6);
    plot(this_VGS, modeled_IDS*1e6);
end

title('I_{DS} vs. V_{GS}');
xlabel('V_{GS} (V)');
ylabel('I_{DS} (\muA)');

% plot long-channel measured and modeled IDS vs VDS
num=37;
figure
hold on
for i = 1:1
    this_VDS = data_D_25_25(num*(i-1)+1:num*i, 1);
    this_IDS = data_D_25_25(num*(i-1)+1:num*i, 4);
    
    this_VGS = data_D_25_25(num*i, 2);
    this_VSB = data_D_25_25(num*i, 3);
    
    modeled_IDS = current_from_VDS(this_W, this_L, this_VGS,...
        this_VDS, this_VSB);
    
    plot(this_VDS, this_IDS*1e6);
    plot(this_VDS, modeled_IDS*1e6);
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

func_psi_sL = @(psi_sL_val) VGB - parameters.VFB -...
    parameters.gamma*sqrt(psi_sL_val + constants.phit*exp(...
    (psi_sL_val-2*parameters.phiF-VDB)/constants.phit)) - psi_sL_val;
psi_sL = fsolve(func_psi_sL, ones(size(VGB))*VDB);

% the difference psi_sL-psi_s0 is used often
delta_psi_s = psi_sL - psi_s0;

% alpha calculation - note that this is a DIFFERENT definition
% from the way alpha is defined in the book
alpha = 1 + parameters.gamma*(sqrt(psi_sL) - sqrt(psi_s0))./...
    (delta_psi_s);

% calculate drain current - components "due to" drift and diffusion
IDS1 = W/L*parameters.u*parameters.Cox * (VGB - parameters.VFB ...
    - psi_s0 - parameters.gamma*sqrt(psi_s0) -...
    alpha.*delta_psi_s/2) .* delta_psi_s;
IDS2 = W/L*parameters.u*parameters.Cox*constants.phit*alpha.*delta_psi_s;
IDS = IDS1 + IDS2;

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

% psi_sL is still a vector, because in this case VDB is a vector
% func_psi_sL = @(psi_sL_val) VGB - parameters.VFB -...
%     parameters.gamma*sqrt(psi_sL_val + constants.phit*exp(...
%     (psi_sL_val-2*parameters.phiF-VDB)/constants.phit)) - psi_sL_val;
% psi_sL = fsolve(func_psi_sL, ones(size(VDS))*max(VDS)*2);
func_psi_sL = @(psi_sL_val) (((VGB-parameters.VFB-psi_sL_val)...
    ./parameters.gamma).^2)...
    -psi_sL_val - constants.phit*exp(...
    (psi_sL_val-2*parameters.phiF-VDB)/constants.phit);
psi_sL = fsolve(func_psi_sL, VDB);

% the difference psi_sL-psi_s0 is used often
delta_psi_s = psi_sL - psi_s0;

% alpha calculation - note that this is a DIFFERENT definition
% from the way alpha is defined in the book
alpha = 1 + parameters.gamma*(sqrt(psi_sL) - sqrt(psi_s0))./...
    (delta_psi_s);

% calculate drain current - components "due to" drift and diffusion
IDS1 = W/L*parameters.u*parameters.Cox * (VGB - parameters.VFB...
    - psi_s0 - parameters.gamma*sqrt(psi_s0) -...
    alpha.*delta_psi_s/2) .* delta_psi_s;
IDS2 = W/L*parameters.u*parameters.Cox*constants.phit*alpha.*delta_psi_s;
IDS = IDS1 + IDS2;

end

%% Short-Channel Effects

% not sure if these will be separate functions or multipliers
% that are turned on/off above,
% so this is just a placeholder