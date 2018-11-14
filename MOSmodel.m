% Caleb Rees Tulloss
% Jagannaath Shiva Letchumanan
% ELEN 6302 MOS
% Project: Simplified All-Region MOSFET Model

% Main project code

clear;
clc;

%% Notes

% Parameters and physical constants are each in their own file

%% Long-Channel Model Construction
% We add a few effects at a time to show model development
close all

% OS-indpendent file separator and directory info
f = filesep;
vg_dir = 'vg_data';
vd_dir = 'vd_data';

% leakage current value (based on the long-channel measured data)
leakage_lim = 1.5e-11;

% name format: data_G/D_W_L
% file columns:  VDS	VGS     VSB     IDS

% long-channel data - used for extraction
data_G_25_25 = dlmread([vg_dir, f, 'W25000_L25000_idvg.txt']);
data_D_25_25 = dlmread([vd_dir, f, 'W25000_L25000_idvd.txt']);
this_W = 25e-4;
this_L = 25e-4;

% other datasets used for deltaL extraction
data_G = ["W25000_L25000_idvg.txt","W25000_L2000_idvg.txt","W25000_L1000_idvg.txt",...
    "W25000_L800_idvg.txt","W25000_L600_idvg.txt"];
data_G = strcat(vg_dir, f, data_G);
data_L = [25e-4;2e-4;1e-4;0.8e-4;0.6e-4];
% other datasets used for velocity sat Ec extraction
data_D = ["W25000_L600_idvd.txt","W25000_L800_idvd.txt",...
    "W25000_L1000_idvd.txt","W25000_L2000_idvd.txt",...
    "W25000_L25000_idvd.txt"];
data_D = strcat(vd_dir, f, data_D);
data_L_D = [0.6e-4;0.8e-4;1e-4;2e-4;25e-4];

% run parameter extraction
num_data_sets = 7;
num = 73;
[gamma, NA, VFB, phi0] = extract_gamma(data_G_25_25, num_data_sets, num);
phiF = constants.phit*log(NA/1e10);

% Extract the required mobility parameters
num=73;
num_data_sets = 7;

[mu0, a_theta, eta_E] = extract_mu(this_W, this_L, gamma, VFB, phiF,...
    constants.roomTemp, data_G_25_25, num_data_sets, num);

% Extract deltaL
delta_L = extract_deltaL(data_G,data_L);

num_data_sets = 7;
num = 73;
[gamma_a, gamma_b, gamma_c, gamma_d] = extract_gamma_short(data_G, data_L,...
    gamma, phi0, VFB+0.07, num_data_sets, num);

[a_theta_a, a_theta_b, a_theta_c, a_theta_d, eta_a, eta_b, eta_c, eta_d]...
    = extract_mu_short(data_G, data_L, this_W, gamma, VFB, phiF, mu0,...
    a_theta, eta_E, delta_L, constants.roomTemp, num_data_sets, num);

[Ec, la, VE] = extract_VDS_params(data_D, data_L_D, gamma,...
    gamma_a, gamma_b, gamma_c, gamma_d, VFB, phiF,...
    mu0, a_theta_a, a_theta_b, a_theta_c, a_theta_d,...
    eta_a, eta_b, eta_c, eta_d, delta_L);


%% Validation: Appendix K and more

% % all these tests for W=L=25um
% this_W = 25e-4;
% this_L = 25e-4;
%
% % rms error - remove leakage points bc we aren't doing leakage?
%
% % remember to add log plots back!
%
% % K1: DC tests
%
% % continuity and smooth behavior - densely spaced plots
% % IDS vs. VDS for fixed VGS
% cont_VGS = 1;
% cont_VSB = 0;
% cont_VDS = linspace(-0.1, 4, 1000)';
% cont_IDS = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
%     gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
%     a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L,...
%     constants.roomTemp,...
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
% cont_log_IDS = log(current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
%     gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
%     a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L,...
%     constants.roomTemp,...
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
% IDS_zeroBias_vgs = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
%     gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
%     a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L,...
%     constants.roomTemp,...
%     0, 0, 0);
% IDS_zeroBias_vds = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
%     gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
%     a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L,...
%     constants.roomTemp,...
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
%     WI_IDS_1 = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
%         gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
%         a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L,...
%         T,...
%         WI_VGS, WI_VDS_1, WI_VSB);
%     WI_IDS_2 = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
%         gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
%         a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L,...
%         T,...
%         WI_VGS, WI_VDS_2, WI_VSB);
%     % calculate dID/dVDS as a function of VGS
%     derivID = (WI_IDS_2 - WI_IDS_1)./delta_VDS;
%     % calculate SR
%     WI_SR = WI_IDS_1 ./ WI_VDS_1 ./ derivID;
%
%     plot(WI_VGS, WI_SR);
% end
% asymp = 2*(sqrt(exp(1)) - 1);
% plot(WI_VGS, asymp*ones(size(WI_VGS)), '--');
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
% sym_IDS = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
%     gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
%     a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L,...
%     constants.roomTemp,...
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
%     gm_IDS = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
%         gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
%         a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L,...
%         constants.roomTemp,...
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
%     g0_IDS = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
%         gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
%         a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L,...
%         constants.roomTemp,...
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


%% Plots for all Channel Lengths

close all

this_W = 25e-4;

num_vg = 73;
num_vd = 37;

num_datasets_vg = 9;
num_datasets_vd = 10;

% VG plots

vg_listing = dir(vg_dir);
vg_listing = vg_listing(~ismember({vg_listing.name}, {'.', '..'}));

rms_error_vgs = zeros(size(vg_listing, 1), num_datasets_vg);

for k = 1:size(vg_listing, 1)
    
    fname = vg_listing(k).name;
    
    splitname = strsplit(fname, '_');
    L_val = str2double(splitname{2}(2:end));
    this_L = L_val * 1e-7;
    
    data = dlmread([vg_dir, f, fname]);
    
    figure(4*k-3)
    hold on
    figure(4*k-2)
    hold on
    VSB_legend = cell(14, 1);
    for m = 1:7
        this_VGS = data(num_vg*(m-1)+1:num_vg*m, 2);
        this_IDS = data(num_vg*(m-1)+1:num_vg*m, 4);
    
        this_VDS = data(num_vg*m, 1);
        this_VSB = data(num_vg*m, 3);
        
        VSB_legend{m*2-1} = strcat('V_{SB} = ', num2str(this_VSB),...
            ', measured');
        VSB_legend{m*2} = strcat('V_{SB} = ', num2str(this_VSB),...
            ', modeled');
    
        modeled_IDS = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
            gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
            a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L, Ec,...
            la, VE, ...
            constants.roomTemp, this_VGS, this_VDS, this_VSB);
        
        ln_this_IDS = log(this_IDS);
        ln_modeled_IDS = log(modeled_IDS);
        
        figure(4*k-3)
        plot(this_VGS, this_IDS*1e6,'*');
        plot(this_VGS, modeled_IDS*1e6);
        figure(4*k-2)
        plot(this_VGS, ln_this_IDS,'*');
        plot(this_VGS, ln_modeled_IDS);

    
        modeled_IDS = modeled_IDS(15:end);
            this_IDS = this_IDS(15:end);
        
        % calculate rms error
        difference = this_IDS-modeled_IDS;
        normalized_difference = difference./this_IDS;
        sum_sq = sum(normalized_difference.^2);
        this_rms_error = sqrt(sum_sq/59);
    
        rms_error_vgs(k, m) = this_rms_error;
    end
    
    this_VDS = data(num_vg*m, 1);
    
    figure(4*k-3)
    title(['I_{DS} vs. V_{GS}, L = ', num2str(L_val/1000),...
        '\mum, V_{DS} = ', num2str(this_VDS), 'V']);
    xlabel('V_{GS} (V)');
    ylabel('I_{DS} (\muA)');
    legend(VSB_legend)

    figure(4*k-2)
    title(['ln(I_{DS}) vs. V_{GS}, L = ', num2str(L_val/1000),...
        '\mum, V_{DS} = ', num2str(this_VDS), 'V']);
    xlabel('V_{GS} (V)');
    ylabel('ln(I_{DS})');
    legend(VSB_legend)
    
    figure(4*k-1)
    hold on
    figure(4*k)
    hold on
    VSB_legend = cell(4, 1);
    for m = 8:9
        this_VGS = data(num_vg*(m-1)+1:num_vg*m, 2);
        this_IDS = data(num_vg*(m-1)+1:num_vg*m, 4);
    
        this_VDS = data(num_vg*m, 1);
        this_VSB = data(num_vg*m, 3);
        
        VSB_legend{(m-7)*2-1} =...
            strcat('V_{SB} = ', num2str(this_VSB),...
            ', measured');
        VSB_legend{(m-7)*2} =...
            strcat('V_{SB} = ', num2str(this_VSB),...
            ', modeled');
    
        modeled_IDS = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
            gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
            a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L, Ec,...
            la, VE, ...
            constants.roomTemp, this_VGS, this_VDS, this_VSB);
        
        ln_this_IDS = log(this_IDS);
        ln_modeled_IDS = log(modeled_IDS);
        
        figure(4*k-1)
        plot(this_VGS, this_IDS*1e6,'*');
        plot(this_VGS, modeled_IDS*1e6);
        figure(4*k)
        plot(this_VGS, ln_this_IDS,'*');
        plot(this_VGS, ln_modeled_IDS);

    
        modeled_IDS = modeled_IDS(15:end);
            this_IDS = this_IDS(15:end);
        
        % calculate rms error
        difference = this_IDS-modeled_IDS;
        normalized_difference = difference./this_IDS;
        sum_sq = sum(normalized_difference.^2);
        this_rms_error = sqrt(sum_sq/59);
    
        rms_error_vgs(k, m) = this_rms_error;
    end
    
    this_VDS = data(num_vg*m, 1);
    
    figure(4*k-1)
    title(['I_{DS} vs. V_{GS}, L = ', num2str(L_val/1000),...
        '\mum, V_{DS} = ', num2str(this_VDS), 'V']);
    xlabel('V_{GS} (V)');
    ylabel('I_{DS} (\muA)');
    legend(VSB_legend)

    figure(4*k)
    title(['ln(I_{DS}) vs. V_{GS}, L = ', num2str(L_val/1000),...
        '\mum, V_{DS} = ', num2str(this_VDS), 'V']);
    xlabel('V_{GS} (V)');
    ylabel('ln(I_{DS})');
    legend(VSB_legend)
end

% VD plots

vd_listing = dir(vd_dir);
vd_listing = vd_listing(~ismember({vd_listing.name}, {'.', '..'}));

rms_error_vds = zeros(size(vd_listing, 1), num_datasets_vd);

for k = 1:size(vd_listing, 1)
    
    fname = vd_listing(k).name;
    
    splitname = strsplit(fname, '_');
    L_val = str2double(splitname{2}(2:end));
    this_L = L_val * 1e-7;
    
    data = dlmread([vd_dir, f, fname]);
    
    figure(4*size(vg_listing, 1) + 4*k-3)
    hold on
    figure(4*size(vg_listing, 1) + 4*k-2)
    VGS_legend = cell(num_datasets_vd, 1);
    for m = 1:num_datasets_vd/2
        this_VDS = data(num_vd*(m-1)+1:num_vd*m, 1);
        this_IDS = data(num_vd*(m-1)+1:num_vd*m, 4);
    
        this_VGS = data(num_vd*m, 2);
        this_VSB = data(num_vd*m, 3);
        
        VGS_legend{m*2-1} = strcat('V_{GS} = ', num2str(this_VGS),...
            ', measured');
        VGS_legend{m*2} = strcat('V_{GS} = ', num2str(this_VGS),...
            ', modeled');
    
        modeled_IDS = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
            gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
            a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L, Ec,...
            la, VE, ...
            constants.roomTemp, this_VGS, this_VDS, this_VSB);
        
        figure(4*size(vg_listing, 1) + 4*k-3)
        plot(this_VDS, this_IDS*1e6,'*');
        plot(this_VDS, modeled_IDS*1e6);
        figure(4*size(vg_listing, 1) + 4*k-2)
        semilogy(this_VDS, this_IDS,'*');
        hold on
        semilogy(this_VDS, modeled_IDS);

        modeled_IDS_notbelow_leakage =...
            ((modeled_IDS < leakage_lim) == 0);
    
        modeled_IDS = modeled_IDS(modeled_IDS_notbelow_leakage);
            this_IDS = this_IDS(modeled_IDS_notbelow_leakage);
        
        % calculate rms error
        difference = this_IDS-modeled_IDS;
        normalized_difference = difference./this_IDS;
        sum_sq = sum(normalized_difference.^2);
        this_rms_error = sqrt(sum_sq/num);
    
        rms_error_vds(k, m) = this_rms_error;
    end
    
    this_VSB = data(num_vd*m, 3);
    
    figure(4*size(vg_listing, 1) + 4*k-3)
    title(['I_{DS} vs. V_{DS}, L = ', num2str(L_val/1000),...
        '\mum, V_{SB} = ', num2str(this_VSB), 'V']);
    xlabel('V_{DS} (V)');
    ylabel('I_{DS} (\muA)');
    legend(VGS_legend)

    figure(4*size(vg_listing, 1) + 4*k-2)
    title(['I_{DS} vs. V_{DS}, L = ', num2str(L_val/1000),...
        '\mum, V_{SB} = ', num2str(this_VSB), 'V']);
    xlabel('V_{DS} (V)');
    ylabel('I_{DS} (A)');
    legend(VGS_legend)
    
    figure(4*size(vg_listing, 1) + 4*k-1)
    hold on
    figure(4*size(vg_listing, 1) + 4*k)
    VGS_legend = cell(num_datasets_vd, 1);
    for m = num_datasets_vd/2+1:num_datasets_vd
        this_VDS = data(num_vd*(m-1)+1:num_vd*m, 1);
        this_IDS = data(num_vd*(m-1)+1:num_vd*m, 4);
    
        this_VGS = data(num_vd*m, 2);
        this_VSB = data(num_vd*m, 3);
        
        VGS_legend{(m-num_datasets_vd/2)*2-1} =...
            strcat('V_{GS} = ', num2str(this_VGS),...
            ', measured');
        VGS_legend{(m-num_datasets_vd/2)*2} =...
            strcat('V_{GS} = ', num2str(this_VGS),...
            ', modeled');
    
        modeled_IDS = current(this_W, this_L, gamma, gamma_a, gamma_b, gamma_c,...
            gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
            a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L, Ec,...
            la, VE, ...
            constants.roomTemp, this_VGS, this_VDS, this_VSB);
        
        
        figure(4*size(vg_listing, 1) + 4*k-1)
        plot(this_VDS, this_IDS*1e6,'*');
        plot(this_VDS, modeled_IDS*1e6);
        figure(4*size(vg_listing, 1) + 4*k)
        semilogy(this_VDS, this_IDS,'*');
        hold on
        semilogy(this_VDS, modeled_IDS);

        modeled_IDS_notbelow_leakage =...
            ((modeled_IDS < leakage_lim) == 0);
    
        modeled_IDS = modeled_IDS(modeled_IDS_notbelow_leakage);
            this_IDS = this_IDS(modeled_IDS_notbelow_leakage);
        
        % calculate rms error
        difference = this_IDS-modeled_IDS;
        normalized_difference = difference./this_IDS;
        sum_sq = sum(normalized_difference.^2);
        this_rms_error = sqrt(sum_sq/num);
    
        rms_error_vds(k, m) = this_rms_error;
    end
    
    this_VSB = data(num_vd*m, 3);
    
    figure(4*size(vg_listing, 1) + 4*k-1)
    title(['I_{DS} vs. V_{DS}, L = ', num2str(L_val/1000),...
        '\mum, V_{SB} = ', num2str(this_VSB), 'V']);
    xlabel('V_{DS} (V)');
    ylabel('I_{DS} (\muA)');
    legend(VGS_legend)

    figure(4*size(vg_listing, 1) + 4*k)
    title(['I_{DS} vs. V_{DS}, L = ', num2str(L_val/1000),...
        '\mum, V_{SB} = ', num2str(this_VSB), 'V']);
    xlabel('V_{DS} (V)');
    ylabel('I_{DS} (A)');
    legend(VGS_legend)
end

%% Current-Calculating Function

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

%% Parameter Extraction

%Extracting gamma from plots of Vt against sqrt(Vsb + phi0)

function [gamma, Na, VFB, phi0] = extract_gamma(data_G, num_data_sets, num)

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
VFB = VFB - 0.07;

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

%Extracting the value of delta L

function delta_L = extract_deltaL(datasets,data_L)

gm_max = zeros(numel(datasets),1);

% for loop determining the maximum transconductances for each channel
% length
for i=1:numel(datasets)
    data_G = dlmread(datasets(i));
    this_VGS = data_G(1:73, 2);
    this_IDS = data_G(1:73, 4);
    gm = diff(this_IDS)./diff(this_VGS);
    [gm_max(i) , ~] = max(gm);
end

%Linear Least Squares fit of 1/gmmax to the channel length
A = [data_L*1e4 ones(numel(datasets),1)];
B = 1./gm_max;
lms_result = ((A'*A)\A')*B;
slope = lms_result(1);
yint = lms_result(2);

%x-intercept is the delta_L
delta_L = -yint/slope*1e-4;

end

% Extract the variations in gamma due to charge sharing
function [a, b, c, d] = extract_gamma_short(datasets, data_L, gamma,...
    phi0, VFB, num_data_sets, num)

gamma_all = zeros(numel(datasets),num_data_sets);

for j=1:numel(datasets)
    data_G = dlmread(datasets(j));
    
    Vt = zeros(num_data_sets,1);
    Vsb = zeros(num_data_sets,1);
    
    for i=1:num_data_sets
        this_VGS = data_G(num*(i-1)+1:num*i, 2);
        this_IDS = data_G(num*(i-1)+1:num*i, 4);
        
        Vsb(i) = data_G(num*i, 3);
        %Using maximum transconductance to linearly extrapolate the threshold
        %voltage
        gm = diff(this_IDS)./diff(this_VGS);
        [gm_max , index_max_gm] = max(gm);
        %Assuming alpha = 1
        Vt(i) = this_VGS(index_max_gm) - this_IDS(index_max_gm)/gm_max - 0.05;
    end
    
    %Using the previously known values of VFB and phi0, gammas are obtained
    %for all VSB and L
    gamma_all(j,:) = (Vt-VFB-phi0)./sqrt(phi0+Vsb);
end

%Linear regression to fit gamma of the form:
%gamma_hat = gamma*(a*sqrt(VSB) + b/L^2 + c*sqrt(VSB)/L^2 + d)
Y = reshape(gamma_all,[],1)/gamma;
X1 = kron(sqrt(Vsb),ones(numel(data_L),1));
X2 = kron(ones(numel(Vsb),1),1./(data_L).^2);
X3 = X1.*X2;
X4 = ones(size(X1,1),1);
X = [X1 X2 X3 X4];
B = X\Y;

%Error of the fit coefficients
error1 = Y - X*B;
error1 = (error1./Y).^2*100;

%Coefficients of the linear regression
a = B(1);
b = B(2);
c = B(3);
d = B(4);

end

%Extracting the variations in a_theta and eta_E for short channels
function [a, b, c, d, a2, b2, c2, d2] = extract_mu_short(datasets,...
    data_L, W, gamma, VFB, phiF, mu0, a_theta_long, eta_E_long, delta_L,...
    temp, num_data_sets, num)

a_theta_L = zeros(numel(datasets),1);
a_theta = zeros(num_data_sets,1);
eta_L = zeros(numel(datasets),1);
eta = zeros(num_data_sets,1);

% only matters for appendix K WI VDS-dependence test,
% where temperature is allowed to change
% otherwise phit = 0.026V
phit = constants.k * temp / constants.q;
phiF = phiF*phit/constants.phit;

for j = 1:numel(datasets)
    data_G = dlmread(datasets(j));
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
        
        % Initial guesses based on long channel values
        x0 = [a_theta_long eta_E_long];
        
        %Function handle for fitting
        IDS_model = @(x, VGB) (W/(data_L(j)-delta_L)*mu0*parameters.Cox *...
            (VGB - VFB - psi_s0 - gamma*sqrt(psi_s0) - alpha.*...
            delta_psi_s/2).*delta_psi_s + W/(data_L(j)-delta_L)*mu0*...
            parameters.Cox*phit*alpha.*delta_psi_s)./(1-x(1)/...
            constants.eps.*(QB_avg + x(2)*QI_avg));
        
        %Least squares curve fit to determine the parameter values for each
        %VSB and L
        x = lsqcurvefit(IDS_model,x0,VGB,IDS);
        a_theta(i) = x(1);
        eta(i) = x(2);
    end
    
    %Parameters averaged over VSB to be just functions of L
    a_theta_L(j) = mean(a_theta);
    eta_L(j) = mean(eta);
end

%Linear regression for a_theta of the form:
%a_theta = a/sqrt(L) + b/L + c/L^2
Y = a_theta_L;
X1 = 1./sqrt(data_L);
X2 = 1./data_L;
X3 = 1./(data_L).^2;
X4 = ones(numel(data_L),1);
X = [X1 X2 X3 X4];
B = X\Y;

%Error of the fit coefficients
error1 = Y - X*B;
error1 = (error1./Y).^2*100;

%Coefficients of the linear regressions
a = B(1);
b = B(2);
c = B(3);
d = B(4);

%Linear regression for eta_E of the form:
%eta_E = a/sqrt(L) + b/L + c/L^2
Y = eta_L;
X1 = 1./sqrt(data_L);
X2 = 1./data_L;
X3 = 1./(data_L).^2;
X4 = ones(numel(data_L),1);
X = [X1 X2 X3 X4];
B = X\Y;

%Error of the fit coefficients
error2 = Y - X*B;
error2 = (error2./Y).^2*100;

%Coefficients of the linear regression
a2 = B(1);
b2 = B(2);
c2 = B(3);
d2 = B(4);

end

% Velocity saturation  and CLM - extracting Ec and lp

function [Ec, la, VE] = extract_VDS_params(datasets, data_L, gamma,...
    gamma_a, gamma_b, gamma_c, gamma_d, VFB, phiF,...
    mu0, a_theta_a, a_theta_b, a_theta_c, a_theta_d,...
    eta_a, eta_b, eta_c, eta_d, delta_L)

this_W = 25e-4;

num_vd = 37;

num_VE = 10;

VE_vec = linspace(0.0001, 1, num_VE)';

la_vec = zeros(num_VE, 1);
Ec_vec = zeros(num_VE, 1);
rms_error_vec = zeros(num_VE, 1);

for m = 1:num_VE
    this_VE = VE_vec(m);
    
    VDS_over_Ec = zeros(numel(datasets)*num_vd, 1);
    VDS_vec = zeros(numel(datasets)*num_vd, 1);
    
    x_sat_vec = zeros(numel(datasets)*num_vd, 1);
    lp_vec = zeros(numel(datasets)*num_vd, 1);
    
    for k=1:numel(datasets)
        this_L = data_L(k);
        L = this_L - delta_L;
        
        data = dlmread(datasets(k));
        this_VDS = data(149:185, 1);
        this_IDS = data(149:185, 4);
        this_VGS = data(149, 2);
        this_VSB = data(149, 3);
        
        % Ec
        
        modeled_current = current(this_W, this_L, gamma, gamma_a, gamma_b,...
            gamma_c, gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
            a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L, parameters.Ec,...
            parameters.la, parameters.VE, constants.roomTemp,...
            this_VGS, this_VDS, this_VSB);
        
        q = modeled_current ./ this_IDS;
        
        VDS_over_Ec(num_vd*(k-1) + 1:num_vd*k, 1) = L .*sqrt(...
            (((2.*q)-1).^2 - 1)/2);
        VDS_vec(num_vd*(k-1) + 1:num_vd*k, 1) = this_VDS;
        
        % la
        
        this_gamma = gamma*(gamma_a.*sqrt(this_VSB) + gamma_b/this_L^2 +...
            gamma_c.*sqrt(this_VSB)...
            /this_L^2 + gamma_d);
        
        % calculate saturation VDS and IDS
        VW = (-this_gamma/2 + sqrt(this_gamma.^2 / 4 + (this_VGS + this_VSB)...
            - VFB)).^2 - 2*phiF;
        VDS_prime = VW - this_VSB;
        
        IDS_prime = current(this_W, this_L, gamma, gamma_a, gamma_b,...
            gamma_c, gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
            a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L, parameters.Ec,...
            parameters.la, parameters.VE, constants.roomTemp,...
            this_VGS, VDS_prime, this_VSB);
        
        % only care about saturation region
        VDS_saturation_part = ((this_VDS > VDS_prime) == 1);
        this_VDS = this_VDS(VDS_saturation_part);
        this_IDS = this_IDS(VDS_saturation_part);
        
        num_sat = size(this_VDS,1);
        
        lp = L * (1 - IDS_prime ./ this_IDS);
        x = log(1 + (this_VDS-VDS_prime)./this_VE);
        
        x_sat_vec(num_vd*(k-1) + 1:num_vd*(k-1) + num_sat, 1) = x;
        lp_sat_vec(num_vd*(k-1) + 1:num_vd*(k-1) + num_sat, 1) = lp;
    end
    
    A = [VDS_vec ones(numel(datasets)*num_vd, 1)];
    B = VDS_over_Ec;
    lms_result = ((A'*A)\A')*B;
    slope = lms_result(1);
    yint = lms_result(2);
    this_Ec = 1/slope;
    this_Ec = real(this_Ec);
    
    Ec_vec(m) = this_Ec;
    
    nonzero_sat = ((x_sat_vec ~= 0) == 1);
    x_sat_vec = x_sat_vec(nonzero_sat);
    lp_sat_vec = lp_sat_vec(nonzero_sat);
    num_sat = size(lp_sat_vec, 1);
    
    A = [x_sat_vec ones(num_sat, 1)];
    B = lp_sat_vec;
    lms_result = ((A'*A)\A')*B;
    slope = lms_result(1);
    yint = lms_result(2);
    this_la = slope;
    
    la_vec(m) = this_la;
    
    sq_error_vec = zeros(numel(datasets)*num_vd, 1);
    % calculate error
    for k=1:numel(datasets)
        this_L = data_L(k);
        L = this_L - delta_L;
        
        data = dlmread(datasets(k));
        this_VDS = data(149:185, 1);
        this_IDS = data(149:185, 4);
        this_VGS = data(149, 2);
        this_VSB = data(149, 3);
        
        modeled_current = current(this_W, this_L, gamma, gamma_a, gamma_b,...
            gamma_c, gamma_d, VFB, phiF, mu0, a_theta_a, a_theta_b, a_theta_c,...
            a_theta_d, eta_a, eta_b, eta_c, eta_d, delta_L, this_Ec,...
            this_la, this_VE, constants.roomTemp,...
            this_VGS, this_VDS, this_VSB);
        
        sq_error = ((modeled_current - this_IDS) ./ this_IDS).^2;
        
        sq_error_vec(num_vd*(k-1) + 1:num_vd*k, 1) = sq_error;
    end
    rms_error = sqrt(sum(sq_error_vec)/(numel(datasets)*num_vd));
    rms_error_vec(m) = rms_error;
    
end

% find the VE with the min rms error
[~,index] = min(rms_error_vec);
Ec = Ec_vec(index);
la = la_vec(index);
VE = VE_vec(index);

end