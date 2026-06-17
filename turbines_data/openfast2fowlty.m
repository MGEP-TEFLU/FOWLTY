%% openfast2fowlty.m
% Generates a FOWLTY-compatible turbine data file (.mat) from OpenFAST input files.
%
% Reads turbine geometry, structural properties, and controller parameters
% from an OpenFAST model directory and builds the wt, env, and pub structs
% used by FOWLTY/SimWindFarm.
%
% If a pre-computed Cp/Ct table is found in the DISCON.IN file, it is used
% directly. Otherwise, a built-in BEM solver is run as a fallback.
%
% AUTHORS: Based on original NREL 5MW script by Jacob Grunnet and Mikael Svenstrup.
%          Extended for general turbines by Yerai Peña-Sanchez.

% =========================================================================
%  USER INPUT SECTION - Edit the variables below
% =========================================================================

% --- Path to the main OpenFAST .fst file ---
% fst_path = 'IEA-15-240-RWT_git/OpenFAST/IEA-15-240-RWT-Monopile/IEA-15-240-RWT-Monopile.fst';
fst_path = 'IEA-3.4-130-RWT-git/openfast/IEA-3.4-130-RWT.fst';

% --- Output name (a file <output_name>.mat will be saved) ---
output_name = 'IEA3p4MW';

% --- Platform options ---
% 'none'   : fixed-bottom / no platform dynamics (outputs zero motion)
% 'scaled' : scale DeepCWind (NREL 5MW) platform to match this turbine
% 'custom' : provide your own platform data in the plat struct below
platform_mode = 'scaled';

% --- Save output .mat file? ---
save_output = true;

% --- Load data into Simulink library block after running? ---
% Set to true to interactively save data into a Simulink block.
% You will be prompted to select the block and press Enter.
interactive = true;

% --- Pitch controller gains [-] ---
% ROSCO gains are tuned for high-fidelity aeroelastic models (OpenFAST)
% and are not used directly here. Instead, fixed gains Kp and Ki are used
% with a flat gain scheduling table (gs = 1), which is more robust in the
% simplified FOWLTY model. Tune these values for your turbine:
%   Kp : proportional gain. More negative = faster response.
%   Ki : integral gain.     More negative = stronger steady-state correction.
%
% Validated values per turbine:
%   IEA 15MW, floating (scaled DeepCWind) : Kp = -0.5,  Ki = -0.05
%   IEA 3.4MW onshore, fixed-bottom : Kp = -0.01, Ki = -0.001
%
% If the pitch oscillates at above-rated winds, reduce the gains.
% If the rotor speed takes too long to settle, increase them.
pitch_Kp = -0.01;   % [-] proportional gain
pitch_Ki = -0.001;  % [-] integral gain

% --- Pitch controller gain scaling factor [-] ---
% Additional scaling applied on top of Kp and Ki. Useful for quick tuning
% without changing the base values above.
%   1.0 : no additional scaling
%   0.5 : halve the gains (more conservative)
pitch_gain_scale = 1.0;   % [-]

% --- Pitch controller design parameters ---
% Only used if no ROSCO table is found (fallback computation from scratch).
% Typical values: damp=0.7, wund=0.6 (same as NREL 5MW baseline)
ctrl_damp = 0.7;    % Damping ratio [-]
ctrl_wund = 0.6;    % Undamped natural frequency [rad/s]

% --- BEM grid (only used if no Cp/Ct file is found) ---
% If a pre-computed table exists, these are ignored.
bem_beta = -5:1:30;     % Pitch angles [deg]
bem_tsr  = 2:0.5:14.5;  % Tip speed ratios [-]

% --- Rated power [W] ---
% This is used for the pub struct. Should match VS_RtPwr in DISCON.IN.
% Leave as NaN to read automatically from DISCON.IN.
P_rated_user = NaN;   % e.g. 15e6 for 15 MW

% --- Tower first fore-aft natural frequency [rad/s] ---
% Cannot be computed reliably from ElastoDyn files alone (monopile soil
% interaction and substructure effects make Rayleigh method inaccurate).
% Set the known value for your turbine here. References:
%   IEA 15MW monopile : 1.55 rad/s (0.247 Hz)  - Gaertner et al. 2020 / NREL BModes
%   IEA 3.4MW onshore : 2.2 rad/s (0.350 Hz)  - estimated from tower geometry (108m, 3-6m diameter)
%   NREL 5MW          : 2.01 rad/s (0.320 Hz)  - Jonkman et al. 2009
tower_eigfreq = 2.2;   % [rad/s] - EDIT THIS FOR YOUR TURBINE

% --- Air density [kg/m³] ---
rho = 1.225;

% =========================================================================
%  END OF USER INPUT - Do not edit below this line
% =========================================================================
clearvars -except fst_path output_name rho ctrl_damp ctrl_wund ...
    bem_beta bem_tsr plat P_rated_user platform_mode deepcwind_path ...
    platform_DoF interactive tower_eigfreq save_output ...
    pitch_Kp pitch_Ki pitch_gain_scale
clc

fprintf('=============================================================\n');
fprintf(' openfast2fowlty - Turbine data file generator for FOWLTY\n');
fprintf('=============================================================\n\n');

%% --- Step 1: Parse OpenFAST files ---
fprintf('[1/5] Parsing OpenFAST input files...\n');

files = parse_fst(fst_path);
ed    = parse_elastodyn(files.elastodyn);
sd    = parse_servodyn(files.servodyn);
ad    = parse_aerodyn(files.aerodyn);
blade = parse_aerodyn_blade(ad.blade_file);
dc    = parse_discon(sd.discon_file);
bl    = parse_elastodyn_blade(ed.blade_file, ed.HubRad, ed.TipRad);
tw    = parse_elastodyn_tower(ed.tower_file, ed.TowerHt);

fprintf('\n');

%% --- Step 2: Aerodynamic tables (Cp/Ct) ---
fprintf('[2/5] Loading aerodynamic performance tables...\n');

perf_file_exists = ~isempty(dc.perf_file) && exist(dc.perf_file, 'file');

if perf_file_exists
    fprintf('  Cp/Ct table found: %s\n', dc.perf_file);
    cpct = parse_cpct_ref(dc.perf_file);
    % Convert to FOWLTY convention: cp(i_beta, i_tsr), ct(i_beta, i_tsr)
    % parse_cpct_ref returns (n_tsr x n_pitch), so transpose
    cp   = cpct.cp';   % (n_pitch x n_tsr)
    ct   = cpct.ct';
    beta = cpct.beta;  % [deg]
    tsr  = cpct.tsr;
else
    fprintf('\n');
    fprintf('  !! WARNING: No Cp/Ct performance file detected.\n');
    fprintf('  !! Running built-in simplified BEM solver.\n');
    fprintf('  !! For higher accuracy, consider using CCBlade (WISDEM):\n');
    fprintf('  !!   pip install git+https://github.com/WISDEM/CCBlade.git\n');
    fprintf('  !!   See: https://github.com/WISDEM/CCBlade\n');
    fprintf('\n');

    % Load and extrapolate airfoils with Viterna method
    fprintf('  Loading %d airfoil polars...\n', ad.num_airfoils);
    blade_length = blade.span(end) - blade.span(1);
    AR           = blade_length / mean(blade.chord);
    airfoils     = cell(1, ad.num_airfoils);
    for i = 1:ad.num_airfoils
        af_raw      = parse_airfoil(ad.airfoil_files{i});
        airfoils{i} = viterna_extrapolate(af_raw, AR);
    end

    % Run BEM
    opts.beta = bem_beta;
    opts.tsr  = bem_tsr;
    opts.rho  = rho;
    bem  = run_bem(blade, airfoils, ed, opts);
    cp   = bem.cp;
    ct   = bem.ct;
    beta = bem.beta;
    tsr  = bem.tsr;

    fprintf('  NOTE: BEM results have not been validated against a reference.\n');
    fprintf('        Consider running compare_bem_vs_ref() if you have a reference table.\n');
end

fprintf('\n');

%% --- Step 3: Build controller ---
fprintf('[3/5] Building controller...\n');

% % Rated values
P_rated  = dc.VS_RtPwr;
if ~isnan(P_rated_user)
    P_rated = P_rated_user;
end
wg_rated = dc.VS_RefSpd;   % [rad/s] rated generator speed
N        = ed.GBRatio;      % gearbox ratio (1 for direct drive)
R        = ed.TipRad;

% % Optimal TSR and Cp from table
[cp_max, idx]            = max(cp(:));
[i_beta_opt, i_tsr_opt]  = ind2sub(size(cp), idx);
tsr_opt  = tsr(i_tsr_opt);
beta_opt = beta(i_beta_opt);
fprintf('  Cp max = %.4f at beta=%.1f deg, TSR=%.2f\n', cp_max, beta_opt, tsr_opt);

% % Rated torque
M_rated = P_rated / wg_rated;
fprintf('  P rated  = %.3e W\n',    P_rated);
fprintf('  M rated  = %.3e N·m\n',  M_rated);
fprintf('  wg rated = %.4f rad/s\n', wg_rated);

% % Region 2 optimal torque coefficient
%   M_opt = k_M * wg^2   where k_M = 0.5*rho*R^3*pi*R^2*Cp_max/(tsr_opt^3*N^3)
A   = pi * R^2;
k_M = 0.5 * rho * R^3 * A * cp_max / (tsr_opt^3 * N^3);

% % Generator speed transition points
wg_min = dc.VS_MinOMSpd;   % [rad/s] minimum generator speed
wg_15  = wg_min;           % end of region 1
wg_2   = wg_min * 1.1;     % start of region 2 (10% above min)
wg_25  = wg_rated * 0.95;  % start of region 2.5
wg_3   = wg_rated;         % rated speed

% % Torque table
wg_res = 0.001;
wg_r   = 0:wg_res:wg_rated;

idx1  = wg_r < wg_15;                          % Region 1
M_r1  = zeros(1, sum(idx1));

idx15      = wg_r >= wg_15 & wg_r < wg_2;      % Region 1.5
M_15_start = 0;
M_15_end   = k_M * wg_2^2;
M_r15      = M_15_start + (M_15_end - M_15_start) / (wg_2 - wg_15) * (wg_r(idx15) - wg_15);

idx2  = wg_r >= wg_2 & wg_r < wg_25;           % Region 2
M_r2  = k_M * wg_r(idx2).^2;

idx25 = wg_r >= wg_25;                          % Region 2.5
M_25  = k_M * wg_25^2;
M_r25 = M_25 + (M_rated - M_25) / (wg_3 - wg_25) * (wg_r(idx25) - wg_25);
M_r25 = min(M_r25, M_rated);

M_r = [M_r1, M_r15, M_r2, M_r25];

% % Pitch controller gains
% Flat gain scheduling (gs = 1): gains Kp and Ki carry all the information.
% ROSCO table is read to confirm it exists but gains are set by the user.
I_DT = N^2 * ed.GenIner + ed.HubIner + 3*bl.inertia;   % drivetrain inertia at LSS

if ~isempty(dc.PC_GS_angles) && ~isempty(dc.PC_GS_KP)
    fprintf('  ROSCO gain scheduling table found (%d points) - using user-defined Kp/Ki\n', ...
            length(dc.PC_GS_angles));
    beta_gs = (-5:0.2:90)';
    gs      = repmat(ones(1, length(beta_gs)), 6, 1);
    Kp      = pitch_Kp * pitch_gain_scale;
    Ki      = pitch_Ki * pitch_gain_scale;
else
    fprintf('  WARNING: No ROSCO gain scheduling found, computing from scratch.\n');
    w_rot_r = wg_rated / N;
    dpdv    = -28.24e6 * (R/63)^2;   % scaled from NREL 5MW
    Kp      = 2 * I_DT * w_rot_r * ctrl_damp * ctrl_wund / (N * dpdv) * pitch_gain_scale;
    Ki      = I_DT * w_rot_r * ctrl_wund^2 / (N * dpdv) * pitch_gain_scale;
    beta_gs = beta;
    gs      = repmat(ones(1, length(beta)), 6, 1);
end

fprintf('  Kp = %.4f,  Ki = %.4f\n', Kp, Ki);

% % Power levels for gain scheduling
lvls = linspace(P_rated*0.2, P_rated*1.06, 6);

fprintf('\n');

%% --- Step 4: Build FOWLTY structs ---
fprintf('[4/5] Building FOWLTY data structures...\n');

% % env
env.rho = rho;

% Refine Cp/Ct tables to finer grid for smoother interpolation in FOWLTY
beta_fine = linspace(beta(1), beta(end), 200);
tsr_fine  = linspace(tsr(1),  tsr(end),  150);
[TSR_c, BETA_c] = meshgrid(tsr,      beta);
[TSR_f, BETA_f] = meshgrid(tsr_fine, beta_fine);
cp = interp2(TSR_c, BETA_c, cp, TSR_f, BETA_f, 'spline');
ct = interp2(TSR_c, BETA_c, ct, TSR_f, BETA_f, 'spline');
beta = beta_fine;
tsr  = tsr_fine;

% % CP and CT tables
wt.cp.table = cp;
wt.cp.beta  = beta;
wt.cp.tsr   = tsr;
wt.ct.table = ct;
wt.ct.beta  = beta;
wt.ct.tsr   = tsr;

% % Blade properties
wt.blade.length  = ed.TipRad - ed.HubRad;   % [m]
wt.blade.mass    = bl.mass;                  % [kg]  integrated from BMassDen
wt.blade.inertia = bl.inertia;               % [kg·m²] integrated from BMassDen

% % Hub properties
wt.hub.height  = ed.TowerHt + ed.TowerBsHt + ed.Twr2Shft;   % [m]
wt.hub.inertia = ed.HubIner;                                  % [kg·m²]
wt.hub.radius  = ed.HubRad;                                   % [m]
wt.hub.mass    = ed.HubMass;                                   % [kg]

% % Nacelle properties
wt.nac.mass = ed.NacMass + ed.YawBrMass;   % [kg]

% % Tower properties
wt.tower.mass    = tw.mass;        % [kg]  integrated from TMassDen
wt.tower.height  = ed.TowerHt;    % [m]
wt.tower.damp    = tw.damp;        % [-]   from TwrFADmp in ElastoDyn tower file
wt.tower.eigfreq = tower_eigfreq;  % [rad/s] user-defined (see USER INPUT SECTION)

% % Generator properties
wt.gen.inertia      = ed.GenIner;       % [kg·m²] around HSS
wt.gen.N            = N;                % [-] gearbox ratio
wt.gen.effeciency   = sd.GenEff / 100;  % [-]
wt.gen.ratedspeed   = wg_rated;         % [rad/s]
wt.gen.timeconstant = 0.1;              % [s]

% % Initial conditions for Simulink integrators
wt.gen.ini      = wg_rated;   % [rad/s] generator speed at rated
wt.rotor.ini    = wg_rated / N;  % [rad/s] rotor speed at rated
wt.gen.torq.ini = M_rated;    % [N·m]   generator torque at rated

% % Rotor properties (LSS)
wt.rotor.radius     = ed.TipRad;
wt.rotor.inertia    = ed.HubIner + 3*bl.inertia;   % [kg·m²]
wt.rotor.mass       = ed.HubMass + 3*bl.mass;       % [kg]
wt.rotor.ratedspeed = wg_rated / N;                 % [rad/s]
wt.rotor.spring     = ed.DTTorSpr;                  % [N·m/rad]
wt.rotor.damp       = ed.DTTorDmp;                  % [N·m/(rad/s)]

% % Drive train state-space definition
wt.dt.ssA = [-wt.rotor.damp/wt.rotor.inertia          wt.rotor.damp/(wt.rotor.inertia*N)   -wt.rotor.spring/wt.rotor.inertia; ...
              wt.rotor.damp/(wt.gen.inertia*N)         -wt.rotor.damp/(wt.gen.inertia*N^2)   wt.rotor.spring/(wt.gen.inertia*N); ...
              1                                        -1/N                                   0];
wt.dt.ssB = [1/wt.rotor.inertia   0; ...
             0                    -1/wt.gen.inertia; ...
             0                     0];
wt.dt.ssC = blkdiag(eye(2), wt.rotor.spring);
wt.dt.ssD = zeros(3,2);

% % Second order pitch actuator properties
if ~isnan(dc.PA_CornerFreq)
    pitch_omega = dc.PA_CornerFreq;
    pitch_damp  = dc.PA_Damping;
else
    pitch_omega = 11.11;   % NREL 5MW default [rad/s]
    pitch_damp  = 0.6;
end
wt.pitch.omega    = pitch_omega;
wt.pitch.damp     = pitch_damp;
wt.pitch.ssA      = [0 1; -pitch_omega^2 -2*pitch_omega*pitch_damp];
wt.pitch.ssB      = [0; pitch_omega^2];
wt.pitch.ssC      = [1 0];
wt.pitch.ssD      = 0;
wt.pitch.sysPitch = ss([0 1; -pitch_omega^2 -2*pitch_omega*pitch_damp], ...
                        [0; pitch_omega^2], [1 0], 0);

% % Tower top properties
wt.top.mass = wt.nac.mass + wt.rotor.mass;

% % Torque controller
wt.ctrl.torq.wg      = wg_r;
wt.ctrl.torq.M       = M_r;
wt.ctrl.torq.r3      = wg_3;
wt.ctrl.torq.lim     = M_rated;
wt.ctrl.torq.ratelim = dc.VS_MaxRat;
wt.ctrl.torq.beta    = 0;

% % Pitch controller
wt.ctrl.pitch.beta    = beta_gs;
wt.ctrl.pitch.pwr     = lvls;
wt.ctrl.pitch.gs      = gs;
wt.ctrl.pitch.Pgain   = Kp;
wt.ctrl.pitch.Igain   = Ki;
wt.ctrl.pitch.ratelim = abs(dc.PC_MaxRat) * 180/pi;   % [deg/s]
wt.ctrl.pitch.ulim    = dc.PC_MaxPit * 180/pi;         % [deg]
wt.ctrl.pitch.llim    = dc.PC_MinPit * 180/pi;         % [deg]

% % Other control parameters
timestep = 0.0125;
fcorner  = 0.25;
alpha    = exp(-2*pi*timestep*fcorner);

wt.ctrl.gen.rated      = wg_rated;
wt.ctrl.p_rated        = P_rated;
wt.ctrl.alpha          = alpha;
wt.ctrl.Ts             = timestep;
wt.ctrl.gen.effeciency = wt.gen.effeciency;

% % Measurement system values
wt.meas.Ts          = timestep;
wt.meas.pitchres    = 0.01;
wt.meas.rotres      = 0.001;
wt.meas.aneres      = 0.001;
wt.meas.genres      = 0.001;
wt.meas.nacres      = 0.01;
wt.meas.delay       = 0.01;
wt.meas.rottau      = 0.01;
wt.meas.gentau      = 0.01;
wt.meas.anenoise    = 0.0071;
wt.meas.rotSpdnoise = 1e-4;
wt.meas.genSpdnoise = 2e-4;
wt.meas.genTrqnoise = 0.9;
wt.meas.genPwrnoise = 10;
wt.meas.pitchnoise  = 1.5e-3;
wt.meas.nacAccnoise = 5e-4;

% % Save public values
pub.rotor.radius  = ed.TipRad;
pub.rated         = P_rated;
pub.rotor.Nslope  = 10;
pub.rotor.Mult    = 10e10;
pub.tower.Nslope  = 4;
pub.tower.Mult    = 10e10;
pub.shaft.Nslope  = 4;
pub.shaft.Mult    = 10e10;

fprintf('\n');

%% --- Step 4b: Platform ---
fprintf('[4b] Building platform struct...\n');

switch lower(platform_mode)
    case 'none'
        fprintf('  Platform mode: none (fixed-bottom, zero motion)\n');
        nDoF_p           = 2;    % surge + pitch
        plat.DoF         = [1 3];
        plat.nDoF        = nDoF_p;
        plat.CGz         = 0;
        plat.d_N2P       = ed.TowerHt + ed.Twr2Shft;
        plat.nF2pT       = [zeros(nDoF_p-1,1); plat.d_N2P];
        plat.pM2nM_surge = [0 1 zeros(1,2*(nDoF_p-1))];
        plat.pM2nM_pitch = [zeros(1,2*(nDoF_p-1)) 0 plat.d_N2P];
        n_states         = 2*nDoF_p;
        plat.HydSys      = ss(zeros(n_states), zeros(n_states,nDoF_p), ...
                              zeros(2*nDoF_p,n_states), zeros(2*nDoF_p,nDoF_p));
        fprintf('  Zero-motion platform SS built (order %d)\n', n_states);

    case 'scaled'
        deepcwind_path = 'DeepCWind_Original.mat';
        platform_DoF   = [1 3];
        if ~exist(deepcwind_path, 'file')
            error('Platform mode is ''scaled'' but DeepCWind file not found:\n  %s\nPlace DeepCWind_Original.mat in the working directory.', deepcwind_path);
        end
        plat = scale_platform(deepcwind_path, ed.TipRad, platform_DoF, ed);

    case 'custom'
        if isempty(fieldnames(plat))
            error('Platform mode is ''custom'' but plat struct is empty.\nFill in the plat struct in the USER INPUT SECTION.');
        end
        fprintf('  Platform mode: custom (user-defined)\n');

    otherwise
        error('Unknown platform_mode: ''%s''. Use ''none'', ''scaled'', or ''custom''.', platform_mode);
end

fprintf('\n');

%% --- Step 5: Save ---
fprintf('[5/5] Saving...\n\n');

fprintf('  Summary:\n');
fprintf('    Blade mass    : %.1f kg\n',                wt.blade.mass);
fprintf('    Blade inertia : %.4e kg·m²\n',            wt.blade.inertia);
fprintf('    Tower mass    : %.1f kg\n',                wt.tower.mass);
fprintf('    Tower eigfreq : %.4f rad/s (%.4f Hz)\n',  wt.tower.eigfreq, wt.tower.eigfreq/(2*pi));
fprintf('    Rotor inertia : %.4e kg·m²\n',            wt.rotor.inertia);
fprintf('    Cp max        : %.4f\n',                    max(wt.cp.table(:)));
fprintf('    Pitch Kp      : %.4f\n',                    wt.ctrl.pitch.Pgain);
fprintf('    Pitch Ki      : %.4f\n',                    wt.ctrl.pitch.Igain);
fprintf('    Platform mode : %s\n',                      platform_mode);

save_path = [output_name '.mat'];
if save_output
    save(save_path, 'wt', 'env', 'pub', 'plat');
    fprintf('  Saved: %s\n', save_path);
else
    fprintf('  Save skipped (save_output = false)\n');
end
fprintf('\n=============================================================\n');
fprintf(' Done! Load with: load(''%s'')\n', save_path);
fprintf('=============================================================\n');

%% --- Optional: Load data into Simulink library block ---
if interactive
    if ~save_output
        warning('interactive=true pero save_output=false: el bloque se cargará en Simulink pero el .mat no se habrá guardado en disco.');
    end
    p.wt     = wt;
    p.env    = env;
    p.public = pub;
    p.plat   = plat;
    disp(' ');
    disp('Select and highlight the turbine block in the Simulink library');
    disp('and press Enter to save the turbine data into the block.');
    disp('Remember to unlock the library first (right-click -> Library -> Unlock).');
    pause
    set_param(gcb, 'w_gen_ini', num2str(wg_rated));
    set_param(gcb, 'auto', 'on');
    set_param(gcb, 'UserData', p);
    set_param(gcb, 'UserDataPersistent', 'on');
    fprintf('Turbine data saved into Simulink block.\n');
end