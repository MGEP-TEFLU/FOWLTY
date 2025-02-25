%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create different wind field realisations from a given wind field
%
%
% Author: Yerai Peña-Sanchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars

%% Load the case of wind field you want to replicate
originalName        = 'data/wind_11ms_1000s';
original            = load(originalName);

windData.Ua         = original.wind.input.Ua;          % [m/s]  Mean wind speed
windData.Ti         = original.wind.input.Ti;          %        Turbulence intensity
windData.d          = original.wind.input.d;           % [m]    Grid size
windData.Lx         = original.wind.input.Lx;          % [m]    Wind field length
windData.Ly         = original.wind.input.Ly;          % [m]    Wind field width
windData.SimTime    = original.wind.input.SimTime;     % [s]    Simulation length
windData.desc       = original.wind.input.desc;        %        Description
windData.Ts         = original.wind.input.Ts;          % [s]    Sampling time

%% Re-define (if needed) a variable
% Note that if Lx and Ly are changed, the simulation will not run.
windData.Ua         = 12;                            % [m/s]  New mean wind speed
% windData.Ts         = 0.01;                         % [s]    New sampling time
% windData.Ti         = 0.1;                          %        New Turbulence intensity
% windData.d          = 10;                           % [m]    New grid size
% windData.SimTime    = 35;                            % [s]    New simulation length

%% Creating a new wind field realisation
% A couple of options
opt.parallel        = 'Yes';                        % [Yes or No] Use parallel computation or not
opt.saveName        = ['wind_' int2str(windData.Ua) 'ms_' int2str(windData.SimTime) 's.mat'];                % [str] Name to save the new wind field file (comment to avoid saving)

wind                = recreate_windfield(windData,opt);


