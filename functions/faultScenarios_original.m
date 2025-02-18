function [faultSignal,faultOptions] = faultScenarios(scenario,t,nT)
%% Definition and selection of pre-defined fault scenarios ------------------------------------------------------------------------
% Inputs:
%   - scenario          number of a pre-defined scenario
%   - t                 time vector of the simulation
%   - nT                number of turbines composing the farm
% Outputs:
%   - fault_signal      signal containing the fault information for simulink
%
% Created by: 
%   Yerai Peña-Sanchez (2022)
% ---------------------------------------------------------------------------------------------------------------------------------

%% Definition of possible faults
% - Pitch actuator faults:
%       1 - stuck at a fixed value
%       2 - constant offset fault
%       3 - drifting fault
%       4 - constant gain fault
%       5 - change of dynamics
% - Pitch sensor faults:
%       1 - stuck at a fixed value
%       2 - constant offset fault
%       3 - drifting fault
%       4 - constant gain fault
%       5 - precision degradation
% - Generator speed sensor fault:
%       1 - stuck at a fixed value
%       2 - constant offset fault
%       3 - drifting fault
%       4 - constant gain fault
%       5 - precision degradation
% - Generator faults:
%       1 - stuck at a fixed value
%       2 - constant offset fault
%       3 - drifting fault
%       4 - constant gain fault
%       5 - loss of effectiveness
%
% 0 is no fault situation for all the cases

%% Definition of fault parameters
% Pitch actuator fault ------------------------------------------------------------------------------------------------------------
faultOptions.PitchActOffset         = -1;      % [deg]      Offset on the actuator Pitch
faultOptions.PitchActGain           = 1.2;     % [-]        Gain on the actuator Pitch
faultOptions.PitchActDriftSlope     = 1e-5;    % [deg/s]    Slope of the measurement drifting 
faultOptions.pitchNatFreq           = 5.73;    % [rad/s]    New natural frequency of the pitch actuator (non faulty = 11.11)
faultOptions.pitchDamp              = 0.45;    % [-]        New damping factor of the pitch actuator (non faulty = 0.6)

% Pitch sensor fault --------------------------------------------------------------------------------------------------------------
faultOptions.PitchSenOffset         = -1;      % [deg]      Offset on the Pitch sensor
faultOptions.PitchSenGain           = 1.2;     % [-]        Gain on the pitch sensor
faultOptions.PitchSenDriftSlope     = 4e-5;    % [deg/s]    Slope of the measurement drifting 
faultOptions.PitchSenPrecision      = 0.1;     % [0-1]      Precision error of the sensor (0 for normal operation)

% Rotor speed sensor fault ----------------------------------------------------------------------------------------------------
faultOptions.RotSenOffset           = .1;      % [rad/s]    Offset on the rotor speed sensor
faultOptions.RotSenGain             = 1.2;     % [-]        Gain on the rotor speed sensor
faultOptions.RotSenDriftSlope       = 1e-5;    % [rad/s^2]  Slope of the measurement drifting 
faultOptions.RotSenPrecision        = 0.05;    % [0-1]      Precision error of the sensor (0 for normal operation)

% Generator speed sensor fault ----------------------------------------------------------------------------------------------------
faultOptions.GenSenOffset           = 1;       % [rad/s]    Offset on the generator speed sensor
faultOptions.GenSenGain             = 1.2;     % [-]        Gain on the generator speed sensor
faultOptions.GenSenDriftSlope       = 1e-3;    % [rad/s^2]  Slope of the measurement drifting 
faultOptions.GenSenPrecision        = 0.1;    % [0-1]       Precision error of the sensor (0 for normal operation)

% Generator fault -----------------------------------------------------------------------------------------------------------------
faultOptions.GenActOffset           = 1e3;     % [Nm]       Offset on the generator torque
faultOptions.GenActGain             = 1.2;     % [-]        Gain on the generator torque
faultOptions.GenActDriftSlope       = 1e-3;    % [rad/s^2]  Slope of the generator torque drifting 
faultOptions.GenActEffectiveness    = 1000;    % [0-1]      Effectiveness of the generator (1 for normal operation)

% Predefinition of the fault signals to 0 (no fault)
fault_PitchActuator                 = zeros(length(t),nT);
fault_PitchSensor                   = zeros(length(t),nT);
fault_RotSensor                     = zeros(length(t),nT);
fault_GenSensor                     = zeros(length(t),nT);
fault_GenActuator                   = zeros(length(t),nT);

%% Pre-defined fault scenarios
% To define a new case change the 0 values (no fault) of the matrices of the different possible faults for the desired fault number 
% at the desired time instants (using columns) and for the desired turbines (using rows). 
% By way of example, fault_PitchActuator((t > 5) & (t < 50),2) = 1 sets a pitch actuator fault of type 1 (stuck at a fixed value) from the 
% second 5 to the second 50 of the simularion on the turbine number 2 of the farm.

switch scenario
    case 0 % No fault

end


%% Put all the fault signals together ----------------------------------------------------------------------------------------------
faultSignal.time                    = t;
faultSignal.signals.values          = kron(fault_PitchActuator,[1 0 0 0 0])+...     % Pitch actuator fault
                                      kron(fault_PitchSensor,[0 1 0 0 0])+...       % Pitch sensor fault
                                      kron(fault_RotSensor,[0 0 1 0 0])+...         % Rotor speed sensor fault
                                      kron(fault_GenSensor,[0 0 0 1 0])+...         % Generator speed sensor fault
                                      kron(fault_GenActuator,[0 0 0 0 1]);          % Drive train fault
end

