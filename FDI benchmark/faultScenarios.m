function [faultSignal,faultOptions] = faultScenarios(scenario,t,nT)
%% Definition and selection of pre-defined fault scenarios (simplified benchmar version) ------------------------------------------
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
% - Generator power sensor fault:
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
faultOptions.PitchActOffset         = 3;      % [deg]      Offset on the actuator Pitch
faultOptions.PitchActGain           = 1.1;     % [-]        Gain on the actuator Pitch
faultOptions.PitchActDriftSlope     = 1e-5;    % [deg/s]    Slope of the measurement drifting 
faultOptions.pitchNatFreq           = 5.73;    % [rad/s]    New natural frequency of the pitch actuator (non faulty = 11.11)
faultOptions.pitchDamp              = 0.45;    % [-]        New damping factor of the pitch actuator (non faulty = 0.6)

% Pitch sensor fault --------------------------------------------------------------------------------------------------------------
faultOptions.PitchSenOffset         = -1;      % [deg]      Offset on the Pitch sensor
faultOptions.PitchSenGain           = 1.1;     % [-]        Gain on the pitch sensor
faultOptions.PitchSenDriftSlope     = 4e-5;    % [deg/s]    Slope of the measurement drifting 
faultOptions.PitchSenPrecision      = 0.1;     % [0-1]      Precision error of the sensor (0 for normal operation)

% Rotor speed sensor fault ----------------------------------------------------------------------------------------------------
faultOptions.RotSenOffset           = .1;      % [rad/s]    Offset on the rotor speed sensor
faultOptions.RotSenGain             = 1.1;     % [-]        Gain on the rotor speed sensor
faultOptions.RotSenDriftSlope       = 1e-5;    % [rad/s^2]  Slope of the measurement drifting 
faultOptions.RotSenPrecision        = 0.05;    % [0-1]      Precision error of the sensor (0 for normal operation)

% Generator speed sensor fault ----------------------------------------------------------------------------------------------------
faultOptions.GenSenOffset           = 2;       % [rad/s]    Offset on the generator speed sensor
faultOptions.GenSenGain             = 1.05;    % [-]        Gain on the generator speed sensor
faultOptions.GenSenDriftSlope       = 1e-3;    % [rad/s^2]  Slope of the measurement drifting 
faultOptions.GenSenPrecision        = 0.1;     % [0-1]      Precision error of the sensor (0 for normal operation)

% Generator power sensor fault ----------------------------------------------------------------------------------------------------
faultOptions.GenPSenOffset          = 200;     % [W]        Offset on the generator power sensor
faultOptions.GenPSenGain            = 1.1;     % [-]        Gain on the generator power sensor
faultOptions.GenPSenDriftSlope      = 1;       % [W/s]      Slope of the measurement drifting 
faultOptions.GenPSenPrecision       = 0.1;     % [0-1]      Precision error of the sensor (0 for normal operation)

% Generator fault -----------------------------------------------------------------------------------------------------------------
faultOptions.GenActOffset           = 1e3;     % [Nm]       Offset on the generator torque
faultOptions.GenActGain             = 1.1;     % [-]        Gain on the generator torque
faultOptions.GenActDriftSlope       = 1e-3;    % [rad/s^2]  Slope of the generator torque drifting 
faultOptions.GenActEffectiveness    = 1000;    % [0-1]      Effectiveness of the generator (1 for normal operation)

% Predefinition of the fault signals to 0 (no fault)
fault_PitchActuator                 = zeros(length(t),nT);
fault_PitchSensor                   = zeros(length(t),nT);
fault_RotSensor                     = zeros(length(t),nT);
fault_GenSensor                     = zeros(length(t),nT);
fault_GenPSensor                    = zeros(length(t),nT);
fault_GenActuator                   = zeros(length(t),nT);

%% Benchmark case fault scenarios
load data/randomFaultVariables.mat           % Load saved random variables

sc = scenario;                          % Just to make the code shorter
if sc <= 0 || sc > 10                   % No fault
else
    faultOptions.PitchActOffset         = randPitchActOffset(sc);
    faultOptions.pitchNatFreq           = randpitchNatFreq(sc);
    faultOptions.pitchDamp              = randpitchDamp(sc);
    faultOptions.PitchSenGain           = randPitchSenGain(sc);
    faultOptions.PitchSenDriftSlope     = randPitchSenDriftSlope(sc);
    faultOptions.GenSenGain             = randGenSenGain(sc);
    faultOptions.GenPSenGain            = randGenPSenGain(sc);
    faultOptions.GenActOffset           = randGenActOffset(sc);

    fault_PitchSensor((  t > fTimes(sc,1)) & (t < fTimes(sc,1)+fDeltaT(sc,1)),rTurb(sc,1)) = 1;
    fault_PitchSensor((  t > fTimes(sc,2)) & (t < fTimes(sc,2)+fDeltaT(sc,2)),rTurb(sc,1)) = 3;
    fault_PitchSensor((  t > fTimes(sc,3)) & (t < fTimes(sc,3)+fDeltaT(sc,3)),rTurb(sc,1)) = 4;
    fault_GenPSensor((   t > fTimes(sc,4)) & (t < fTimes(sc,4)+fDeltaT(sc,4)),rTurb(sc,1)) = 4;
    fault_GenSensor((    t > fTimes(sc,5)) & (t < fTimes(sc,5)+fDeltaT(sc,5)),rTurb(sc,1)) = 4;
    fault_PitchActuator((t > fTimes(sc,6)) & (t < fTimes(sc,6)+fDeltaT(sc,6)),rTurb(sc,1)) = 1;
    fault_PitchActuator((t > fTimes(sc,7)) & (t < fTimes(sc,7)+fDeltaT(sc,7)),rTurb(sc,1)) = 2;
    fault_PitchActuator((t > fTimes(sc,8)) & (t < fTimes(sc,8)+fDeltaT(sc,8)),rTurb(sc,1)) = 4;
    fault_GenActuator((  t > fTimes(sc,9)) & (t < fTimes(sc,9)+fDeltaT(sc,9)),rTurb(sc,1)) = 2;
end

%% Put all the fault signals together ----------------------------------------------------------------------------------------------
faultSignal.time                    = t;
faultSignal.signals.values          = kron(fault_PitchActuator,[1 0 0 0 0 0])+...     % Pitch actuator fault
                                      kron(fault_PitchSensor,[0 1 0 0 0 0])+...       % Pitch sensor fault
                                      kron(fault_RotSensor,[0 0 1 0 0 0])+...         % Rotor speed sensor fault
                                      kron(fault_GenSensor,[0 0 0 1 0 0])+...         % Generator speed sensor fault
                                      kron(fault_GenPSensor,[0 0 0 0 1 0])+...        % Generator power sensor fault
                                      kron(fault_GenActuator,[0 0 0 0 0 1]);          % Drive train fault
end

