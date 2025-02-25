%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation for a faulty offshore wind farm simulation
%
%
% Author: Yerai Peña-Sanchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
warning('off')

%% Definition of the simulation case
dataSet                     = 1;                    % From 1 to 10, each one representing a different wind 
                                                    % speed from 5m/s to 23m/s and a different fault scenario

faultFlag                   = 0;                    % 0 or 1 to select between a healthy (0) or a faulty (1) simulation
noiseFlag                   = 0;                    % 0 or 1 to activate (1) or not (0) measurement noise


SimulinkModel               = 'benchmarkWindfarm';  % Name of the simulink wind farm file

windSpeeds                  = (5:2:23);
WindCase                    = ['data/wind_' int2str(windSpeeds(dataSet)) 'ms_1000s'];

[turb,wind]                 = setSimulation(SimulinkModel,WindCase);
nT                          = length(turb.farm.x);  % Number of turbines composing the farm

%% Definition of simulation variables
tMax                        = wind.SimTime;         % [s] maximum simulation time
dt                          = wind.Ts;              % [s] Sampling time for the definition of faults


tsim                        = (dt:dt:tMax).';       % [s] Vector of the simulation time
    
%% Wave generation
Platform                    = load('DeepCWind/SS/DeepCWind.mat');

Hs                          = 2;                    % [m] Significant wave height
Tw                          = 10;                   % [s] Tipical wave period

[eta,Fex,~,~,~]             = whitenoiseWave(dt,tMax,Hs,Tw,Platform.w,Platform.Fe.',turb,nT,1);

%% Definition of faults
if faultFlag == 1
    scenario                = dataSet;
else
    scenario                = 0;
end

[faultSignal,faultOptions]  = faultScenarios(scenario,tsim,nT); 

%% Run simulink model
tic
sim(SimulinkModel);
toc

%% Results comparisson
% All the logged signals from simulink can be obtained via the logsout variable in matlab. Here, some of them 
% are saved and ploted, however, there are more variables from the simulation can be found in logsout.

t                           = logsout.getElement('pitch').Values.Time;
beta_all                    = logsout.getElement('pitch').Values.Data;                  % [deg]
nBlades                     = (size(beta_all,2)/nT)/3;                                  % To check if the definition of the blades is made independently 
                                                                                        % (faulty case) or the three blades together (healthy case)
C_betaRef                   = kron(eye(nT),[zeros(nBlades) eye(nBlades) zeros(nBlades)]);
C_betaMeas                  = kron(eye(nT),[zeros(nBlades) zeros(nBlades) eye(nBlades)]);
C_beta                      = kron(eye(nT),[eye(nBlades) zeros(nBlades) zeros(nBlades)]);

beta                        = (C_beta*beta_all.').';                                    % [deg] Real pitch angle
betaRef                     = (C_betaRef*beta_all.').';                                 % [deg] Reference pitch angle
betaMeas                    = (C_betaMeas*beta_all.').';                                % [deg] Measured pitch angle

power_gen                   = logsout.getElement('P_farm').Values.Data/1000;            % [kW]  (originally [W])
torque_gen                  = logsout.getElement('M_gen').Values.Data/1000;             % [kNm] (originally [Nm]) 
torque_rot                  = logsout.getElement('M_shaft').Values.Data/1000;           % [kNm] (originally [Nm])
w_gen                       = logsout.getElement('w_gen').Values.Data*30/pi;            % [rpm] (originally [rad/s])

v_nac                       = logsout.getElement('V_meas').Values.Data;                 % [m/s]

%% Plots
% Some of the variables from the Simulink simulation are plotted here in a different plot for each turbine. 
% The selection of the variables to plot is just an example and a different set can be considered.

plotFlag = 1;       % 0 or 1 to plot the results (1) or not (0)
if plotFlag == 1
    newFig = 0;     % 0 or 1 to plot the results on top of the same figure (0) or in a new figure (1)
    
    for i = 1:nT
        if newFig == 0      
            f = figure(100+i);
        else
            f = figure;
        end
        f.WindowState = 'maximized';
        
        vPlot = 2;
        hPlot = 3;
        
        subplot(vPlot,hPlot,1)
        plot(t,beta(:,((i-1)*nBlades+1:nBlades*i)),'linewidth',2)
        title('Blade pitch [deg]'); grid on; hold on
        
        subplot(vPlot,hPlot,2)
        plot(t,torque_rot(:,i),'linewidth',2)
        title('Rotor torque [kNm]'); grid on; hold on
        
        subplot(vPlot,hPlot,3)
        plot(t,power_gen(:,i),'linewidth',2)
        title('Generator power [kW]'); grid on; hold on
        
        subplot(vPlot,hPlot,hPlot+1)
        plot(t,torque_gen(:,i),'linewidth',2)
        title('Generator torque [kNm]'); grid on; hold on
        
        subplot(vPlot,hPlot,hPlot+2)
        plot(t,w_gen(:,i),'linewidth',2)
        title('Generator speed [rpm]'); grid on; hold on
        
        subplot(vPlot,hPlot,hPlot+3)
        plot(t,v_nac(:,i),'linewidth',2)
        title('Wind speed [m/s]'); grid on; hold on
    end
end







