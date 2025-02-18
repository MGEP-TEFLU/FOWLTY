%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation for a faulty offshore wind farm simulation
%
%
% Author: Yerai Peña-Sanchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
warning('off')

%% Definition of the simulation case
SimulinkModel               = 'farmName';
WindCase                    = 'windName'; 

[turb,wind]                 = setSimulation(SimulinkModel,WindCase);
nT                          = length(turb.farm.x);

%% Definition of simulation variables
tMax                        = wind.SimTime;         % [s] maximum simulation time
dt                          = wind.Ts;              % [s] Sampling time for the definition of faults

noiseFlag                   = 0;                    % [-] Flag to activate or not measurement noise

tsim                        = (dt:dt:tMax).';       % [s] Vector of the simulation time
    
%% Wave generation
Platform                    = load('DeepCWind/SS/DeepCWind.mat');

Hs                          = 2;                    % Significant wave height [m]
Tw                          = 10;                   % Tipical wave period [s]

[eta,Fex,~,~,~]             = whitenoiseWave(dt,tMax,Hs,Tw,Platform.w,Platform.Fe.',turb,nT,1);

%% Definition of faults
scenario                    = 0;

[faultSignal,faultOptions]  = faultScenarios(scenario,tsim,nT); 

%% Run simulink model
tic
sim(SimulinkModel);
toc

%% Results comparisson

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
plotFlag = 1;
if plotFlag == 1
    newFig = 0;     % 0 - Plot on top of the same figure
                    % 1 - Plot using new figure
    
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









