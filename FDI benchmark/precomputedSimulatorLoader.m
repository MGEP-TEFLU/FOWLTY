%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation for a faulty offshore wind farm simulation
%
%
% Author: Yerai Peña-Sanchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
warning('off')

windSpeeds                      = (5:2:23);
nT                              = 7;                    % Number of turbines composing the farm
faultFlag                       = 1;                    % 0 or 1 to select between a healthy (0) or a faulty (1) simulation
noiseFlag                       = 1;                    % 0 or 1 to activate (1) or not (0) measurement noise

for dataSet = 1:10
    
    if faultFlag == 1
        if noiseFlag == 1
            name                = ['precomputed simulations/faulty_noisy_' int2str(windSpeeds(dataSet)) 'ms'];
        else
            name                = ['precomputed simulations/faulty_noNoise_' int2str(windSpeeds(dataSet)) 'ms'];
        end
    else
        if noiseFlag == 1
            name                = ['precomputed simulations/noFault_noisy_' int2str(windSpeeds(dataSet)) 'ms'];
        else
            name                = ['precomputed simulations/noFault_noNoise_' int2str(windSpeeds(dataSet)) 'ms'];
        end
    end

    load(name)

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
            plot(t,betaRef(:,((i-1)*3+1:3*i)),t,betaMeas(:,((i-1)*3+1:3*i)),'--','linewidth',2)
            title('Blade pitch [deg]'); grid on; hold on
    
            subplot(vPlot,hPlot,2)
            plot(t,w_rot(:,i),'linewidth',2)
            title('Rotor speed [rpm]'); grid on; hold on
            
            subplot(vPlot,hPlot,3)
            plot(t,power_gen(:,i),'linewidth',2)
            title('Generator power [kW]'); grid on; hold on
            
            subplot(vPlot,hPlot,hPlot+1)
            plot(t,torque_gen(:,i),t,torque_genRef(:,i),'--','linewidth',2)
            title('Generator torque [kNm]'); grid on; hold on
            
            subplot(vPlot,hPlot,hPlot+2)
            plot(t,w_gen(:,i),'linewidth',2)
            title('Generator speed [rpm]'); grid on; hold on
            
            subplot(vPlot,hPlot,hPlot+3)
            plot(t,v_wind(:,i),'linewidth',2)
            title('Wind speed [m/s]'); grid on; hold on
        end
    end
end








