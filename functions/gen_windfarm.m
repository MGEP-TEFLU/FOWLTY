function mainsys = gen_windfarm(mdlfile,path,farm,wind)
% MDLCREATE - creates wind farm model based on programmatic input
%
% Copyright 2009 - Aalborg University
% Author Jacob Grunnet
% Modified Martin Kragelund
% Modified Yerai Peña-Sanchez
% 

% % Extract grid information from wind struct
grid       = wind.grid;
% % Extract turbine grid positions
% % % % % % % % % % % % % % % % % % wind.xgrids = floor(farm.pos(1,:)/grid.size)+1;
posTot     = [farm.pos farm.posMast];
wind.xgrids = floor(posTot(1,:)/grid.size)+1;

% % Remove mdlfile extension and store in mdlname
mdlname    = mdlfile(1:end-4);

% % Find the number of turbines
Nstr       = int2str(length(farm.turbines));

% % Load the aeolus library
load_system('lib_FOWLTY');
% % Close the model if it is already open
bdclose(mdlname);
% % Create the new model based on the template
mainsys = new_system(mdlname,'model','lib_FOWLTY/Farm Template - Taylor');
% % Save and break all links in order to modify the new library
save_system(mainsys,[path mdlfile],'BreakUserLinks',true)

% % Setup muxes according to number of turbines
set_param([mdlname '/Turbines/Pdemanddemux'],'Outputs',Nstr);
set_param([mdlname '/Turbines/Vrotdemux'],'Outputs',Nstr);
set_param([mdlname '/Turbines/Vnacdemux'],'Outputs',Nstr);
set_param([mdlname '/Turbines/Fexdemux'],'Outputs',Nstr);
set_param([mdlname '/Turbines/Faultdemux'],'Outputs',Nstr);

set_param([mdlname '/Turbines/Pmux'],'Inputs',Nstr);
set_param([mdlname '/Turbines/CTmux'],'Inputs',Nstr);
set_param([mdlname '/Turbines/wgenmux'],'Inputs',Nstr);
set_param([mdlname '/Turbines/vwindmux'],'Inputs',Nstr);
set_param([mdlname '/Turbines/pitchmux'],'Inputs',Nstr);
set_param([mdlname '/Turbines/mshaftmux'],'Inputs',Nstr);
set_param([mdlname '/Turbines/mtowmux'],'Inputs',Nstr);
set_param([mdlname '/Turbines/mgenmux'],'Inputs',Nstr);
set_param([mdlname '/Turbines/anacmux'],'Inputs',Nstr);
set_param([mdlname '/Turbines/auxmux'],'Inputs',Nstr);

% % Set toolbox vesion
% % Fake version atm
create_windfarm_paths = which('create_windfarm','-all');

if(isunix)
    create_windfarm_paths = create_windfarm_paths{end}(2:(end-18));
elseif(ispc)
    create_windfarm_paths = create_windfarm_paths{end}(4:(end-18));
end
try
    v          = ver(create_windfarm_paths);
    p.version  = v.Version;
catch me
    disp('Error in getting FOWLTY version - setting default version 0.0.0')
    disp('You might have just opened the library instead of installing it. Use install.m')
    p.version = '0.0.0';
end
set_param([mdlname '/Turbines'],'UserData',p);
set_param([mdlname '/Turbines'],'UserDataPersistent','on');



%% Add turbines
% % Variables for block positioning
y = 0;
x = 150;
% % Initial x and y
tl1x = 50;
tl1y = 50;
ptotal = 0;
% % Load turbine library
load_system('libturbines');


for i = 1:length(farm.turbines)
% %     Add a turbine block
    wt(i).block = add_block(['libturbines/' farm.turbines{i}],[mdlname '/Turbines/WT' int2str(i)]);
    
% %     Extract turbine data
    partmp = get_param(wt(i).block,'UserData');
% %     Save the public data
    wtpar(i) = partmp.public;
    
% %     Set block position
    wtp = get_param(wt(i).block,'Position');
    set_param(wt(i).block,'Position',[tl1x+x tl1y+y,wtp(3)+x-wtp(1) wtp(4)+y-wtp(2)]);
       
% %     Increase y position
    y = y+175;

% %     Input - connect muxes
    add_line([mdlname '/Turbines'],['Vrotdemux/' int2str(i)],['WT' int2str(i) '/2'],'Autorouting','on')
    add_line([mdlname '/Turbines'],['Vnacdemux/' int2str(i)],['WT' int2str(i) '/3'],'Autorouting','on')
    add_line([mdlname '/Turbines'],['Pdemanddemux/' int2str(i)],['WT' int2str(i) '/1'],'Autorouting','on')
    add_line([mdlname '/Turbines'],['Fexdemux/' int2str(i)],['WT' int2str(i) '/4'],'Autorouting','on')
    add_line([mdlname '/Turbines'],['Faultdemux/' int2str(i)],['WT' int2str(i) '/5'],'Autorouting','on')
    
    
% %     Output - connect muxes
    add_line([mdlname '/Turbines'],['WT' int2str(i) '/1'],['Pmux/' int2str(i)],'Autorouting','on')
    add_line([mdlname '/Turbines'],['WT' int2str(i) '/2'],['CTmux/' int2str(i)],'Autorouting','on')
    add_line([mdlname '/Turbines'],['WT' int2str(i) '/3'],['wgenmux/' int2str(i)],'Autorouting','on')
    add_line([mdlname '/Turbines'],['WT' int2str(i) '/4'],['vwindmux/' int2str(i)],'Autorouting','on')
    add_line([mdlname '/Turbines'],['WT' int2str(i) '/5'],['pitchmux/' int2str(i)],'Autorouting','on')
    add_line([mdlname '/Turbines'],['WT' int2str(i) '/6'],['mshaftmux/' int2str(i)],'Autorouting','on')
    add_line([mdlname '/Turbines'],['WT' int2str(i) '/7'],['mtowmux/' int2str(i)],'Autorouting','on')
    add_line([mdlname '/Turbines'],['WT' int2str(i) '/8'],['mgenmux/' int2str(i)],'Autorouting','on')
    add_line([mdlname '/Turbines'],['WT' int2str(i) '/9'],['anacmux/' int2str(i)],'Autorouting','on')
    add_line([mdlname '/Turbines'],['WT' int2str(i) '/10'],['auxmux/' int2str(i)],'Autorouting','on')
    
    
% %     Farm controller parameters
    pfarm.radius(i) = wtpar(i).rotor.radius;
    pfarm.rated(i) = wtpar(i).rated;
    pfarm.Cp(i) = 0.45; %Not currently available in turbines
    ptotal = ptotal+wtpar(i).rated;
end

% % %Setup farm control block
pfarm.rho = 1.2231;
pfarm.N = length(farm.turbines);

set_param([mdlname '/Farm Control'],'UserData',pfarm);
set_param([mdlname '/Farm Control'],'UserDataPersistent','on');

% % Setup network operator block
set_param([mdlname '/Network'],'pref',num2str(ptotal));
set_param([mdlname '/Network'],'maxp',num2str(ptotal));
set_param([mdlname '/Network'],'minp',num2str(ptotal/5));

% % Setup network operator block
set_param([mdlname '/Network/Grid'],'fgain',num2str(ptotal*1e-16));


% % Setup Power load block
set_param([mdlname '/Network/Network Load'],'maxp',num2str(ptotal));

%% Prepare wake and wind blocks
% % Calculate common parameters
farmx = farm.pos(1,:);

% % Calculate static wake expansion
for i = 1:length(farm.turbines)
% %     Meandering radius at all grid spacings down wind
    meand(i,:) = sqrt((wtpar(i).rotor.radius)^2+(0:grid.size:grid.xsize).*wtpar(i).rotor.radius/2);
end

% % % % Extract wake expansion at turbines
% % for i = 1:length(farm.turbines)
% %     for j = 1:length(farm.turbines)         % Length from j to i
% %         gdelta = max(min(farm.pos(1,i)-farm.pos(1,j),grid.xsize),0);
% %         
% %         if (gdelta > 0)                     % For downwind turbines
% %             gdeltag        = floor(gdelta/grid.size);      % in grid points
% % % %             Extract the wake expansion meadn(j,1) is at x = 0. So gdeltag+1
% %             exprad{i}.r(j) = meand(j,gdeltag+1);
% % % %             Calculate wake delays
% %             delays(j,i)    = gdeltag;
% %         else % Upwind turbines
% %             exprad{i}.r(j) = 0;
% %             delays(j,i)    = 0;
% %         end
% %     end
% % end

[~,nMast] = size(farm.posMast);
[~,nTurb] = size(farm.pos);
% % Extract wake expansion at turbines
for i = 1:length(posTot(1,:))
    for j = 1:length(posTot(1,:))         % Length from j to i
        gdelta = max(min(posTot(1,i)-posTot(1,j),grid.xsize),0);
        
        if(gdelta > 0)                   % For downwind turbines/masts
            gdeltag        = floor(gdelta/grid.size);      % in grid points
% %             Extract the wake expansion meadn(j,1) is at x = 0. So gdeltag+1
            if j <= nTurb
                exprad{i}.r(j) = meand(j,gdeltag+1);
            end
% %             Calculate wake delays
            delays(j,i)    = gdeltag;
        else % Upwind turbines
            if j <= nTurb
                exprad{i}.r(j) = 0;
            end
            delays(j,i)    = 0;
        end
    end
end



%% Setup wind field block
% % Name of block
wfname = [mdlname '/Wind Field'];

% % Add ambient data
set_param([wfname '/Ambient Field'],'UserData',wind);
set_param([wfname '/Ambient Field'],'UserDataPersistent','on');

% % Setup concatenators
% % Set output logging
set_param([wfname '/aconcat'],'NumInputs',Nstr);
h = get_param([wfname '/aconcat'],'PortHandles');
set_param(h.Outport(1),'DataLogging','on');
set_param(h.Outport(1),'DataLoggingName','deficit');
set_param(h.Outport(1),'DataLoggingNameMode','Custom');
    
set_param([wfname '/wcconcat'],'NumInputs',Nstr);
h = get_param([wfname '/wcconcat'],'PortHandles');
set_param(h.Outport(1),'DataLogging','on');
set_param(h.Outport(1),'DataLoggingName','wakecenter');
set_param(h.Outport(1),'DataLoggingNameMode','Custom');


% % Setup muxes
set_param([wfname '/CTdemux'],'Outputs',Nstr);
%set_param([wfname '/Enabledemux'],'Outputs',Nstr);
set_param([wfname '/Vrotmux'],'Inputs',Nstr);
set_param([wfname '/Vnacmux'],'Inputs',Nstr);
if nMast > 0
    set_param([wfname '/Vmastmux'],'Inputs',num2str(nMast));
else
    set_param([wfname '/Vmastmux'],'Inputs','1');
end

% % Setup sample time of rate transition blocks
set_param([wfname '/CTrate'],'OutPortSampleTime','wind.Ts');
%set_param([wfname '/ENrate'],'OutPortSampleTime',num2str(wind.Ts));

%% Add wake and wind speed blocks
y = 0;
x = 0;
% % Load wind library
load_system('libwind');
for i = 1:length(farm.turbines)
% %     Add wake block
    wtwk(i).block = add_block(['libwind/Wake effects'],[mdlname '/Wind Field/Wake' int2str(i)]);
    
% %     Set block position
    wkp = get_param(wtwk(i).block,'Position');
    set_param(wtwk(i).block,'Position',[wkp(1)+x wkp(2)+y,wkp(3)+x wkp(4)+y]);
    
% %     Unlink block for modification
    set_param(wtwk(i).block,'LinkStatus','none');
    
% %     Add wind speed block
    wtsp(i).block = add_block(['libwind/Wind speed'],[mdlname '/Wind Field/Wind Speed' int2str(i)]);
% %     Set block position
    wsp = get_param(wtsp(i).block,'Position');
    set_param(wtsp(i).block,'Position',[wsp(1)+x+300 wsp(2)+y,wsp(3)+x+300 wsp(4)+y]);
% %     Unlink block for modification
    set_param(wtsp(i).block,'LinkStatus','none');
    
        
% %     Add to y position
    y = y+175;
    
% %     Input wake
    add_line(wfname,['CTdemux/' int2str(i)],['Wake' int2str(i) '/2'],'Autorouting','on')
    add_line(wfname,'Ambient Field/2',['Wake' int2str(i) '/1'],'Autorouting','on')
    
% %     Input wind speed
    add_line(wfname,'aconcat/1',['Wind Speed' int2str(i) '/3'],'Autorouting','on')
    add_line(wfname,'wcconcat/1',['Wind Speed' int2str(i) '/2'],'Autorouting','on')
    add_line(wfname,'Ambient Field/1',['Wind Speed' int2str(i) '/1'],'Autorouting','on')
    
    
% %     Output wake
    add_line(wfname,['Wake' int2str(i) '/1'],['wcconcat/' int2str(i)],'Autorouting','on')
    add_line(wfname,['Wake' int2str(i) '/2'],['aconcat/' int2str(i)],'Autorouting','on')
    
% %     Output wind speed
    add_line(wfname,['Wind Speed' int2str(i) '/1'],['Vrotmux/' int2str(i)],'Autorouting','on')
    add_line(wfname,['Wind Speed' int2str(i) '/2'],['Vnacmux/' int2str(i)],'Autorouting','on')

    
% %     Configure wake
    clear wp;
    
    wp.farm.x      = farmx;
    wp.farm.xTot   = posTot(1,:);
    wp.grid        = grid;
    wp.wt.x        = farm.pos(1,i);
    wp.wt.y        = farm.pos(2,i);
    wp.wt.rotor.r  = wtpar(i).rotor.radius;
    wp.wt.num      = i;
    wp.wt.exprad   = exprad{i}.r;
    wp.wt.meand.r  = meand(i,:);
    wp.ts          = wind.Ts;
    
% %     Set wake parameters
    set_param(wtwk(i).block,'UserData',wp);
    set_param(wtwk(i).block,'UserDataPersistent','on');
    
% %     Set wind speed parameters
    set_param(wtsp(i).block,'UserData',wp);
    set_param(wtsp(i).block,'UserDataPersistent','on');
    

% %     Insert delays (wind travelling speed) 
% %     Name of wake block
    wkstr          = [wfname '/Wake' int2str(i)];
% %     Setup muxes
% % % % % % % % % %     set_param([wkstr '/ademux'],'Outputs',Nstr);
% % % % % % % % % %     set_param([wkstr '/wcdemux'],'Outputs',Nstr);
% % % % % % % % % %     set_param([wkstr '/amux'],'Inputs',Nstr);
% % % % % % % % % %     set_param([wkstr '/wcmux'],'Inputs',Nstr);
    set_param([wkstr '/ademux'],'Outputs',int2str(length(delays(i,:))));
    set_param([wkstr '/wcdemux'],'Outputs',int2str(length(delays(i,:))));
    set_param([wkstr '/amux'],'Inputs',int2str(length(delays(i,:))));
    set_param([wkstr '/wcmux'],'Inputs',int2str(length(delays(i,:))));
    
% %     Load simulink to get delay block
    load_system('simulink')
    ydel           = 0;
    xdel           = 0;
    ydel0          = 0;
    xdel0          = 500;
    for j = 1:length(delays(i,:))
% %         If downwind
        if(delays(i,j)>0)
% %             CT delay block
            ctd.block = add_block('simulink/Discrete/Integer Delay',[wfname '/Wake' int2str(i) '/aDelay' int2str(j)], ...
                        'MakeNameUnique','on' ,...
                          'samptime','wind.Ts', ...
                          'NumDelays',int2str(delays(i,j)), ...
                          'vinit','1');
% %             Set block position
            delPos = get_param(ctd.block,'Position');
            set_param(ctd.block,'Position',[xdel0+xdel ydel0+ydel,xdel0+delPos(3)-delPos(1)+xdel ydel0+delPos(4)-delPos(2)+ydel]);
            set_param(ctd.block,'LinkStatus','none');
                      
% %             Wake center delay block          
            wcd.block = add_block('simulink/Discrete/Integer Delay',[wfname '/Wake' int2str(i) '/WCDelay' int2str(j)], ...
                          'MakeNameUnique','on' ,...
                          'samptime','wind.Ts', ...
                          'NumDelays',int2str(delays(i,j)), ...
                          'vinit',int2str(farm.pos(2,i)));
% %             Set block position
            delPos = get_param(wcd.block,'Position');
            set_param(wcd.block,'Position',[xdel0+xdel ydel0+ydel+60,xdel0+delPos(3)-delPos(1)+xdel ydel0+delPos(4)-delPos(2)+ydel+60]);
            set_param(wcd.block,'LinkStatus','none');
            
         
% %             Connect to muxes
            add_line([wfname '/Wake' int2str(i)],['ademux/' int2str(j)],['aDelay' int2str(j) '/1'],'Autorouting','on')
            add_line([wfname '/Wake' int2str(i)],['aDelay' int2str(j) '/1'],['amux/' int2str(j)],'Autorouting','on')
            add_line([wfname '/Wake' int2str(i)],['wcdemux/' int2str(j)],['WCDelay' int2str(j) '/1'],'Autorouting','on')
            add_line([wfname '/Wake' int2str(i)],['WCDelay' int2str(j) '/1'],['wcmux/' int2str(j)],'Autorouting','on')
            
        else
            %Hackish - avoids algebraic loop. No delay - Connect constant output 1 to muxes
            
            ter1.block = add_block('simulink/Sinks/Terminator',[wfname '/Wake' int2str(i) '/TermAde' int2str(j)], ...
                                   'MakeNameUnique','on');
            %Set block position
            delPos = get_param(ter1.block,'Position');
            set_param(ter1.block,'Position',[xdel0+xdel ydel0+ydel,xdel0+delPos(3)-delPos(1)+xdel ydel0+delPos(4)-delPos(2)+ydel]);
            set_param(ter1.block,'LinkStatus','none');

            ter2.block = add_block('simulink/Sinks/Terminator',[wfname '/Wake' int2str(i) '/TermWcde' int2str(j)], ...
                                   'MakeNameUnique','on');
            %Set block position
            delPos = get_param(ter2.block,'Position');
            set_param(ter2.block,'Position',[xdel0+xdel ydel0+ydel+60,xdel0+delPos(3)-delPos(1)+xdel ydel0+delPos(4)-delPos(2)+ydel+60]);
            set_param(ter2.block,'LinkStatus','none');

            add_line([wfname '/Wake' int2str(i)],['ademux/' int2str(j)],['TermAde' int2str(j) '/1'],'Autorouting','on')
            add_line([wfname '/Wake' int2str(i)],['wcdemux/' int2str(j)],['TermWcde' int2str(j) '/1'],'Autorouting','on') 
            add_line([wfname '/Wake' int2str(i)],'Nodelay/1',['amux/' int2str(j)],'Autorouting','on')
            add_line([wfname '/Wake' int2str(i)],'Nodelay/1',['wcmux/' int2str(j)],'Autorouting','on')
        end
        ydel = ydel + 120;
    end
end

%% Add the wind speed blocks of the measurement masts
for i = 1:nMast
% %     Add wind speed block
    wtsp(i).block = add_block(['libwind/Measurement mast'],[mdlname '/Wind Field/Measurement mast' int2str(i)]);
% %     Set block position
    wsp = get_param(wtsp(i).block,'Position');
    set_param(wtsp(i).block,'Position',[wsp(1)+x+300 wsp(2)+y,wsp(3)+x+300 wsp(4)+y]);
% %     Unlink block for modification
    set_param(wtsp(i).block,'LinkStatus','none');
    
% %     Input wind speed
    add_line(wfname,'aconcat/1',['Measurement mast' int2str(i) '/3'],'Autorouting','on')
    add_line(wfname,'wcconcat/1',['Measurement mast' int2str(i) '/2'],'Autorouting','on')
    add_line(wfname,'Ambient Field/1',['Measurement mast' int2str(i) '/1'],'Autorouting','on')

% %     Output wind speed
    add_line(wfname,['Measurement mast' int2str(i) '/1'],['Vmastmux/' int2str(i)],'Autorouting','on')

% %     Add to y position
    y = y+175;

    clear wp;
    wp.farm.x      = farm.posMast(1,:);
    wp.grid        = grid;
    wp.mast.x      = farm.posMast(1,i);
    wp.mast.y      = farm.posMast(2,i);
    wp.mast.num    = length(farmx)+i;
    wp.ts          = wind.Ts;

    wp.farm.x      = farmx;
    wp.grid        = grid;
    wp.wt.rotor.r  = wtpar(1).rotor.radius;
    wp.wt.exprad   = exprad{length(farmx)+i}.r;
    wp.ts          = wind.Ts;
    
% %     Set wind speed parameters
    set_param(wtsp(i).block,'UserData',wp);
    set_param(wtsp(i).block,'UserDataPersistent','on');
    
end



% if(isunix)
%     add_block('simulink/Sources/Clock',[mdlname '/Clock']);
%     add_block('simulink/Sinks/Display',[mdlname '/Display']);
%     set_param([mdlname '/Display'],'Position',[230 545 370 575]);
%     set_param([mdlname '/Clock'],'Position',[40 550 60 570]);
%     add_line(mdlname,['Clock/1'],['Display/1'],'Autorouting','on');
%     set_param([mdlname '/Display'],'format','long');
%     %set_param([mdlname '/Display'],'position',get_param([mdlname '/Display'],'position')+[0 0 45 0]);
% end

%Add display when generated model is opened in unix
%Is not included yet as it needs debugging/validation
%set_param(mainsys,'PostLoadFcn',sprintf('if(isunix)\n\tmdlname = gcs;\n\tadd_block(''simulink/Sources/Clock'',[mdlname ''/Clock'']);\n\tadd_block(''simulink/Sinks/Display'',[mdlname ''/Display'']);\n\tset_param([mdlname ''/Display''],''Position'',[230 545 370 575]);\n\tset_param([mdlname ''/Clock''],''Position'',[40 550 60 570]);\n\tadd_line(mdlname,[''Clock/1''],[''Display/1''],''Autorouting'',''on'');\n\tset_param([mdlname ''/Display''],''format'',''long'');\nend'));
%set_param(mainsys,'CloseFcn',sprintf('if(isunix)\n\tmdlname = gcs;\n\tdelete_line(mdlname,[''Clock/1''],[''Display/1'']);\n\tdelete_block([mdlname ''/Clock'']);\n\tdelete_block([mdlname ''/Display'']);\nend'))

% %Setting Default parameters in simulation configuration

% if(isfield(wind,'SimTime'))
%     set_param(mainsys,'StopTime',num2str(wind.SimTime));
% else
%     set_param(mainsys,'StopTime',num2str((size(wind.Ux,2)-round(wind.grid.xsize/wind.grid.size))*wind.Ts));
% end
set_param(mainsys,'StopTime','tMax');
set_param(mainsys,'InlineParams','on');
set_param(mainsys,'EnhancedBackFolding','on')
set_param(mainsys,'SimulationMode','normal')



%Save model
save_system(mainsys);

bdclose(mdlname);

end
