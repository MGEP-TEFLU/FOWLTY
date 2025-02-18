function create_windfarm()
% Interactive script for creating a wind farm model
%
% Copyright 2009 - Aalborg University
% Author Jacob Grunnet
% Modified Martin Kragelund
% Modified Yerai Peña-Sanchez
%

%% Ask if a wind farm wants to be created
go                      = questdlg('Create a new wind farm model?','Wind farm creator','Continue','Cancel','Continue');
if strcmp(go,'Cancel')
    return
end

%% Select the turbine model
load_system('libturbines');
str1                    = find_system('libturbines');
str                     = cell(length(str1)-1,1);
for i = 2:length(str1)
 st                     = str1{i};
 str{i-1}               = st(length('libturbines/')+1:end);
end
[s,v]                   = listdlg('PromptString','Select a turbine model',...
                                  'SelectionMode','single',...
                                  'ListString',str);
if(~v)
    return
end

parm                    = get_param(['libturbines/' str{s}],'UserData');
rotorRadius             = parm.public.rotor.radius;
turbine                 = str{s};

% Number of rows and columns of turbines. Minimum 4 turbines in row
num                     = 3;    % Default value
nGrid                   = max(ceil(sqrt(num)),4);
turbDist                = 400;
[yvec,xvec]             = meshgrid(0:(nGrid-1),0:(nGrid-1));
xvec                    = 100 + xvec(1:num)*turbDist;
yvec                    = yvec(1:num)*turbDist + 100;
pos                     = [xvec(:) yvec(:)].';

numMast                 = 1;
posMast                 = [0 mean(yvec)].';

%% Select number of turbines/measuring posts and their location
options.Resize          = 'on';
options.WindowStyle     = 'normal';
options.Interpreter     = 'none';

while(true)
    prompt              = {'Enter number of turbines','Enter turbine position matrix (x1 y1 ; x2 y2 ; ... ; xn yn)',...
                           'Enter number of wind measuring masts','Enter masts position matrix (not used if number of masts is 0)'};
    dlgtitle            = 'Wind Farm Creator';
    dims                = [1 70];
    definput            = {mat2str(num), mat2str(pos.'), mat2str(numMast), mat2str(posMast.')};
    answer              = inputdlg(prompt,dlgtitle,dims,definput);
    if(isempty(answer))
        return
    end
    num                 = str2double(answer{1});
    pos                 = str2num(answer{2}).';
    numMast             = str2double(answer{3});
    posMast             = str2num(answer{4}).';

    % Check for errors
    if(num ~= length(pos(1,:)))
        err             = errordlg('Specified turbine number and locations do not coincide.','Wind Farm Creator');
        uiwait(err);
        nGrid           = max(ceil(sqrt(num)),4);
        turbDist        = 400;
        [yvec,xvec]     = meshgrid(0:(nGrid-1),0:(nGrid-1));
        xvec            = 100 + xvec(1:num)*turbDist;
        yvec            = yvec(1:num)*turbDist + 100;
        pos             = [xvec(:) yvec(:)]';
    elseif(isempty(num) || mod(num,1) ~= 0 || num < 0)
        err=errordlg('Not a valid turbine number.','Wind Farm Creator');
        uiwait(err);
    elseif(isempty(pos) || ~isempty(find((size(pos) == [2,num]) == 0,1)) || ~isempty(find(pos < 0,1)) || ~isempty(find(pos(2,:) < rotorRadius)))
        err=errordlg('Not a valid position matrix. Ensure that all turbines have positive coordinates and y position is at least equal to rotor radius.','Wind Farm Creator');
        uiwait(err);
        nGrid           = max(ceil(sqrt(num)),4);
        turbDist        = 400;
        [yvec,xvec]     = meshgrid(0:(nGrid-1),0:(nGrid-1));
        xvec            = 100 + xvec(1:num)*turbDist;
        yvec            = yvec(1:num)*turbDist + 100;
        pos             = [xvec(:) yvec(:)]';
    elseif(isempty(numMast) || mod(numMast,1) ~= 0 || numMast < 0)
        err=errordlg('Not a valid wind measurement mast number.','Wind Farm Creator');
        uiwait(err);
    elseif(numMast~=0 && (isempty(posMast) || ~isempty(find((size(posMast) == [2,numMast]) == 0,1)) || ~isempty(find(posMast < 0,1))))
        err=errordlg('Not a valid position matrix. Ensure that all measuring masts have positive coordinates.','Wind Farm Creator');
        uiwait(err);
        posMast         = [zeros(numMast,1) linspace(min(yvec),max(yvec),numMast).'].';
    else
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Put graph with the turbines in order to make sure it is what the user wants (put the break inside another if) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        break
    end
end

if numMast == 0
    posMast             = [];
end
%% Create the wind field
[wind,wfile,wpath] = create_windfield([pos posMast]);

if(isempty(wind))
    [wfile,wpath,wfidx] = uigetfile('*.mat','Select a Wind Data File');

    if isequal(wfile,0) || isequal(wpath,0)
       % no dice then
       return;
    end

    % Loads wind profile data
    l                   = load([wpath wfile]);
    wind                = l.wind;
end

[file,path]             = uiputfile('*.mdl','Save Wind Farm Model As',[wpath 'windfarm.mdl']);

if isequal(file,0) || isequal(path,0)
   return
end

if(exist([path file]))
    disp('file existed and is removed')
    delete([path file])
end

%% Create data structures for mdlcreate
for i = 1:num
    farm.turbines{i}    = turbine;
end

farm.pos                = pos;
farm.posMast            = posMast;

% close system stuff
gen_windfarm(file,path,farm,wind);

%% Create an init function at the desired path if it doesn't exist
initLoc                 = which('init_.m');

S                       = readlines(initLoc);
S(11,:)                 = sprintf("SimulinkModel               = '%s';",file(1:end-4));
S(12,:)                 = sprintf("WindCase                    = '%s';",wfile(1:end-4));

if strcmp(turbine(end-5:end),'faulty') ~= 1
    S(33:35,:)          = [];
    S(33,:)             = sprintf('faultSignal                 = [tsim zeros(length(tsim),nT)];');
end

if exist([path 'init_' file(1:end-2)]) == 0
    fileID              = fopen([path 'init_' file(1:end-2)],'w');
    fwrite(fileID, strjoin(S, '\n'));
    fclose(fileID);
end

if exist([path 'faultScenarios.m']) == 0 && strcmp(turbine(end-5:end),'faulty')
    copyfile([initLoc(1:end-7) 'faultScenarios_original.m'],[path 'faultScenarios.m']);
end

if exist([path 'windRealisationGenerator.m']) == 0
    S                   = readlines([initLoc(1:end-7) 'windRealisationGenerator_original.m']);
    S(10,:)             = sprintf("originalName        = '%s';",wfile(1:end-4));
    fileID              = fopen([path 'windRealisationGenerator.m'],'w');
    fwrite(fileID, strjoin(S, '\n'));
    fclose(fileID);
end
end