function dc = parse_discon(discon_path)
    % PARSE_DISCON  Lee el fichero DISCON.IN de ROSCO y extrae los parámetros
    %               del controlador necesarios para FOWLTY.
    %
    % Formato DISCON.IN:  VALOR(ES)   ! KEYWORD   - descripcion
    % El keyword siempre aparece como primer token TRAS el símbolo !
    %
    % INPUT:
    %   discon_path : path completo al fichero DISCON.IN (string)
    %
    % OUTPUT:
    %   dc : struct con parámetros del controlador. Ver cabecera para lista completa.
    %
    % EJEMPLO:
    %   dc = parse_discon(sd.discon_file)
    
    base_dir = fileparts(discon_path);
    
    fid = fopen(discon_path, 'r');
    if fid == -1
        error('parse_discon: No se puede abrir: %s', discon_path);
    end
    
    % Keywords escalares: {keyword, campo}
    scalar_keys = {
        'VS_RtPwr',       'VS_RtPwr';
        'VS_RtTq',        'VS_RtTq';
        'VS_RefSpd',      'VS_RefSpd';
        'VS_MinOMSpd',    'VS_MinOMSpd';
        'VS_Rgn2K',       'VS_Rgn2K';
        'VS_TSRopt',      'VS_TSRopt';
        'VS_MaxTq',       'VS_MaxTq';
        'VS_MaxRat',      'VS_MaxRat';
        'VS_GenEff',      'VS_GenEff';
        'PC_RefSpd',      'PC_RefSpd';
        'PC_MaxPit',      'PC_MaxPit';
        'PC_MinPit',      'PC_MinPit';
        'PC_MaxRat',      'PC_MaxRat';
        'PC_MinRat',      'PC_MinRat';
        'PC_FinePit',     'PC_FinePit';
        'PC_GS_n',        'PC_GS_n';
        'PA_CornerFreq',  'PA_CornerFreq';
        'PA_Damping',     'PA_Damping';
        'WE_BladeRadius', 'WE_BladeRadius';
        'WE_GearboxRatio','WE_GearboxRatio';
        'WE_Jtot',        'WE_Jtot';
        'WE_RhoAir',      'WE_RhoAir';
    };
    
    % Inicializar struct
    dc = struct();
    for i = 1:size(scalar_keys,1)
        dc.(scalar_keys{i,2}) = NaN;
    end
    dc.PC_GS_angles    = [];
    dc.PC_GS_KP        = [];
    dc.PC_GS_KI        = [];
    dc.perf_file       = '';
    dc.perf_table_size = [0 0];
    
    kw_map = containers.Map(scalar_keys(:,1), scalar_keys(:,2));
    
    while ~feof(fid)
        raw_line = fgetl(fid);
        if ~ischar(raw_line), continue; end
    
        % --- Fichero de performance (entre comillas, keyword PerfFileName) ---
        % Formato: "path/al/fichero"   ! PerfFileName - descripcion
        if isempty(dc.perf_file)
            tok = regexp(raw_line, '"([^"]+)"\s*!?\s*PerfFileName', 'tokens');
            if ~isempty(tok)
                rel = strtrim(tok{1}{1});
                if ~strcmpi(rel, 'unused')
                    dc.perf_file = resolve_path(fullfile(base_dir, rel));
                end
            end
        end
    
        % --- Separar parte izquierda (valores) y parte derecha (keyword+desc) ---
        % El ! separa los valores del keyword
        bang_idx = strfind(raw_line, '!');
        if isempty(bang_idx), continue; end
    
        left  = strtrim(raw_line(1:bang_idx(1)-1));      % valores
        right = strtrim(raw_line(bang_idx(1)+1:end));     % keyword - descripcion
    
        if isempty(left) || isempty(right), continue; end
    
        % Extraer keyword: primer token de la parte derecha
        right_parts = strsplit(right);
        right_parts = right_parts(~cellfun(@isempty, right_parts));
        if isempty(right_parts), continue; end
        keyword = right_parts{1};
    
        % Extraer valores de la parte izquierda
        left_parts = strsplit(left);
        left_parts = left_parts(~cellfun(@isempty, left_parts));
        if isempty(left_parts), continue; end
    
        % --- PerfTableSize: dos valores ---
        if strcmp(keyword, 'PerfTableSize')
            if length(left_parts) >= 2
                v1 = str2double(left_parts{1});
                v2 = str2double(left_parts{2});
                if ~isnan(v1) && ~isnan(v2)
                    dc.perf_table_size = [v1, v2];
                end
            end
            continue;
        end
    
        % --- Vectores de gain scheduling (múltiples valores) ---
        if strcmp(keyword, 'PC_GS_angles')
            dc.PC_GS_angles = parse_vector(left_parts);
            continue;
        end
        if strcmp(keyword, 'PC_GS_KP')
            dc.PC_GS_KP = parse_vector(left_parts);
            continue;
        end
        if strcmp(keyword, 'PC_GS_KI')
            dc.PC_GS_KI = parse_vector(left_parts);
            continue;
        end
    
        % --- Escalares simples ---
        if isKey(kw_map, keyword) && isnan(dc.(kw_map(keyword)))
            val = str2double(left_parts{1});
            if ~isnan(val)
                dc.(kw_map(keyword)) = val;
            end
        end
    end
    
    fclose(fid);
    
    % --- Verificar esenciales ---
    essential = {'VS_RtPwr','VS_RtTq','VS_RefSpd','VS_TSRopt','PC_MaxRat'};
    all_ok = true;
    for i = 1:length(essential)
        if isnan(dc.(essential{i}))
            warning('parse_discon: No se encontró: %s', essential{i});
            all_ok = false;
        end
    end
    if isempty(dc.perf_file)
        warning('parse_discon: No se encontró PerfFileName');
        all_ok = false;
    end
    if isempty(dc.PC_GS_angles)
        warning('parse_discon: No se encontraron PC_GS_angles');
        all_ok = false;
    end
    
    % --- Resumen ---
    if all_ok
        fprintf('parse_discon: OK\n');
    else
        fprintf('parse_discon: completado con warnings\n');
    end
    fprintf('  Rated power    : %.3e W\n',    dc.VS_RtPwr);
    fprintf('  Rated torque   : %.3e N·m\n',  dc.VS_RtTq);
    fprintf('  Rated gen spd  : %.4f rad/s\n',dc.VS_RefSpd);
    fprintf('  Min gen spd    : %.4f rad/s\n',dc.VS_MinOMSpd);
    fprintf('  TSR optimal    : %.2f\n',       dc.VS_TSRopt);
    fprintf('  Rgn2K          : %.4e\n',       dc.VS_Rgn2K);
    fprintf('  PC GS points   : %d\n',         length(dc.PC_GS_angles));
    fprintf('  Pitch max rate : %.4f rad/s\n', dc.PC_MaxRat);
    fprintf('  Actuator freq  : %.4f rad/s\n', dc.PA_CornerFreq);
    fprintf('  Blade radius   : %.3f m\n',     dc.WE_BladeRadius);
    fprintf('  Drivetrain J   : %.3e kg·m²\n', dc.WE_Jtot);
    fprintf('  Perf file      : %s\n',         dc.perf_file);
    fprintf('  Perf table size: %d x %d\n',    dc.perf_table_size(1), dc.perf_table_size(2));
end

% -------------------------------------------------------------------------
function vec = parse_vector(tokens)
    vec = cellfun(@str2double, tokens);
    vec = vec(~isnan(vec));
end

% -------------------------------------------------------------------------
function out = resolve_path(p)
    parts = strsplit(strrep(p, '\', '/'), '/');
    resolved = {};
    for i = 1:length(parts)
        if strcmp(parts{i}, '..')
            if ~isempty(resolved)
                resolved(end) = [];
            end
        elseif strcmp(parts{i}, '.') || isempty(parts{i})
        else
            resolved{end+1} = parts{i}; %#ok<AGROW>
        end
    end
    if ~isempty(p) && (p(1) == '/' || p(1) == '\')
        out = ['/' strjoin(resolved, '/')];
    else
        out = strjoin(resolved, filesep);
    end
end