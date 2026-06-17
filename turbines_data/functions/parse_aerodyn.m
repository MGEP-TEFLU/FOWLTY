function ad = parse_aerodyn(ad_path)
    % PARSE_AERODYN  Lee el fichero AeroDyn de OpenFAST y extrae la lista de
    %                ficheros de perfiles aerodinámicos y el fichero de pala.
    %
    % INPUT:
    %   ad_path : path completo al fichero AeroDyn .dat (string)
    %
    % OUTPUT:
    %   ad      : struct con los campos:
    %               .num_airfoils  - Número de perfiles [-]
    %               .airfoil_files - Cell array con paths a cada fichero polar
    %               .blade_file    - Path al fichero AeroDyn blade
    %               .rho           - Densidad del aire (si no es "default") [kg/m³]
    %
    % EJEMPLO:
    %   ad = parse_aerodyn(files.aerodyn)
    
    base_dir = fileparts(ad_path);
    
    fid = fopen(ad_path, 'r');
    if fid == -1
        error('parse_aerodyn: No se puede abrir el fichero: %s', ad_path);
    end
    
    ad.num_airfoils  = 0;
    ad.airfoil_files = {};
    ad.blade_file    = '';
    ad.rho           = NaN;
    
    reading_airfoils = false;
    airfoil_count    = 0;
    
    while ~feof(fid)
        raw_line = fgetl(fid);
        if ~ischar(raw_line), continue; end
    
        % --- Densidad del aire (puede ser "default") ---
        if isnan(ad.rho)
            tok = regexp(raw_line, '"default"\s+AirDens', 'tokens');
            if ~isempty(tok)
                ad.rho = NaN; % se cogerá del .fst
            else
                line = regexprep(raw_line, '\s+-\s+.*$', '');
                line = regexprep(line, '!.*$', '');
                parts = strsplit(strtrim(line));
                parts = parts(~cellfun(@isempty, parts));
                if length(parts) >= 2 && strcmp(parts{2}, 'AirDens')
                    val = str2double(parts{1});
                    if ~isnan(val)
                        ad.rho = val;
                    end
                end
            end
        end
    
        % --- Número de ficheros de perfiles ---
        if ad.num_airfoils == 0
            line = regexprep(raw_line, '\s+-\s+.*$', '');
            line = regexprep(line, '!.*$', '');
            parts = strsplit(strtrim(line));
            parts = parts(~cellfun(@isempty, parts));
            if length(parts) >= 2 && strcmp(parts{2}, 'NumAFfiles')
                ad.num_airfoils = str2double(parts{1});
                reading_airfoils = true;
                continue;
            end
        end
    
        % --- Leer paths de perfiles ---
        % Las líneas de perfiles son strings entre comillas
        % La primera tiene además el keyword AFNames, las demás solo el path
        if reading_airfoils && airfoil_count < ad.num_airfoils
            tok = regexp(raw_line, '"([^"]+)"', 'tokens');
            if ~isempty(tok)
                rel = strtrim(tok{1}{1});
                full_p = resolve_path(fullfile(base_dir, rel));
                airfoil_count = airfoil_count + 1;
                ad.airfoil_files{airfoil_count} = full_p;
                if airfoil_count == ad.num_airfoils
                    reading_airfoils = false;
                end
                continue;
            end
        end
    
        % --- Blade file (ADBlFile(1)) ---
        if isempty(ad.blade_file)
            tok = regexp(raw_line, '"([^"]+)"\s+ADBlFile\(1\)', 'tokens');
            if ~isempty(tok)
                rel = strtrim(tok{1}{1});
                ad.blade_file = resolve_path(fullfile(base_dir, rel));
            end
        end
    end
    
    fclose(fid);
    
    % Verificar
    if ad.num_airfoils == 0
        warning('parse_aerodyn: No se encontró NumAFfiles');
    end
    if isempty(ad.airfoil_files)
        warning('parse_aerodyn: No se encontraron ficheros de perfiles');
    elseif length(ad.airfoil_files) ~= ad.num_airfoils
        warning('parse_aerodyn: Se esperaban %d perfiles, se encontraron %d', ...
            ad.num_airfoils, length(ad.airfoil_files));
    end
    if isempty(ad.blade_file)
        warning('parse_aerodyn: No se encontró ADBlFile');
    end
    
    fprintf('parse_aerodyn: OK\n');
    fprintf('  Num airfoils : %d\n',    ad.num_airfoils);
    fprintf('  Blade file   : %s\n',    ad.blade_file);
    fprintf('  Airfoil 1    : %s\n',    ad.airfoil_files{1});
    fprintf('  Airfoil end  : %s\n',    ad.airfoil_files{end});
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