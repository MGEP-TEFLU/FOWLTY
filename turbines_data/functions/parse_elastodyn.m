function ed = parse_elastodyn(ed_path)
    % PARSE_ELASTODYN  Lee el fichero ElastoDyn de OpenFAST y extrae los
    %                  parámetros estructurales necesarios para FOWLTY.
    %
    % INPUT:
    %   ed_path : path completo al fichero ElastoDyn .dat (string)
    %
    % OUTPUT:
    %   ed      : struct con campos de parámetros estructurales
    %
    % EJEMPLO:
    %   ed = parse_elastodyn(files.elastodyn)
    
    % Directorio base para resolver paths relativos
    base_dir = fileparts(ed_path);
    
    % Abrir fichero
    fid = fopen(ed_path, 'r');
    if fid == -1
        error('parse_elastodyn: No se puede abrir el fichero: %s', ed_path);
    end
    
    % Keywords escalares que buscamos: {keyword_en_fichero, campo_en_struct}
    scalar_keys = {
        'TipRad',    'TipRad';
        'HubRad',    'HubRad';
        'HubMass',   'HubMass';
        'HubIner',   'HubIner';
        'GenIner',   'GenIner';
        'NacMass',   'NacMass';
        'YawBrMass', 'YawBrMass';
        'GBRatio',   'GBRatio';
        'GBoxEff',   'GBoxEff';
        'DTTorSpr',  'DTTorSpr';
        'DTTorDmp',  'DTTorDmp';
        'RotSpeed',  'RotSpeed';
        'TowerHt',   'TowerHt';
        'TowerBsHt', 'TowerBsHt';
        'NumBl',     'NumBl';
        'OverHang',  'OverHang';
        'ShftTilt',  'ShftTilt';
        'Twr2Shft',  'Twr2Shft';
    };
    
    % Inicializar struct con NaN
    ed = struct();
    for i = 1:size(scalar_keys,1)
        ed.(scalar_keys{i,2}) = NaN;
    end
    ed.PreCone    = NaN;
    ed.blade_file = '';
    ed.tower_file = '';
    
    % Construir mapa keyword->campo para búsqueda rápida
    kw_map = containers.Map(scalar_keys(:,1), scalar_keys(:,2));
    
    % Leer línea a línea
    while ~feof(fid)
        raw_line = fgetl(fid);
        if ~ischar(raw_line), continue; end
        
        % --- Paths entre comillas (blade file, tower file) ---
        if isempty(ed.blade_file)
            tok = regexp(raw_line, '"([^"]+)"\s+BldFile1', 'tokens');
            if ~isempty(tok)
                ed.blade_file = resolve_path(fullfile(base_dir, strtrim(tok{1}{1})));
            end
        end
        if isempty(ed.tower_file)
            tok = regexp(raw_line, '"([^"]+)"\s+TwrFile', 'tokens');
            if ~isempty(tok)
                ed.tower_file = resolve_path(fullfile(base_dir, strtrim(tok{1}{1})));
            end
        end
        
        % --- Escalares: formato "VALOR   KEYWORD   - descripcion" ---
        % Quitar todo tras el primer " - " (descripción)
        line = regexprep(raw_line, '\s+-\s+.*$', '');
        % Quitar comentarios con !
        line = regexprep(line, '!.*$', '');
        line = strtrim(line);
        if isempty(line), continue; end
        
        % Split en tokens por espacios/tabs
        parts = strsplit(line);
        parts = parts(~cellfun(@isempty, parts));
        if length(parts) < 2, continue; end
        
        val_str = parts{1};
        keyword  = parts{2};
        
        % Limpiar keyword de posibles paréntesis tipo "PreCone(1)"
        kw_clean = regexprep(keyword, '\(.*\)', '');
        
        % ¿Es un keyword que buscamos?
        if isKey(kw_map, keyword) && isnan(ed.(kw_map(keyword)))
            ed.(kw_map(keyword)) = str2double(val_str);
            
        elseif isKey(kw_map, kw_clean) && isnan(ed.(kw_map(kw_clean)))
            ed.(kw_map(kw_clean)) = str2double(val_str);
            
        % PreCone: coger PreCone(1) solamente
        elseif strcmp(keyword, 'PreCone(1)') && isnan(ed.PreCone)
            ed.PreCone = str2double(val_str);
        end
    end
    
    fclose(fid);
    
    % --- Verificar campos esenciales ---
    essential = {'TipRad','HubRad','HubMass','HubIner','GenIner', ...
                 'NacMass','GBRatio','RotSpeed','TowerHt','NumBl'};
    all_ok = true;
    for i = 1:length(essential)
        if isnan(ed.(essential{i}))
            warning('parse_elastodyn: No se encontró el campo: %s', essential{i});
            all_ok = false;
        end
    end
    
    % --- Resumen ---
    if all_ok
        fprintf('parse_elastodyn: OK\n');
    else
        fprintf('parse_elastodyn: completado con warnings\n');
    end
    fprintf('  Rotor radius  : %.3f m\n',       ed.TipRad);
    fprintf('  Hub radius    : %.3f m\n',       ed.HubRad);
    fprintf('  Num blades    : %d\n',           ed.NumBl);
    fprintf('  Rotor speed   : %.2f rpm\n',     ed.RotSpeed);
    fprintf('  GB ratio      : %.1f\n',         ed.GBRatio);
    fprintf('  Tower height  : %.3f m\n',       ed.TowerHt);
    fprintf('  Hub mass      : %.0f kg\n',      ed.HubMass);
    fprintf('  Nacelle mass  : %.0f kg\n',      ed.NacMass);
    fprintf('  Gen inertia   : %.1f kg·m²\n',   ed.GenIner);
    fprintf('  DT spring     : %.3e N·m/rad\n', ed.DTTorSpr);
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
            % ignorar
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