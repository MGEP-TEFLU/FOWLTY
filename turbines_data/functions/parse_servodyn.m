function sd = parse_servodyn(sd_path)
    % PARSE_SERVODYN  Lee el fichero ServoDyn de OpenFAST y extrae los
    %                 parámetros de control y generador necesarios para FOWLTY.
    %
    % INPUT:
    %   sd_path : path completo al fichero ServoDyn .dat (string)
    %
    % OUTPUT:
    %   sd      : struct con los campos:
    %               .GenEff       - Eficiencia del generador [%]
    %               .GBRatio      - Ratio de la caja (confirmación vs ElastoDyn)
    %               .discon_file  - Path al fichero DISCON.IN del controlador
    %
    % EJEMPLO:
    %   sd = parse_servodyn(files.servodyn)
    
    base_dir = fileparts(sd_path);
    
    fid = fopen(sd_path, 'r');
    if fid == -1
        error('parse_servodyn: No se puede abrir el fichero: %s', sd_path);
    end
    
    % Keywords escalares
    scalar_keys = {
        'GenEff',  'GenEff';
    };
    
    sd = struct();
    for i = 1:size(scalar_keys,1)
        sd.(scalar_keys{i,2}) = NaN;
    end
    sd.discon_file = '';
    
    kw_map = containers.Map(scalar_keys(:,1), scalar_keys(:,2));
    
    while ~feof(fid)
        raw_line = fgetl(fid);
        if ~ischar(raw_line), continue; end
    
        % Path al DISCON.IN (entre comillas, seguido de DLL_InFile)
        if isempty(sd.discon_file)
            tok = regexp(raw_line, '"([^"]+)"\s+DLL_InFile', 'tokens');
            if ~isempty(tok)
                rel = strtrim(tok{1}{1});
                if ~strcmpi(rel, 'unused')
                    sd.discon_file = resolve_path(fullfile(base_dir, rel));
                end
            end
        end
    
        % Escalares
        line = regexprep(raw_line, '\s+-\s+.*$', '');
        line = regexprep(line, '!.*$', '');
        line = strtrim(line);
        if isempty(line), continue; end
    
        parts = strsplit(line);
        parts = parts(~cellfun(@isempty, parts));
        if length(parts) < 2, continue; end
    
        keyword = parts{2};
        if isKey(kw_map, keyword) && isnan(sd.(kw_map(keyword)))
            sd.(kw_map(keyword)) = str2double(parts{1});
        end
    end
    
    fclose(fid);
    
    % Verificar
    if isnan(sd.GenEff)
        warning('parse_servodyn: No se encontró GenEff');
    end
    if isempty(sd.discon_file)
        warning('parse_servodyn: No se encontró DLL_InFile (DISCON.IN)');
    end
    
    fprintf('parse_servodyn: OK\n');
    fprintf('  Gen efficiency : %.3f %%\n', sd.GenEff);
    fprintf('  DISCON file    : %s\n',      sd.discon_file);

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