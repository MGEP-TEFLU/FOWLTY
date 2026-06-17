function files = parse_fst(fst_path)
% PARSE_FST  Lee el fichero principal .fst de OpenFAST y devuelve los paths
%            a los subfiles necesarios para FOWLTY.
%
% INPUT:
%   fst_path  : path completo al fichero .fst (string)
%
% OUTPUT:
%   files     : struct con los campos:
%                 .elastodyn  - path al fichero ElastoDyn
%                 .aerodyn    - path al fichero AeroDyn
%                 .servodyn   - path al fichero ServoDyn
%                 .base_dir   - directorio base del .fst
%
% EJEMPLO:
%   files = parse_fst('github_files/OpenFAST/IEA-15-240-RWT-Monopile/IEA-15-240-RWT-Monopile.fst')

% Directorio base para resolver paths relativos
files.base_dir = fileparts(fst_path);

% Abrir fichero
fid = fopen(fst_path, 'r');
if fid == -1
    error('parse_fst: No se puede abrir el fichero: %s', fst_path);
end

% Keywords que buscamos y campo destino
keywords = {'EDFile', 'AeroFile', 'ServoFile'};
fields   = {'elastodyn', 'aerodyn', 'servodyn'};

% Inicializar campos a empty
for i = 1:length(fields)
    files.(fields{i}) = '';
end

% Leer línea a línea
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), continue; end
    
    % Eliminar comentarios (todo lo que va después de !)
    comment_idx = strfind(line, '!');
    if ~isempty(comment_idx)
        line = line(1:comment_idx(1)-1);
    end
    
    % Buscar cada keyword
    for i = 1:length(keywords)
        if contains(line, keywords{i})
            % Extraer el string entre comillas
            quoted = regexp(line, '"([^"]*)"', 'tokens');
            if ~isempty(quoted)
                rel_path = strtrim(quoted{1}{1});
                % Resolver path relativo respecto al base_dir
                full_path = fullfile(files.base_dir, rel_path);
                % Normalizar separadores (por si acaso hay ../ etc.)
                full_path = resolve_path(full_path);
                files.(fields{i}) = full_path;
            end
        end
    end
end

fclose(fid);

% Verificar que hemos encontrado los ficheros esenciales
essential = {'elastodyn', 'aerodyn', 'servodyn'};
for i = 1:length(essential)
    if isempty(files.(essential{i}))
        warning('parse_fst: No se encontró el fichero para: %s', essential{i});
    else
        if ~exist(files.(essential{i}), 'file')
            warning('parse_fst: Fichero no existe en disco: %s', files.(essential{i}));
        end
    end
end

end

% -------------------------------------------------------------------------
function out = resolve_path(p)
% Resuelve ../ y ./ en un path sin necesitar que el fichero exista en disco
% (a diferencia de realpath/cd que requieren existencia)
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
% Reconstruir (preservar slash inicial si era path absoluto)
if p(1) == '/' || p(1) == '\'
    out = ['/' strjoin(resolved, '/')];
else
    out = strjoin(resolved, filesep);
end
end