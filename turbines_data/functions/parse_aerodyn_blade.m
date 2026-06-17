function blade = parse_aerodyn_blade(blade_path)
    % PARSE_AERODYN_BLADE  Lee el fichero de geometría de pala de AeroDyn y
    %                      extrae la distribución radial, cuerda, torsión y
    %                      el índice de perfil para cada nodo.
    %
    % INPUT:
    %   blade_path : path completo al fichero AeroDyn blade .dat (string)
    %
    % OUTPUT:
    %   blade      : struct con los campos:
    %                  .num_nodes - Número de nodos [-]
    %                  .span      - Posición radial de cada nodo [m]  (num_nodes x 1)
    %                  .chord     - Longitud de cuerda [m]            (num_nodes x 1)
    %                  .twist     - Ángulo de torsión [deg]           (num_nodes x 1)
    %                  .af_id     - Índice de perfil (1-based) [-]    (num_nodes x 1)
    %
    % NOTA: Se ignoran BlCrvAC, BlSwpAC, BlCrvAng (sweep y coning) ya que
    %       el BEM los trata como pala recta.
    %
    % EJEMPLO:
    %   blade = parse_aerodyn_blade(ad.blade_file)
    
    fid = fopen(blade_path, 'r');
    if fid == -1
        error('parse_aerodyn_blade: No se puede abrir: %s', blade_path);
    end
    
    blade.num_nodes = 0;
    blade.span      = [];
    blade.chord     = [];
    blade.twist     = [];
    blade.af_id     = [];
    
    % Buscar NumBlNds para saber cuántos nodos hay
    num_nodes = 0;
    while ~feof(fid)
        raw_line = fgetl(fid);
        if ~ischar(raw_line), continue; end
    
        line = regexprep(raw_line, '\s+-\s+.*$', '');
        line = regexprep(line, '!.*$', '');
        parts = strsplit(strtrim(line));
        parts = parts(~cellfun(@isempty, parts));
    
        if length(parts) >= 2 && strcmp(parts{2}, 'NumBlNds')
            num_nodes = str2double(parts{1});
            break;
        end
    end
    
    if num_nodes == 0
        fclose(fid);
        error('parse_aerodyn_blade: No se encontró NumBlNds');
    end
    
    % Buscar la cabecera de la tabla de nodos
    % La tabla empieza después de las dos líneas de cabecera:
    %   BlSpn  BlCrvAC  BlSwpAC  BlCrvAng  BlTwist  BlChord  BlAFID ...
    %   (m)    (m)      (m)      (deg)     (deg)    (m)      (-)    ...
    header_found = false;
    while ~feof(fid)
        raw_line = fgetl(fid);
        if ~ischar(raw_line), continue; end
        % La línea de cabecera contiene 'BlSpn'
        if contains(raw_line, 'BlSpn')
            % Leer la siguiente línea de unidades y ya estamos listos
            fgetl(fid); % línea de unidades (m) (m) etc.
            header_found = true;
            break;
        end
    end
    
    if ~header_found
        fclose(fid);
        error('parse_aerodyn_blade: No se encontró la cabecera de la tabla de nodos');
    end
    
    % Leer los nodos
    % Formato de cada línea:
    % BlSpn  BlCrvAC  BlSwpAC  BlCrvAng  BlTwist  BlChord  BlAFID  BlCb  BlCenBn  BlCenBt
    % col 1  col 2    col 3    col 4     col 5    col 6    col 7   ...
    span  = zeros(num_nodes, 1);
    chord = zeros(num_nodes, 1);
    twist = zeros(num_nodes, 1);
    af_id = zeros(num_nodes, 1);
    
    node = 0;
    while ~feof(fid) && node < num_nodes
        raw_line = fgetl(fid);
        if ~ischar(raw_line), continue; end
    
        % Ignorar líneas vacías y comentarios
        line = regexprep(raw_line, '!.*$', '');
        line = strtrim(line);
        if isempty(line), continue; end
    
        parts = strsplit(line);
        parts = parts(~cellfun(@isempty, parts));
    
        % Necesitamos al menos 7 columnas
        if length(parts) < 7, continue; end
    
        % Intentar parsear la primera columna como número
        val = str2double(parts{1});
        if isnan(val), continue; end
    
        node = node + 1;
        span(node)  = str2double(parts{1}); % BlSpn   col 1
        % parts{2} = BlCrvAC  (ignorado)
        % parts{3} = BlSwpAC  (ignorado)
        % parts{4} = BlCrvAng (ignorado)
        twist(node) = str2double(parts{5}); % BlTwist col 5
        chord(node) = str2double(parts{6}); % BlChord col 6
        af_id(node) = str2double(parts{7}); % BlAFID  col 7
    end
    
    fclose(fid);
    
    if node ~= num_nodes
        warning('parse_aerodyn_blade: Se esperaban %d nodos, se leyeron %d', ...
            num_nodes, node);
    end
    
    blade.num_nodes = node;
    blade.span      = span(1:node);
    blade.chord     = chord(1:node);
    blade.twist     = twist(1:node);
    blade.af_id     = af_id(1:node);
    
    fprintf('parse_aerodyn_blade: OK\n');
    fprintf('  Num nodes : %d\n',         blade.num_nodes);
    fprintf('  Span range: %.2f - %.2f m\n', blade.span(1), blade.span(end));
    fprintf('  Chord root: %.3f m\n',     blade.chord(1));
    fprintf('  Chord tip : %.3f m\n',     blade.chord(end));
    fprintf('  Twist root: %.3f deg\n',   blade.twist(1));
    fprintf('  Twist tip : %.3f deg\n',   blade.twist(end));
    fprintf('  AF IDs    : %d to %d\n',   min(blade.af_id), max(blade.af_id));
end