function af = parse_airfoil(af_path)
% PARSE_AIRFOIL  Lee un fichero de perfil aerodinámico en formato
%                AirfoilInfo (usado por AeroDyn v15) y extrae las tablas
%                de Cl y Cd en función del ángulo de ataque.
%
% INPUT:
%   af_path : path completo al fichero polar .dat (string)
%
% OUTPUT:
%   af      : struct con los campos:
%               .alpha - Ángulo de ataque [deg]  (N x 1)
%               .cl    - Coeficiente de sustentación [-] (N x 1)
%               .cd    - Coeficiente de resistencia [-]  (N x 1)
%               .path  - Path del fichero leído
%
% NOTA: Si el fichero tiene múltiples tablas (NumTabs > 1) se usa la
%       primera tabla únicamente (Re nominal).
%
% EJEMPLO:
%   af = parse_airfoil(ad.airfoil_files{1})

fid = fopen(af_path, 'r');
if fid == -1
    error('parse_airfoil: No se puede abrir: %s', af_path);
end

af.alpha = [];
af.cl    = [];
af.cd    = [];
af.path  = af_path;

num_alpha   = 0;
in_table    = false;
alpha_count = 0;
alpha_arr   = [];
cl_arr      = [];
cd_arr      = [];

while ~feof(fid)
    raw_line = fgetl(fid);
    if ~ischar(raw_line), continue; end

    % --- Detección de cabecera de tabla ANTES de limpiar comentarios ---
    % La cabecera es: "!    Alpha      Cl      Cd        Cm"
    % Puede empezar con ! pero contiene Alpha y Cl
    if ~in_table && contains(raw_line, 'Alpha') && contains(raw_line, 'Cl')
        in_table = true;
        continue;
    end

    % --- Leer filas de datos si estamos dentro de la tabla ---
    if in_table && alpha_count < num_alpha
        % Las filas de datos son números, no comentarios
        line = strtrim(raw_line);
        if isempty(line), continue; end

        parts = strsplit(line);
        parts = parts(~cellfun(@isempty, parts));
        if length(parts) < 3, continue; end

        v1 = str2double(parts{1});
        v2 = str2double(parts{2});
        v3 = str2double(parts{3});
        if isnan(v1) || isnan(v2) || isnan(v3), continue; end

        alpha_count = alpha_count + 1;
        alpha_arr(alpha_count) = v1; %#ok<AGROW>
        cl_arr(alpha_count)    = v2; %#ok<AGROW>
        cd_arr(alpha_count)    = v3; %#ok<AGROW>

        if num_alpha > 0 && alpha_count == num_alpha
            break; % tabla completa, ignorar tablas adicionales
        end
        continue;
    end

    % --- Parseo de keywords (solo si no estamos leyendo tabla) ---
    % Quitar comentarios
    line = regexprep(raw_line, '[!#].*$', '');
    line = strtrim(line);
    if isempty(line), continue; end

    parts = strsplit(line);
    parts = parts(~cellfun(@isempty, parts));
    if length(parts) < 2, continue; end

    % Número de puntos en la tabla
    if num_alpha == 0 && strcmp(parts{2}, 'NumAlf')
        num_alpha = str2double(parts{1});
        alpha_arr = zeros(num_alpha, 1);
        cl_arr    = zeros(num_alpha, 1);
        cd_arr    = zeros(num_alpha, 1);
    end
end

fclose(fid);

if alpha_count == 0
    warning('parse_airfoil: No se encontraron datos en: %s', af_path);
    return;
end

if num_alpha > 0 && alpha_count ~= num_alpha
    warning('parse_airfoil: Se esperaban %d puntos, se leyeron %d en: %s', ...
        num_alpha, alpha_count, af_path);
end

af.alpha = alpha_arr(1:alpha_count);
af.cl    = cl_arr(1:alpha_count);
af.cd    = cd_arr(1:alpha_count);

end