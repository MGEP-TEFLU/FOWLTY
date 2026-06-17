function bl = parse_elastodyn_blade(blade_path, R_hub, R_tip)
% PARSE_ELASTODYN_BLADE  Lee el fichero ElastoDyn blade y calcula la masa
%                        total y la inercia de la pala.
%
% INPUT:
%   blade_path : path al fichero ElastoDyn blade .dat
%   R_hub      : radio del buje [m] (para calcular r absoluto)
%   R_tip      : radio total del rotor [m]
%
% OUTPUT:
%   bl         : struct con:
%                  .mass    - masa total de la pala [kg]
%                  .inertia - inercia respecto al eje del rotor [kg·m²]
%                  .BlFract - fracción de longitud de pala (0 a 1)
%                  .BMassDen - densidad lineal de masa [kg/m]
%
% La masa se calcula integrando BMassDen a lo largo del span:
%   mass = integral(BMassDen * dr)
% La inercia se calcula respecto al apex del rotor:
%   inertia = integral(BMassDen * r_abs^2 * dr)
% donde r_abs = R_hub + BlFract*(R_tip - R_hub)
%
% EJEMPLO:
%   bl = parse_elastodyn_blade(ed.blade_file, ed.HubRad, ed.TipRad)

fid = fopen(blade_path, 'r');
if fid == -1
    error('parse_elastodyn_blade: No se puede abrir: %s', blade_path);
end

bl.BlFract   = [];
bl.BMassDen  = [];
bl.mass      = NaN;
bl.inertia   = NaN;

num_stations = 0;
in_table     = false;
station_count = 0;

BlFract_arr  = [];
BMassDen_arr = [];

while ~feof(fid)
    raw_line = fgetl(fid);
    if ~ischar(raw_line), continue; end

    % Número de estaciones
    if num_stations == 0
        line = regexprep(raw_line, '\s+-\s+.*$', '');
        line = regexprep(line, '!.*$', '');
        parts = strsplit(strtrim(line));
        parts = parts(~cellfun(@isempty, parts));
        if length(parts) >= 2 && strcmp(parts{2}, 'NBlInpSt')
            num_stations = str2double(parts{1});
        end
    end

    % Detectar cabecera de tabla (contiene BlFract y BMassDen)
    if ~in_table && contains(raw_line, 'BlFract') && contains(raw_line, 'BMassDen')
        in_table = true;
        fgetl(fid); % saltar línea de unidades
        continue;
    end

    % Leer filas de la tabla
    if in_table && station_count < num_stations
        line = strtrim(raw_line);
        if isempty(line) || line(1) == '!', continue; end

        parts = strsplit(line);
        parts = parts(~cellfun(@isempty, parts));
        if length(parts) < 4, continue; end

        v1 = str2double(parts{1});  % BlFract
        v4 = str2double(parts{4});  % BMassDen
        if isnan(v1) || isnan(v4), continue; end

        station_count = station_count + 1;
        BlFract_arr(station_count)  = v1; %#ok<AGROW>
        BMassDen_arr(station_count) = v4; %#ok<AGROW>

        if station_count == num_stations
            break;
        end
    end
end

fclose(fid);

if station_count == 0
    warning('parse_elastodyn_blade: No se encontraron datos en: %s', blade_path);
    return;
end

bl.BlFract  = BlFract_arr(:);
bl.BMassDen = BMassDen_arr(:);

% Radio absoluto de cada estación
blade_length = R_tip - R_hub;
r_abs = R_hub + bl.BlFract * blade_length;  % [m] desde el centro del rotor

% Integrar con trapecios
bl.mass    = trapz(r_abs, bl.BMassDen);                    % [kg]
bl.inertia = trapz(r_abs, bl.BMassDen .* r_abs.^2);        % [kg·m²]

fprintf('parse_elastodyn_blade: OK\n');
fprintf('  Blade mass    : %.1f kg\n',    bl.mass);
fprintf('  Blade inertia : %.4e kg·m²\n', bl.inertia);

end