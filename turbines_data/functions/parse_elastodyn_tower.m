function tw = parse_elastodyn_tower(tower_path, tower_height)
% PARSE_ELASTODYN_TOWER  Lee el fichero ElastoDyn tower y calcula la masa
%                        total y la frecuencia natural del primer modo FA.
%
% INPUT:
%   tower_path   : path al fichero ElastoDyn tower .dat
%   tower_height : altura total de la torre [m] (de ElastoDyn principal)
%
% OUTPUT:
%   tw : struct con:
%          .mass     - masa total de la torre [kg]
%          .eigfreq  - frecuencia natural 1er modo FA amortiguada [rad/s]
%          .damp     - amortiguamiento estructural del 1er modo FA [-]
%          .HtFract  - fracción de altura (0 a 1)
%          .TMassDen - densidad lineal de masa [kg/m]
%          .TwFAStif - rigidez FA [N·m²]
%          .height   - altura total de la torre [m]
%
% La frecuencia natural se calcula con el cociente de Rayleigh usando
% los mode shapes del fichero y la distribución de masa y rigidez.
%
% EJEMPLO:
%   tw = parse_elastodyn_tower(ed.tower_file, ed.TowerHt)

fid = fopen(tower_path, 'r');
if fid == -1
    error('parse_elastodyn_tower: No se puede abrir: %s', tower_path);
end

tw.HtFract  = [];
tw.TMassDen = [];
tw.TwFAStif = [];
tw.mass     = NaN;
tw.eigfreq  = NaN;
tw.damp     = NaN;
tw.height   = NaN;

if nargin >= 2 && ~isnan(tower_height)
    tw.height = tower_height;
end

num_stations  = 0;
in_table      = false;
station_count = 0;
HtFract_arr   = [];
TMassDen_arr  = [];
TwFAStif_arr  = [];

% Mode shape coefficients (polynomial: φ = c2*x² + c3*x³ + c4*x⁴ + c5*x⁵ + c6*x⁶)
% where x = HtFract (0 to 1)
FA_mode1_sh = zeros(1,5);  % [c2 c3 c4 c5 c6]
damp_pct    = 1.0;         % default 1%

while ~feof(fid)
    raw_line = fgetl(fid);
    if ~ischar(raw_line), continue; end

    % Damping
    line = regexprep(raw_line, '\s+-\s+.*$', '');
    line = regexprep(line, '!.*$', '');
    parts = strsplit(strtrim(line));
    parts = parts(~cellfun(@isempty, parts));

    if length(parts) >= 2
        if strcmp(parts{2}, 'NTwInpSt') && num_stations == 0
            num_stations = str2double(parts{1});
        end
        if strcmp(parts{2}, 'TwrFADmp(1)')
            damp_pct = str2double(parts{1});
        end
        % Mode shape coefficients
        for k = 2:6
            kw = sprintf('TwFAM1Sh(%d)', k);
            if strcmp(parts{2}, kw)
                FA_mode1_sh(k-1) = str2double(parts{1});
            end
        end
    end

    % Detectar cabecera de tabla
    if ~in_table && contains(raw_line, 'HtFract') && contains(raw_line, 'TMassDen')
        in_table = true;
        fgetl(fid); % saltar unidades
        continue;
    end

    % Leer filas de la tabla
    if in_table && station_count < num_stations
        line2 = strtrim(raw_line);
        if isempty(line2) || line2(1) == '!', continue; end
        pts = strsplit(line2);
        pts = pts(~cellfun(@isempty, pts));
        if length(pts) < 3, continue; end
        v1 = str2double(pts{1});
        v2 = str2double(pts{2});
        v3 = str2double(pts{3});
        if isnan(v1) || isnan(v2) || isnan(v3), continue; end
        station_count = station_count + 1;
        HtFract_arr(station_count)  = v1;
        TMassDen_arr(station_count) = v2;
        TwFAStif_arr(station_count) = v3;
        if station_count == num_stations
            in_table = false;
        end
    end
end

fclose(fid);

if station_count == 0
    warning('parse_elastodyn_tower: No se encontraron datos en: %s', tower_path);
    return;
end

tw.HtFract  = HtFract_arr(:);
tw.TMassDen = TMassDen_arr(:);
tw.TwFAStif = TwFAStif_arr(:);

% Altura absoluta
if ~isnan(tw.height)
    h_abs = tw.HtFract * tw.height;
else
    h_abs = tw.HtFract;
    warning('parse_elastodyn_tower: TowerHt no encontrado');
end

% Masa total
tw.mass = trapz(h_abs, tw.TMassDen);

% --- Frecuencia natural por cociente de Rayleigh ---
% Mode shape: φ(x) = c2*x² + c3*x³ + c4*x⁴ + c5*x⁵ + c6*x⁶
% donde x = HtFract ∈ [0,1]
x = tw.HtFract;
c = FA_mode1_sh;  % [c2 c3 c4 c5 c6]

% Evaluar mode shape en cada estación
phi = c(1)*x.^2 + c(2)*x.^3 + c(3)*x.^4 + c(4)*x.^5 + c(5)*x.^6;

% Segunda derivada del mode shape (curvatura)
% φ'' = 2c2 + 6c3*x + 12c4*x² + 20c5*x³ + 30c6*x⁴
% dividido por L² para convertir de fracción a metros
if ~isnan(tw.height) && tw.height > 0
    L = tw.height;
    phi2 = (2*c(1) + 6*c(2)*x + 12*c(3)*x.^2 + 20*c(4)*x.^3 + 30*c(5)*x.^4) / L^2;
else
    phi2 = 2*c(1) + 6*c(2)*x + 12*c(3)*x.^2 + 20*c(4)*x.^3 + 30*c(5)*x.^4;
end

% Cociente de Rayleigh: ω² = ∫EI·φ''² dx / ∫m·φ² dx
numerator   = trapz(h_abs, tw.TwFAStif .* phi2.^2);
denominator = trapz(h_abs, tw.TMassDen .* phi.^2);

if denominator > 0
    omega_n = sqrt(numerator / denominator);  % [rad/s] frecuencia natural no amortiguada
else
    omega_n = 0.32;  % fallback
    warning('parse_elastodyn_tower: No se pudo calcular eigfreq, usando 0.32 rad/s');
end

% Amortiguamiento estructural
tw.damp    = damp_pct / 100;   % fracción crítica

% Frecuencia amortiguada (lo que usa el script NREL 5MW)
tw.eigfreq = omega_n / sqrt(1 - tw.damp^2);

fprintf('parse_elastodyn_tower: OK\n');
fprintf('  Tower mass    : %.1f kg\n',   tw.mass);
fprintf('  Tower eigfreq : %.4f rad/s (%.4f Hz)\n', tw.eigfreq, tw.eigfreq/(2*pi));
fprintf('  Tower damp    : %.2f %%\n',   damp_pct);

end