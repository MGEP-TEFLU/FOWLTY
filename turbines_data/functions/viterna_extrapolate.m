function af_out = viterna_extrapolate(af_in, AR)
% VITERNA_EXTRAPOLATE  Extiende una polar aerodinámica hasta ±180° usando
%                      el método de Viterna-Corrigan para ángulos de ataque
%                      en stall profundo.
%
% El método de Viterna asume que para alpha > alpha_stall:
%   Cl = A1*sin(2*alpha) + A2*cos(alpha)^2/sin(alpha)
%   Cd = B1*sin(alpha)^2 + B2*cos(alpha)
%
% donde A1, A2, B1, B2 se calculan a partir de los coeficientes en el
% punto de máximo Cl (stall) y el aspect ratio de la pala.
%
% INPUT:
%   af_in : struct con campos .alpha [deg], .cl [-], .cd [-]
%   AR    : aspect ratio de la pala (longitud/cuerda media) [-]
%           Para turbinas eólicas típicamente AR ≈ 17
%           Si no se conoce, usar AR = 17 (valor por defecto)
%
% OUTPUT:
%   af_out : struct con mismos campos pero extendido a [-180, 180] deg
%
% REFERENCIA:
%   Viterna & Corrigan (1981), "Fixed pitch rotor performance of large
%   horizontal axis wind turbines", DOE/NASA Workshop on Large HAWTs.
%
% EJEMPLO:
%   af_ext = viterna_extrapolate(af_raw, 17);

if nargin < 2 || isempty(AR)
    AR = 17;
end

alpha = af_in.alpha;  % [deg]
cl    = af_in.cl;
cd    = af_in.cd;

% --- Encontrar el punto de stall (máximo Cl en zona positiva) ---
% Solo miramos el rango -20° a 30° para evitar artefactos
mask_pos = alpha > -20 & alpha < 40;
if ~any(mask_pos)
    af_out = af_in;
    warning('viterna_extrapolate: no hay datos en rango -20 a 40 deg');
    return;
end

[cl_max, i_stall] = max(cl(mask_pos));
idx_pos = find(mask_pos);
i_stall = idx_pos(i_stall);

alpha_stall = alpha(i_stall);   % [deg]
cl_stall    = cl(i_stall);
cd_stall    = cd(i_stall);

% Limitar alpha_stall a máximo 30° (perfiles cilíndricos tienen stall raro)
alpha_stall = min(alpha_stall, 30);
alpha_stall = max(alpha_stall, 5);

% --- Coeficiente de drag máximo (a 90°) ---
% Cd_max = 1.11 + 0.018*AR  (Viterna-Corrigan)
cd_max = min(1.11 + 0.018*AR, 2.0);

% Asegurarse de que cd_stall no sea mayor que cd_max
cd_stall = min(cd_stall, cd_max * 0.9);

% --- Coeficientes de Viterna ---
as  = alpha_stall * pi/180;   % [rad]
B1  = cd_max;
B2  = (cd_stall - cd_max * sin(as)^2) / cos(as);
A1  = B1/2;
A2  = (cl_stall - cd_max * sin(as)*cos(as)) * sin(as) / cos(as)^2;

% --- Grid de salida: 180 puntos de -180 a 180 deg ---
alpha_ext = (-180:2:180)';   % [deg]
cl_ext    = zeros(size(alpha_ext));
cd_ext    = zeros(size(alpha_ext));

for k = 1:length(alpha_ext)
    a_deg = alpha_ext(k);
    a_rad = a_deg * pi/180;

    if a_deg >= -alpha_stall && a_deg <= alpha_stall
        % Zona lineal: interpolar directamente de los datos originales
        cl_ext(k) = interp1(alpha, cl, a_deg, 'linear', 'extrap');
        cd_ext(k) = interp1(alpha, cd, a_deg, 'linear', 'extrap');

    elseif a_deg > alpha_stall && a_deg <= 180 - alpha_stall
        % Zona de stall profundo positiva: Viterna
        cl_ext(k) = A1*sin(2*a_rad) + A2*cos(a_rad)^2/max(sin(a_rad), 1e-6);
        cd_ext(k) = B1*sin(a_rad)^2 + B2*cos(a_rad);

    elseif a_deg < -alpha_stall && a_deg >= -(180 - alpha_stall)
        % Zona de stall profundo negativa: Viterna con simetría
        cl_ext(k) = -(A1*sin(-2*a_rad) + A2*cos(a_rad)^2/max(sin(-a_rad), 1e-6));
        cd_ext(k) = B1*sin(a_rad)^2 - B2*cos(a_rad);

    else
        % Zona de flujo invertido (cerca de ±180°)
        % Usar simetría: comportamiento similar a alpha cercano a 0° pero invertido
        if a_deg > 0
            a_ref = 180 - a_deg;  % [deg]
        else
            a_ref = -(180 + a_deg);
        end
        a_ref_rad = a_ref * pi/180;
        cl_ext(k) = A1*sin(2*a_ref_rad) + A2*cos(a_ref_rad)^2/max(sin(a_ref_rad), 1e-6);
        cl_ext(k) = -cl_ext(k);
        cd_ext(k) = B1*sin(a_ref_rad)^2 + B2*cos(a_ref_rad);
    end
end

% Cd mínimo físico
cd_ext = max(cd_ext, 0.001);

af_out.alpha       = alpha_ext;
af_out.cl          = cl_ext;
af_out.cd          = cd_ext;
af_out.path        = af_in.path;
af_out.alpha_stall = alpha_stall;
af_out.cl_max      = cl_max;
af_out.cd_max      = cd_max;

end