function plat = scale_platform(deepcwind_path, R_new, DoF, ed)
% SCALE_PLATFORM  Escala los coeficientes hidrodinámicos de la plataforma
%                 DeepCWind (diseñada para NREL 5MW, R=63m) a una turbina
%                 de radio R_new.
%
% El escalado es geométrico basado en el ratio de radios de rotor:
%   scale_R  = R_new / R_NREL5MW   (escala lineal)
%   Fuerzas y rigideces horizontales ~ R²
%   Masas, inercias, volúmenes      ~ R³
%
% INPUT:
%   deepcwind_path : path al fichero DeepCWind_Original.mat
%   R_new          : radio del rotor de la nueva turbina [m]
%   DoF            : grados de libertad a considerar (e.g. [1 3] = surge+pitch)
%   ed             : struct de parse_elastodyn (para TowerHt, Twr2Shft, etc.)
%
% OUTPUT:
%   plat : struct compatible con FOWLTY con campos:
%            .DoF, .nDoF, .Fe, .Mu, .Kh, .Mass, .K_moo, .B_drag
%            .CGz, .d_N2P, .nF2pT, .pM2nM_surge, .pM2nM_pitch
%            .HydSys (state-space del sistema hidrodinámico)
%
% EJEMPLO:
%   plat = scale_platform('DeepCWind_Original.mat', ed.TipRad, [1 3], ed)

% --- Parámetros de referencia (NREL 5MW / DeepCWind) ---
R_ref   = 63.0;      % [m] radio rotor NREL 5MW
CGz_ref = -13.46;    % [m] centro de gravedad plataforma (desde SWL)
TowerHt_ref = 87.6;  % [m] altura torre NREL 5MW

% --- Factor de escala ---
scale_R  = R_new / R_ref;
scale_R2 = scale_R^2;
scale_R3 = scale_R^3;

fprintf('scale_platform: scaling DeepCWind by factor %.3f (R=%.1fm -> R=%.1fm)\n', ...
    scale_R, R_ref, R_new);

% --- Cargar datos DeepCWind ---
raw = load(deepcwind_path);

% Verificar campos necesarios
required = {'Fe', 'Mu', 'Kh', 'Mass', 'w'};
for i = 1:length(required)
    if ~isfield(raw, required{i})
        error('scale_platform: campo "%s" no encontrado en %s', required{i}, deepcwind_path);
    end
end

% --- Escalar coeficientes ---
% Fe: fuerza de excitación [N o N·m] en frecuencia → ~ R²
Fe_scaled = raw.Fe * scale_R2;

% Mu: masa añadida infinita [kg o kg·m²] → ~ R³
Mu_scaled = raw.Mu * scale_R3;

% Kh: rigidez hidrostática [N/m o N·m/rad] → ~ R³
Kh_scaled = raw.Kh * scale_R3;

% Mass: masa e inercia del sistema [kg o kg·m²] → ~ R³
Mass_scaled = raw.Mass * scale_R3;

% K_moo: rigidez de amarre → ~ R² (fuerza horizontal)
% Si existe en el fichero, escalarlo; si no, usar referencia escalada
if isfield(raw, 'K_moo')
    K_moo_scaled = raw.K_moo * scale_R2;
else
    % Valores de referencia DeepCWind para NREL 5MW
    K_moo_ref = [7.08e4   0      -1.08e5; ...
                 0        1.91e4  0     ; ...
                -1.07e5   0       8.73e7];
    K_moo_scaled = K_moo_ref * scale_R2;
    fprintf('  NOTE: K_moo not found in file, using scaled DeepCWind reference values.\n');
end

% B_drag: arrastre viscoso → ~ R²
if isfield(raw, 'B_drag')
    B_drag_scaled = raw.B_drag * scale_R2;
else
    B_drag_ref = diag([3.95e5, 3.88e6, 3.7e10]);
    B_drag_scaled = B_drag_ref * scale_R2;
    fprintf('  NOTE: B_drag not found in file, using scaled DeepCWind reference values.\n');
end

% --- Seleccionar DoFs ---
plat.DoF  = DoF;
plat.nDoF = length(DoF);

plat.Fe   = Fe_scaled(:, DoF);
plat.Mu   = Mu_scaled(DoF, DoF);
plat.Kh   = Kh_scaled(DoF, DoF);
plat.Mass = Mass_scaled(DoF, DoF);
plat.K_moo  = K_moo_scaled(DoF, DoF);
plat.B_drag = B_drag_scaled(DoF, DoF);

% --- Parámetros geométricos escalados ---
% CGz: distancia desde SWL al centro de gravedad → escala lineal
plat.CGz = CGz_ref * scale_R;

% Distancia desde la góndola al centro de rotación de la plataforma
hub_height_MSL = ed.TowerHt + ed.Twr2Shft;
plat.d_N2P = hub_height_MSL - plat.CGz;

% Mappings de movimiento plataforma → góndola
plat.nF2pT       = [zeros(plat.nDoF-1, 1); plat.d_N2P];
plat.pM2nM_surge = [0 1 zeros(1, 2*(plat.nDoF-1))];
plat.pM2nM_pitch = [zeros(1, 2*(plat.nDoF-1)) 0 plat.d_N2P];

% --- State-space hidrodinámico ---
% Construir SS con matrices escaladas
MassMu = pinv(plat.Mass + plat.Mu);

Ass_t = kron(eye(plat.nDoF), [0 1; 0 0]) + ...
        kron(MassMu*(plat.K_moo + plat.Kh), [0 0; -1 0]) + ...
        kron(MassMu*plat.B_drag, [0 0; 0 -1]);
Bss_t = kron(MassMu, [0; 1]);
CssV_t = kron(eye(plat.nDoF), [0 1]);

% Término de radiación: escalar A y B si existen
if isfield(raw, 'A') && isfield(raw, 'B') && isfield(raw, 'T_rad')
    % Aproximación SS de la fuerza de radiación
    % Usar los DoF seleccionados
    A_sc = raw.A(DoF, DoF, :) * scale_R3;
    B_sc = raw.B(DoF, DoF, :) * scale_R3;
    w_rad = raw.w;

    % Ajuste de SS desde función de transferencia (simplificado)
    % Usamos los valores del fichero original si existe K_rad y T_rad
    if isfield(raw, 'K_rad') && isfield(raw, 'T_rad')
        % Cargar SS de radiación pre-calculado y escalar
        % K_rad es la función impulso de radiación
        % Aproximar con SS de orden reducido usando el método de Prony
        try
            % Intentar cargar sysRad si existe en el mismo directorio
            rad_file = fullfile(fileparts(deepcwind_path), 'DeepCWind_13DoF_RadSS.mat');
            if exist(rad_file, 'file')
                rad = load(rad_file);
                sysRad = rad.sysRad;
                % Escalar matrices del SS de radiación
                sysRad.b = sysRad.b * scale_R3;
                sysRad.c = sysRad.c * scale_R3;
                nR = length(sysRad.a);
                Ass = [Ass_t Bss_t*sysRad.c; sysRad.b*CssV_t sysRad.a];
                Bss = [Bss_t; zeros(nR, plat.nDoF)];
                CssP  = [kron(eye(plat.nDoF),[1 0]) zeros(plat.nDoF,nR)];
                CssV  = [CssV_t zeros(plat.nDoF,nR)];
                CssPV = [eye(2*plat.nDoF) zeros(2*plat.nDoF,nR)];
                plat.HydSys = ss(Ass, Bss, CssPV, 0);
                fprintf('  Radiation SS loaded and scaled from: %s\n', rad_file);
            else
                warning('scale_platform: DeepCWind_13DoF_RadSS.mat not found, using simplified SS (no radiation term).');
                plat.HydSys = build_simple_ss(Ass_t, Bss_t, CssV_t, plat.nDoF);
            end
        catch ME
            warning('scale_platform: Error loading radiation SS: %s', ME.message);
            plat.HydSys = build_simple_ss(Ass_t, Bss_t, CssV_t, plat.nDoF);
        end
    else
        plat.HydSys = build_simple_ss(Ass_t, Bss_t, CssV_t, plat.nDoF);
    end
else
    plat.HydSys = build_simple_ss(Ass_t, Bss_t, CssV_t, plat.nDoF);
end

fprintf('scale_platform: OK\n');
fprintf('  Scale factor  : %.3f\n',  scale_R);
fprintf('  CGz           : %.2f m\n', plat.CGz);
fprintf('  d_N2P         : %.2f m\n', plat.d_N2P);
fprintf('  HydSys order  : %d\n',     order(plat.HydSys));

end

% -------------------------------------------------------------------------
function HydSys = build_simple_ss(Ass_t, Bss_t, CssV_t, nDoF)
% SS simplificado sin término de radiación
CssPV = eye(2*nDoF);
HydSys = ss(Ass_t, Bss_t, CssPV, 0);
end