function bem = run_bem(blade, airfoils, ed, opts)
% RUN_BEM  Calcula las tablas Cp(beta,tsr) y Ct(beta,tsr) de una turbina
%          mediante Blade Element Momentum theory.
%
% INPUT:
%   blade    : struct de parse_aerodyn_blade (span, chord, twist, af_id)
%   airfoils : cell array de structs de parse_airfoil {af1, af2, ...}
%   ed       : struct de parse_elastodyn (TipRad, HubRad, NumBl)
%   opts     : (opcional) struct con opciones:
%                .beta  - vector de pitch [deg]    (default: -5:1:30)
%                .tsr   - vector de TSR [-]        (default: 2:0.25:14.5)
%                .rho   - densidad aire [kg/m³]    (default: 1.225)
%                .tol   - tolerancia conv. [-]     (default: 1e-6)
%                .maxit - iteraciones máximas      (default: 500)
%
% OUTPUT:
%   bem      : struct con:
%                .cp    - matriz Cp (n_beta x n_tsr)
%                .ct    - matriz Ct (n_beta x n_tsr)
%                .beta  - vector pitch [deg]
%                .tsr   - vector TSR [-]
%
% El BEM ignora coning, sweep y curvatura de pala (pala recta).
% Incluye: pérdidas de punta y raíz de Prandtl, corrección de Glauert
% (Buhl) para a > 0.4, e inducción tangencial.
%
% EJEMPLO:
%   airfoils = cell(1, ad.num_airfoils);
%   for i = 1:ad.num_airfoils
%       airfoils{i} = parse_airfoil(ad.airfoil_files{i});
%   end
%   bem = run_bem(blade, airfoils, ed);

% ----------------- Opciones por defecto -----------------
if nargin < 4, opts = struct(); end
if ~isfield(opts, 'beta'),  opts.beta  = -5:1:30;     end
if ~isfield(opts, 'tsr'),   opts.tsr   = 2:0.25:14.5; end
if ~isfield(opts, 'rho'),   opts.rho   = 1.225;       end
if ~isfield(opts, 'tol'),   opts.tol   = 1e-6;        end
if ~isfield(opts, 'maxit'), opts.maxit = 500;         end

beta_v = opts.beta(:)';   % pitch [deg], fila
tsr_v  = opts.tsr(:)';    % TSR, fila

R    = ed.TipRad;         % radio total [m]
Rhub = ed.HubRad;         % radio buje [m]
B    = double(ed.NumBl);  % número de palas

% ----------------- Geometría de pala -----------------
% La posición radial de los nodos AeroDyn es relativa a la raíz de pala,
% hay que sumar el radio del buje para obtener r absoluto
r     = blade.span + Rhub;     % [m] posición radial absoluta
chord = blade.chord;           % [m]
twist = blade.twist;           % [deg]
af_id = blade.af_id;           % [-]
n_nodes = length(r);

% Evitar singularidades en los extremos exactos
r(1)   = max(r(1),   Rhub + 1e-3);
r(end) = min(r(end), R    - 1e-3);

% Anchos de los elementos (integración por trapecios necesita dr)
% Usamos diferencias centradas
dr = zeros(n_nodes,1);
dr(1)       = (r(2)-r(1))/2 + (r(1)-Rhub);
dr(end)     = (r(end)-r(end-1))/2 + (R-r(end));
dr(2:end-1) = (r(3:end)-r(1:end-2))/2;

% ----------------- Preparar interpoladores de polares -----------------
% griddedInterpolant es mucho más rápido que interp1 en bucles
n_af = length(airfoils);
af_interp_cl = cell(n_af,1);
af_interp_cd = cell(n_af,1);
for i = 1:n_af
    [alpha_u, iu] = unique(airfoils{i}.alpha); % por si hay duplicados
    af_interp_cl{i} = griddedInterpolant(alpha_u, airfoils{i}.cl(iu), 'linear', 'nearest');
    af_interp_cd{i} = griddedInterpolant(alpha_u, airfoils{i}.cd(iu), 'linear', 'nearest');
end

% ----------------- Bucle principal -----------------
n_beta = length(beta_v);
n_tsr  = length(tsr_v);
cp = zeros(n_beta, n_tsr);
ct = zeros(n_beta, n_tsr);

fprintf('run_bem: calculando %d x %d = %d combinaciones...\n', ...
    n_beta, n_tsr, n_beta*n_tsr);
t_start = tic;

% Para Cp/Ct no necesitamos la velocidad de viento real:
% trabajamos en variables adimensionales con U=1 y omega=tsr*U/R
U = 1; % [m/s] velocidad de referencia (adimensional)

for i_tsr = 1:n_tsr
    tsr   = tsr_v(i_tsr);
    omega = tsr * U / R;    % [rad/s]

    for i_beta = 1:n_beta
        beta = beta_v(i_beta);  % [deg]

        % Acumuladores de potencia y empuje
        P_total = 0;
        T_total = 0;

        for k = 1:n_nodes
            rk = r(k);
            ck = chord(k);
            % Ángulo total de la sección: twist + pitch
            theta = (twist(k) + beta) * pi/180; % [rad]

            % Solidez local
            sigma = B * ck / (2*pi*rk);

            % Velocidad tangencial local
            lambda_r = omega * rk / U;

            % ---- Iteración BEM para este anillo ----
            a  = 0.3;   % inducción axial inicial
            ap = 0.0;   % inducción tangencial inicial
            
            converged = false;
            for it = 1:opts.maxit
                % Ángulo de flujo
                phi = atan2(U*(1-a), omega*rk*(1+ap)); % [rad]
                if phi <= 0
                    % Flujo invertido: anillo sin contribución útil
                    break;
                end
                sphi = sin(phi);
                cphi = cos(phi);

                % Pérdidas de Prandtl (tip + hub)
                f_tip = B/2 * (R - rk) / (rk * sphi);
                F_tip = 2/pi * acos(min(1, exp(-abs(f_tip))));
                f_hub = B/2 * (rk - Rhub) / (Rhub * sphi);
                F_hub = 2/pi * acos(min(1, exp(-abs(f_hub))));
                F = max(F_tip * F_hub, 1e-4);

                % Ángulo de ataque y coeficientes
                alpha_deg = (phi*180/pi) - (theta*180/pi);
                id = af_id(k);
                cl = af_interp_cl{id}(alpha_deg);
                cd = af_interp_cd{id}(alpha_deg);

                % Coeficientes normal y tangencial
                cn = cl*cphi + cd*sphi;
                ctan = cl*sphi - cd*cphi;

                % ---- Nuevo factor de inducción axial ----
                kappa = sigma * cn / (4 * F * sphi^2);
                if kappa <= 2/3
                    % Momentum theory estándar
                    a_new = kappa / (1 + kappa);
                else
                    % Corrección de Glauert empírica (Buhl)
                    % CT = 8/9 + (4F - 40/9)a + (50/9 - 4F)a²
                    g1 = 2*F*kappa - (10/9 - F);
                    g2 = 2*F*kappa - F*(4/3 - F);
                    g3 = 2*F*kappa - (25/9 - 2*F);
                    if abs(g3) < 1e-6
                        a_new = 1 - 1/(2*sqrt(g2));
                    else
                        a_new = (g1 - sqrt(max(g2,0))) / g3;
                    end
                end
                a_new = max(min(a_new, 0.95), -0.5); % limitar

                % ---- Nuevo factor de inducción tangencial ----
                kappap = sigma * ctan / (4 * F * sphi * cphi);
                ap_new = kappap / (1 - kappap);
                ap_new = max(min(ap_new, 0.5), -0.5);

                % Convergencia con relajación
                relax = 0.5;
                da  = a_new  - a;
                dap = ap_new - ap;
                a  = a  + relax*da;
                ap = ap + relax*dap;

                if abs(da) < opts.tol && abs(dap) < opts.tol
                    converged = true;
                    break;
                end
            end

            if ~converged && it == opts.maxit
                % no convergió: usamos el último valor (suele ser aceptable)
            end

            % ---- Fuerzas del elemento ----
            phi = atan2(U*(1-a), omega*rk*(1+ap));
            if phi <= 0, continue; end
            sphi = sin(phi); cphi = cos(phi);

            alpha_deg = (phi - theta)*180/pi;
            id = af_id(k);
            cl = af_interp_cl{id}(alpha_deg);
            cd = af_interp_cd{id}(alpha_deg);
            cn = cl*cphi + cd*sphi;
            ctan = cl*sphi - cd*cphi;

            % Velocidad relativa
            W2 = (U*(1-a))^2 + (omega*rk*(1+ap))^2;

            % Fuerza normal (empuje) y tangencial por unidad de longitud
            qdyn = 0.5 * opts.rho * W2 * ck;
            dT = B * qdyn * cn   * dr(k);          % [N]
            dQ = B * qdyn * ctan * rk * dr(k);     % [N·m]

            T_total = T_total + dT;
            P_total = P_total + dQ * omega;
        end

        % Coeficientes adimensionales
        A = pi * R^2;
        cp(i_beta, i_tsr) = P_total / (0.5 * opts.rho * A * U^3);
        ct(i_beta, i_tsr) = T_total / (0.5 * opts.rho * A * U^2);
    end

    if mod(i_tsr, 5) == 0
        fprintf('  TSR %.2f (%d/%d) - %.1f s\n', tsr, i_tsr, n_tsr, toc(t_start));
    end
end

fprintf('run_bem: completado en %.1f s\n', toc(t_start));

% ----------------- Salida -----------------
bem.cp   = cp;
bem.ct   = ct;
bem.beta = beta_v;
bem.tsr  = tsr_v;

fprintf('  Cp max = %.4f (beta=%.1f deg, tsr=%.2f)\n', max(cp(:)), ...
    beta_v(find(any(cp==max(cp(:)),2),1)), ...
    tsr_v(find(any(cp==max(cp(:)),1),1)));

end