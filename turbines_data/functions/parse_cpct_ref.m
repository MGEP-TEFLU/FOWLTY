function ref = parse_cpct_ref(ref_path)
% PARSE_CPCT_REF  Lee el fichero Cp_Ct_Cq generado por ROSCO toolbox.
%
% Formato del fichero:
%   # comentarios
%   # Pitch angle vector, N entries (deg)
%   val1  val2  ...
%   # TSR vector, M entries (-)
%   val1  val2  ...
%   # Power coefficient
%   (matriz M x N, una fila por TSR)
%   # Thrust coefficient
%   (matriz M x N)
%
% OUTPUT:
%   ref.cp   - matriz Cp (n_tsr x n_pitch)
%   ref.ct   - matriz Ct (n_tsr x n_pitch)
%   ref.beta - vector pitch [deg]  (1 x n_pitch)
%   ref.tsr  - vector TSR [-]      (1 x n_tsr)

fid = fopen(ref_path, 'r');
if fid == -1
    error('parse_cpct_ref: No se puede abrir: %s', ref_path);
end

ref.beta = [];
ref.tsr  = [];
ref.cp   = [];
ref.ct   = [];

state = 'find_pitch';  % estados: find_pitch, find_tsr, find_cp, find_ct, read_cp, read_ct

cp_rows = {};
ct_rows = {};

while ~feof(fid)
    raw = fgetl(fid);
    if ~ischar(raw), continue; end
    line = strtrim(raw);
    if isempty(line), continue; end

    % Detectar secciones por comentarios #
    if line(1) == '#'
        low = lower(line);
        if contains(low, 'pitch angle')
            state = 'read_pitch';
        elseif contains(low, 'tsr vector')
            state = 'read_tsr';
        elseif contains(low, 'wind speed')
            state = 'skip';  % línea de wind speed, ignorar
        elseif contains(low, 'power coefficient')
            state = 'read_cp';
        elseif contains(low, 'thrust coefficient')
            state = 'read_ct';
        elseif contains(low, 'torque coefficient')
            state = 'read_cq';  % ignoramos Cq
        end
        continue;
    end

    % Leer datos según estado
    vals = str2double(strsplit(line));
    vals = vals(~isnan(vals));
    if isempty(vals), continue; end

    switch state
        case 'read_pitch'
            ref.beta = vals;
            state = 'idle';
        case 'read_tsr'
            ref.tsr = vals;
            state = 'idle';
        case 'read_cp'
            cp_rows{end+1} = vals; %#ok<AGROW>
        case 'read_ct'
            ct_rows{end+1} = vals; %#ok<AGROW>
        case {'skip','idle','read_cq'}
            % nada
    end
end

fclose(fid);

% Montar matrices
if ~isempty(cp_rows)
    ref.cp = cell2mat(cp_rows');   % (n_tsr x n_pitch)
end
if ~isempty(ct_rows)
    ref.ct = cell2mat(ct_rows');
end

fprintf('parse_cpct_ref: OK\n');
fprintf('  Pitch: %.1f to %.1f deg (%d puntos)\n', ...
    ref.beta(1), ref.beta(end), length(ref.beta));
fprintf('  TSR  : %.1f to %.1f (%d puntos)\n', ...
    ref.tsr(1), ref.tsr(end), length(ref.tsr));
fprintf('  Cp max ref = %.4f\n', max(ref.cp(:)));
fprintf('  Ct max ref = %.4f\n', max(ref.ct(:)));

end