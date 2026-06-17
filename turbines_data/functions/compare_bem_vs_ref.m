function compare_bem_vs_ref(bem, ref)
% COMPARE_BEM_VS_REF  Compara visualmente las tablas Cp/Ct del BEM
%                     contra el fichero de referencia de ROSCO.

% --- Transponer BEM a (n_tsr x n_beta) para igualar convención de ref ---
cp_bem_T = bem.cp';
ct_bem_T = bem.ct';

% --- Interpolar BEM al grid de referencia ---
[BETA_ref, TSR_ref] = meshgrid(ref.beta, ref.tsr);
cp_bem_interp = interp2(bem.beta, bem.tsr, cp_bem_T, BETA_ref, TSR_ref, 'linear', NaN);
ct_bem_interp = interp2(bem.beta, bem.tsr, ct_bem_T, BETA_ref, TSR_ref, 'linear', NaN);

% --- Métricas de error ---
mask = ~isnan(cp_bem_interp) & ~isnan(ref.cp);
cp_err = cp_bem_interp(mask) - ref.cp(mask);
ct_err = ct_bem_interp(mask) - ref.ct(mask);

fprintf('\n--- Comparación BEM vs Referencia ---\n');
fprintf('Cp:  RMSE = %.4f,  MaxErr = %.4f,  Bias = %.4f\n', ...
    rms(cp_err), max(abs(cp_err)), mean(cp_err));
fprintf('Ct:  RMSE = %.4f,  MaxErr = %.4f,  Bias = %.4f\n', ...
    rms(ct_err), max(abs(ct_err)), mean(ct_err));
fprintf('Cp max BEM = %.4f,  Cp max REF = %.4f\n', max(bem.cp(:)), max(ref.cp(:)));

% --- Figura 1: Cp contour BEM vs REF ---
figure('Name','Cp: BEM vs Referencia');
subplot(1,3,1);
contourf(bem.beta, bem.tsr, cp_bem_T, 20);
colorbar; xlabel('Pitch (deg)'); ylabel('TSR (-)'); title('Cp BEM'); clim([0 0.5]);
subplot(1,3,2);
contourf(ref.beta, ref.tsr, ref.cp, 20);
colorbar; xlabel('Pitch (deg)'); ylabel('TSR (-)'); title('Cp Referencia'); clim([0 0.5]);
subplot(1,3,3);
contourf(ref.beta, ref.tsr, cp_bem_interp - ref.cp, 20);
colorbar; xlabel('Pitch (deg)'); ylabel('TSR (-)'); title('Cp BEM - Ref');

% --- Figura 2: Ct contour BEM vs REF ---
figure('Name','Ct: BEM vs Referencia');
subplot(1,3,1);
contourf(bem.beta, bem.tsr, ct_bem_T, 20);
colorbar; xlabel('Pitch (deg)'); ylabel('TSR (-)'); title('Ct BEM'); clim([0 1]);
subplot(1,3,2);
contourf(ref.beta, ref.tsr, ref.ct, 20);
colorbar; xlabel('Pitch (deg)'); ylabel('TSR (-)'); title('Ct Referencia'); clim([0 1]);
subplot(1,3,3);
contourf(ref.beta, ref.tsr, ct_bem_interp - ref.ct, 20);
colorbar; xlabel('Pitch (deg)'); ylabel('TSR (-)'); title('Ct BEM - Ref');

% --- Figura 3: Cortes Cp/Ct vs pitch a TSR fijo ---
figure('Name','Cortes Cp/Ct vs pitch a TSR fijo');
tsr_cuts = [6, 8, 9, 11];
cols = lines(length(tsr_cuts));

subplot(1,2,1); hold on;
h_bem = gobjects(length(tsr_cuts),1);
h_ref = gobjects(length(tsr_cuts),1);
for i = 1:length(tsr_cuts)
    t = tsr_cuts(i);
    [~,ib] = min(abs(bem.tsr - t));
    [~,ir] = min(abs(ref.tsr - t));
    h_bem(i) = plot(bem.beta, bem.cp(:,ib), '-',  'Color', cols(i,:), 'LineWidth', 2);
    h_ref(i) = plot(ref.beta, ref.cp(ir,:), '--', 'Color', cols(i,:), 'LineWidth', 1.5);
end
xlabel('Pitch (deg)'); ylabel('Cp (-)'); title('Cp vs pitch a TSR fijo'); grid on;
leg_str = arrayfun(@(t) sprintf('TSR=%.0f', t), tsr_cuts, 'UniformOutput', false);
legend(h_bem, leg_str, 'Location','northeast');
text(0.05, 0.05, '— BEM   - - Ref', 'Units','normalized');

subplot(1,2,2); hold on;
for i = 1:length(tsr_cuts)
    t = tsr_cuts(i);
    [~,ib] = min(abs(bem.tsr - t));
    [~,ir] = min(abs(ref.tsr - t));
    plot(bem.beta, bem.ct(:,ib), '-',  'Color', cols(i,:), 'LineWidth', 2);
    plot(ref.beta, ref.ct(ir,:), '--', 'Color', cols(i,:), 'LineWidth', 1.5);
end
xlabel('Pitch (deg)'); ylabel('Ct (-)'); title('Ct vs pitch a TSR fijo'); grid on;

% --- Figura 4: Cortes Cp/Ct vs TSR a pitch fijo ---
figure('Name','Cortes Cp/Ct vs TSR a pitch fijo');
beta_cuts = [0, 5, 10, 15];
cols = lines(length(beta_cuts));

subplot(1,2,1); hold on;
h_bem = gobjects(length(beta_cuts),1);
for i = 1:length(beta_cuts)
    b = beta_cuts(i);
    [~,ib] = min(abs(bem.beta - b));
    [~,ir] = min(abs(ref.beta - b));
    h_bem(i) = plot(bem.tsr, bem.cp(ib,:), '-',  'Color', cols(i,:), 'LineWidth', 2);
               plot(ref.tsr, ref.cp(:,ir), '--', 'Color', cols(i,:), 'LineWidth', 1.5);
end
xlabel('TSR (-)'); ylabel('Cp (-)'); title('Cp vs TSR a pitch fijo'); grid on;
leg_str = arrayfun(@(b) sprintf('\\beta=%d°', b), beta_cuts, 'UniformOutput', false);
legend(h_bem, leg_str, 'Location','northeast');
text(0.05, 0.05, '— BEM   - - Ref', 'Units','normalized');

subplot(1,2,2); hold on;
for i = 1:length(beta_cuts)
    b = beta_cuts(i);
    [~,ib] = min(abs(bem.beta - b));
    [~,ir] = min(abs(ref.beta - b));
    plot(bem.tsr, bem.ct(ib,:), '-',  'Color', cols(i,:), 'LineWidth', 2);
    plot(ref.tsr, ref.ct(:,ir), '--', 'Color', cols(i,:), 'LineWidth', 1.5);
end
xlabel('TSR (-)'); ylabel('Ct (-)'); title('Ct vs TSR a pitch fijo'); grid on;

end