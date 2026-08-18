function results = study_convergence(opts)
   % runs a few convergence studies and plots them
   arguments
       opts.Airfoil = '2412'
       opts.Nvalues = [50 100 200 400 800]        
       opts.AlphaTest = 5                         
       opts.Uinf = 10                               
       opts.WakeLengths = logspace(-2, 1, 15)       
       opts.AlphaRange = -5:1:15                    
       opts.Save = false
       opts.OutDir = 'figures'
       opts.XfoilFile = 'data/Xfoil.txt'
   end

   if opts.Save && ~isfolder(opts.OutDir); mkdir(opts.OutDir); end

   fprintf('\n=== Convergence studies: NACA %s ===\n', opts.Airfoil);

   %% study 0
   fprintf('\n[0/4] Effect of the gauge condition at alpha = %.1f deg\n', opts.AlphaTest);

   nN0 = numel(opts.Nvalues);
   Cl_nogauge = zeros(nN0, 1);
   Cl_gauge   = zeros(nN0, 1);
   cond_ng    = zeros(nN0, 1);

   wState = warning('off', 'panel_strengths:noGauge');
   for i = 1:nN0
       sFree = solve_case(opts.Airfoil, opts.Nvalues(i), opts.AlphaTest, opts.Uinf, 'FixGauge', false);
       sFix  = solve_case(opts.Airfoil, opts.Nvalues(i), opts.AlphaTest, opts.Uinf, 'FixGauge', true);
       Cl_nogauge(i) = sFree.coef.Cl;
       Cl_gauge(i)   = sFix.coef.Cl;
       cond_ng(i)    = sFree.info.cond;
       fprintf('   N = %4d   Cl no gauge = %.5f   Cl gauge = %.5f   cond = %.2e\n', ...
           opts.Nvalues(i), Cl_nogauge(i), Cl_gauge(i), cond_ng(i));
   end
   warning(wState);

   fig0 = figure('Name', 'Effect of the gauge condition', 'Color', 'w');
   plot(opts.Nvalues, Cl_nogauge, 'o--', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'DisplayName', 'Original (singular)');
   hold on
   plot(opts.Nvalues, Cl_gauge, 's-', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'DisplayName', 'Gauge fixed');
   xlabel('Number of panels, N'); ylabel('C_L');
   title('Removing the constant-doublet null mode');
   legend('Location', 'best'); grid on; box on; hold off

   %% study 1 - panel count
   fprintf('\n[1/4] Panel count convergence at alpha = %.1f deg\n', opts.AlphaTest);

   nN = numel(opts.Nvalues);
   Cl_N = zeros(nN, 1); Cd_N = zeros(nN, 1); Cm_N = zeros(nN, 1);
   time_N = zeros(nN, 1); cond_N = zeros(nN, 1);

   for i = 1:nN
       sol = solve_case(opts.Airfoil, opts.Nvalues(i), opts.AlphaTest, opts.Uinf);
       Cl_N(i) = sol.coef.Cl;
       Cd_N(i) = sol.coef.Cd_press;
       Cm_N(i) = sol.coef.Cm_quarter;
       time_N(i) = sol.runtime;
       cond_N(i) = sol.info.cond;
       fprintf('   N = %4d   Cl = %8.5f   Cd(err) = %9.2e   %6.3f s\n', opts.Nvalues(i), Cl_N(i), Cd_N(i), time_N(i));
   end

   Cl_ref = Cl_N(end);
   err = abs(Cl_N(1:end-1) - Cl_ref);
   Nfit = opts.Nvalues(1:end-1).';
   valid = err > 0;
   if nnz(valid) >= 2
       pFit = polyfit(log(Nfit(valid)), log(err(valid)), 1);
       order = -pFit(1);
   else
       order = NaN;
   end
   fprintf('   Observed order of convergence: %.2f\n', order);

   fig1 = figure('Name', 'Panel count convergence', 'Color', 'w');
   tiledlayout(1, 2, 'TileSpacing', 'compact');
   nexttile;
   plot(opts.Nvalues, Cl_N, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
   xlabel('Number of panels, N'); ylabel('C_L');
   title('C_L convergence'); grid on; box on
   nexttile;
   loglog(Nfit(valid), err(valid), 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
   xlabel('Number of panels, N'); ylabel('|C_L - C_L^{ref}|');
   title(sprintf('Error (order \\approx %.2f)', order)); grid on; box on

   %% study 2 - wake length
   fprintf('\n[2/4] Wake length sensitivity at alpha = %.1f deg\n', opts.AlphaTest);

   nW = numel(opts.WakeLengths);
   Cl_W = zeros(nW, 1);
   Nw = 200;
   for i = 1:nW
       sol = solve_case(opts.Airfoil, Nw, opts.AlphaTest, opts.Uinf, 'WakeLength', opts.WakeLengths(i));
       Cl_W(i) = sol.coef.Cl;
   end
   Cl_W_ref = Cl_W(end);

   fig2 = figure('Name', 'Wake length sensitivity', 'Color', 'w');
   semilogx(opts.WakeLengths, Cl_W, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
   hold on
   yline(Cl_W_ref, '--', 'Converged value', 'Color', [0.4 0.4 0.4]);
   xline(0.9, ':', 'Original hardcoded 0.9c', 'Color', [0.85 0.15 0.15]);
   xlabel('Wake panel length [chords]'); ylabel('C_L');
   title('Sensitivity to wake length'); grid on; box on

   %% study 3 - lift curve
   fprintf('\n[3/4] Lift curve vs thin-aerofoil theory and XFOIL\n');

   Nsweep = [100, 200, 500];
   Cl_vals = zeros(numel(Nsweep), numel(opts.AlphaRange));
   Cm_vals = zeros(numel(Nsweep), numel(opts.AlphaRange));
   Cn_vals = zeros(numel(Nsweep), numel(opts.AlphaRange));
   momentRef = 0.25; % has to match solve_case's default

   for i = 1:numel(Nsweep)
       for j = 1:numel(opts.AlphaRange)
           sol = solve_case(opts.Airfoil, Nsweep(i), opts.AlphaRange(j), opts.Uinf);
           Cl_vals(i,j) = sol.coef.Cl;
           Cm_vals(i,j) = sol.coef.Cm_quarter;
           Cn_vals(i,j) = sol.coef.Cn;
       end
   end

   slopeSolver = polyfit(opts.AlphaRange(:)*pi/180, Cl_vals(end,:).', 1);
   alpha_L0 = -slopeSolver(2) / slopeSolver(1) * 180/pi;
   Cl_thin = 2*pi * (opts.AlphaRange*pi/180 - alpha_L0*pi/180);

   fprintf('   Solver slope : %.4f per rad\n', slopeSolver(1));
   fprintf('   Thin theory  : %.4f per rad\n', 2*pi);

   fig3 = figure('Name', 'Lift curve', 'Color', 'w');
   hold on
   cols = lines(numel(Nsweep));
   for i = 1:numel(Nsweep)
       plot(opts.AlphaRange, Cl_vals(i,:), '-', 'LineWidth', 1.4, 'Color', cols(i,:), ...
           'DisplayName', sprintf('Panel method, N = %d', Nsweep(i)));
   end
   plot(opts.AlphaRange, Cl_thin, 'k--', 'LineWidth', 1.4, 'DisplayName', 'Thin-aerofoil theory (2\pi)');

   haveXfoil = false;
   try
       polar = load_xfoil(opts.XfoilFile);
       inRange = polar.alpha >= min(opts.AlphaRange) & polar.alpha <= max(opts.AlphaRange);
       plot(polar.alpha(inRange), polar.Cl(inRange), 'k-', 'LineWidth', 1.6, 'DisplayName', 'XFOIL (viscous)');
       haveXfoil = true;
   catch ME
       fprintf('   XFOIL comparison skipped: %s\n', ME.message);
   end

   xlabel('\alpha [\circ]', 'Interpreter', 'tex');
   ylabel('C_L', 'Interpreter', 'tex');
   title(sprintf('Lift curve, NACA %s', opts.Airfoil));
   legend('Location', 'northwest');
   grid on; box on; hold off

   %% study 4 - aero centre and centre of pressure
   fprintf('\n[4/4] Aerodynamic centre and centre of pressure\n');

   Cm_ac = Cm_vals(end,:);
   Cn_ac = Cn_vals(end,:);

   pAC = polyfit(Cn_ac, Cm_ac, 1);
   x_ac = momentRef - pAC(1);

   x_cp = momentRef - Cm_ac ./ Cn_ac;

   fprintf('   Aerodynamic centre: x/c = %.4f\n', x_ac);
   fprintf('   Centre of pressure: %.4f at alpha=%.1f -> %.4f at alpha=%.1f\n', ...
       x_cp(1), opts.AlphaRange(1), x_cp(end), opts.AlphaRange(end));

   fig4 = figure('Name', 'Aerodynamic centre and centre of pressure', 'Color', 'w');
   tiledlayout(1, 2, 'TileSpacing', 'compact');

   nexttile;
   plot(opts.AlphaRange, x_cp, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', ...
       'Color', [0.10 0.35 0.75], 'DisplayName', 'Centre of pressure');
   hold on
   yline(x_ac, '--', sprintf('AC, x/c = %.3f', x_ac), 'Color', [0.85 0.15 0.15], 'LineWidth', 1.3);
   xlabel('\alpha [\circ]', 'Interpreter', 'tex'); ylabel('x/c');
   title('CP vs AC'); legend('Location', 'best'); grid on; box on; hold off

   nexttile;
   [xShape, zShape] = panelgen(opts.Airfoil, 200, 0);
   plot(xShape(1:end-1), zShape(1:end-1), 'k', 'LineWidth', 1.3);
   hold on
   plot(x_ac, 0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.85 0.15 0.15], 'MarkerEdgeColor', 'k');
   yline(0, ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
   text(x_ac, 0.06, sprintf('AC, x/c = %.3f', x_ac), 'HorizontalAlignment', 'center', 'FontSize', 9);
   axis equal; xlim([-0.05, 1.05]);
   xlabel('x/c'); ylabel('z/c');
   title(sprintf('NACA %s section', opts.Airfoil)); box on; hold off

   %% results
   results = struct( ...
       'Cl_nogauge', Cl_nogauge, 'Cl_gauge', Cl_gauge, 'cond_nogauge', cond_ng, ...
       'Nvalues', opts.Nvalues, 'Cl_N', Cl_N, 'Cd_error_N', Cd_N, 'Cm_N', Cm_N, ...
       'cond_N', cond_N, 'runtime_N', time_N, 'order', order, ...
       'wakeLengths', opts.WakeLengths, 'Cl_wake', Cl_W, ...
       'alphaRange', opts.AlphaRange, 'Nsweep', Nsweep, 'Cl_alpha', Cl_vals, ...
       'Cm_alpha', Cm_vals, 'Cn_alpha', Cn_vals, ...
       'liftSlope', slopeSolver(1), 'alpha_L0', alpha_L0, 'haveXfoil', haveXfoil, ...
       'x_ac', x_ac, 'x_cp', x_cp);

   if opts.Save
       exportgraphics(fig0, fullfile(opts.OutDir, 'gauge_effect.png'), 'Resolution', 200);
       exportgraphics(fig1, fullfile(opts.OutDir, 'convergence_panels.png'), 'Resolution', 200);
       exportgraphics(fig2, fullfile(opts.OutDir, 'convergence_wake.png'), 'Resolution', 200);
       exportgraphics(fig3, fullfile(opts.OutDir, 'lift_curve.png'), 'Resolution', 200);
       exportgraphics(fig4, fullfile(opts.OutDir, 'ac_and_cp.png'), 'Resolution', 200);
       fprintf('\n   Figures written to %s\n', opts.OutDir);
   end

   fprintf('\n=== Studies complete ===\n');
end

function polar = load_xfoil(filename)
   % reads an XFOIL polar file
   arguments
       filename (1,:) char
   end

   if ~isfile(filename)
       error('load_xfoil:notFound', 'could not find %s', filename);
   end

   raw = readmatrix(filename, 'FileType', 'text');
   if isempty(raw) || size(raw, 2) < 7
       error('load_xfoil:badFormat', '%s doesnt look like an XFOIL polar', filename);
   end
   raw = raw(:, 1:7);

   valid = all(isfinite(raw), 2);
   valid = valid & abs(raw(:,1)) <= 90;
   raw = raw(valid, :);

   if isempty(raw)
       error('load_xfoil:noData', 'no usable rows found in %s', filename);
   end

   [~, order] = sort(raw(:,1));
   raw = raw(order, :);
   [~, unique_idx] = unique(raw(:,1), 'stable');
   raw = raw(unique_idx, :);

   polar.alpha  = raw(:,1);
   polar.Cl     = raw(:,2);
   polar.Cd     = raw(:,3);
   polar.Cdp    = raw(:,4);
   polar.Cm     = raw(:,5);
   polar.topXtr = raw(:,6);
   polar.botXtr = raw(:,7);
   polar.file   = filename;
end
