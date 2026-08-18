function figs = plot_results(sol, opts)
   arguments
       sol
       opts.Save = false
       opts.OutDir = 'figures'
       opts.Xlim = [-0.5, 1.5]      
       opts.Zlim = [-0.75, 0.75]    
       opts.GridSize = 150         
       opts.Skip = 6                
   end

   if opts.Save && ~isfolder(opts.OutDir)
       mkdir(opts.OutDir);
   end

   figsField = plot_flowfield(sol, opts);
   figPressure = plot_pressure(sol, opts);
   figs = [figsField, figPressure];
end

function figs = plot_flowfield(sol, opts)
   N = sol.N;

   xg = linspace(opts.Xlim(1), opts.Xlim(2), opts.GridSize);
   zg = linspace(opts.Zlim(1), opts.Zlim(2), opts.GridSize);
   [xp, zp] = meshgrid(xg, zg);

   [u, w, Cp] = velocity(xp, zp, sol.mu, sol.x, sol.z, sol.u_inf, sol.alpha, N);

   label = sprintf(['NACA-%s\n\\alpha = %.1f\\circ\nU_\\infty = %.1f m s^{-1}\n' ...
                    'N = %d\nC_L = %.4f'], ...
                    sol.airfoil, sol.alpha, sol.u_inf, N, sol.coef.Cl);

   figs = gobjects(1,3);

   % arrow plot
   figs(1) = figure('Name', 'Velocity vectors', 'Color', 'w');
   k = opts.Skip;
   quiver(xp(1:k:end, 1:k:end), zp(1:k:end, 1:k:end), ...
          u(1:k:end, 1:k:end),  w(1:k:end, 1:k:end), 'Color', [0.85 0.15 0.15]);
   hold on
   fill(sol.x(1:N+1), sol.z(1:N+1), [0.15 0.15 0.15], 'EdgeColor', 'k', 'LineWidth', 1.2);
   finish_axes(sprintf('Velocity vectors around NACA-%s', sol.airfoil), opts, label);

   % streamlines
   figs(2) = figure('Name', 'Streamlines', 'Color', 'w');
   h = streamslice(xp, zp, u, w, 2);
   set(h, 'Color', [0.10 0.35 0.75]);
   hold on
   fill(sol.x(1:N+1), sol.z(1:N+1), [0.15 0.15 0.15], 'EdgeColor', 'k', 'LineWidth', 1.2);
   finish_axes(sprintf('Streamlines around NACA-%s', sol.airfoil), opts, label);

   % Cp contour
   figs(3) = figure('Name', 'Pressure coefficient', 'Color', 'w');
   contourf(xp, zp, Cp, 40, 'LineStyle', 'none');
   hold on
   caxis([max(min(Cp(:)), -3), 1]);
   colormap(parula);
   cb = colorbar;
   cb.Label.String = 'C_p';
   fill(sol.x(1:N+1), sol.z(1:N+1), [0.15 0.15 0.15], 'EdgeColor', 'k', 'LineWidth', 1.2);
   finish_axes(sprintf('Pressure coefficient around NACA-%s', sol.airfoil), opts, label);

   if opts.Save
       names = {'quiver', 'streamlines', 'cp_contour'};
       for i = 1:3
           fname = fullfile(opts.OutDir, sprintf('NACA%s_a%.1f_U%.1f_N%d_%s.png', ...
               sol.airfoil, sol.alpha, sol.u_inf, N, names{i}));
           exportgraphics(figs(i), fname, 'Resolution', 200);
       end
   end
end

function fig = plot_pressure(sol, opts)

   s = sol.surf;

   fig = figure('Name', 'Surface pressure', 'Color', 'w');

   plot(s.x(s.upper), s.Cp(s.upper), '-', 'LineWidth', 1.5, ...
       'Color', [0.10 0.35 0.75], 'DisplayName', 'Upper surface');
   hold on
   plot(s.x(s.lower), s.Cp(s.lower), '-', 'LineWidth', 1.5, ...
       'Color', [0.85 0.15 0.15], 'DisplayName', 'Lower surface');
   yline(0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');

   set(gca, 'YDir', 'reverse');
   xlabel('x/c');
   ylabel('C_p');
   title(sprintf('Pressure distribution, NACA-%s at \\alpha = %.1f\\circ (N = %d)', ...
       sol.airfoil, sol.alpha, sol.N), 'Interpreter', 'tex');
   legend('Location', 'best');
   grid on
   box on
   xlim([-0.02, 1.02]);

   ax2 = axes('Position', [0.55 0.18 0.32 0.16]);
   plot(ax2, sol.x(1:sol.N+1), sol.z(1:sol.N+1), 'k', 'LineWidth', 1.2);
   axis(ax2, 'equal', 'off');

   if opts.Save
       exportgraphics(fig, fullfile(opts.OutDir, ...
           sprintf('NACA%s_a%.1f_N%d_cp.png', sol.airfoil, sol.alpha, sol.N)), ...
           'Resolution', 200);
   end
end

function finish_axes(titleText, opts, label)
   axis equal
   xlim(opts.Xlim);
   ylim(opts.Zlim);
   xlabel('x/c');
   ylabel('z/c');
   title(titleText);
   box on

   text(opts.Xlim(2) - 0.03*diff(opts.Xlim), opts.Zlim(2) - 0.03*diff(opts.Zlim), ...
       label, ...
       'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
       'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'tex', ...
       'FontSize', 8);

   hold off
end

function [x_velocity, z_velocity, Cp_field] = velocity(x_points, z_points, mu, airfoil_x, airfoil_z, u_inf, alpha, N)

   if ~isequal(size(x_points), size(z_points))
       error('velocity:sizeMismatch', 'x_points and z_points must match in size');
   end

   alpha = alpha * pi / 180;

   x_velocity = repmat(u_inf * cos(alpha), size(x_points));
   z_velocity = repmat(u_inf * sin(alpha), size(z_points));

   for k = 1:N+1
       p1 = [airfoil_x(k),   airfoil_z(k)];
       p2 = [airfoil_x(k+1), airfoil_z(k+1)];

       [u_doublet, v_doublet] = cdoublet_vec(x_points, z_points, p1, p2);

       x_velocity = x_velocity + mu(k) * u_doublet;
       z_velocity = z_velocity + mu(k) * v_doublet;
   end

   inside = inpolygon(x_points, z_points, airfoil_x(1:N+1), airfoil_z(1:N+1));
   x_velocity(inside) = NaN;
   z_velocity(inside) = NaN;

   if nargout > 2
       Cp_field = 1 - (x_velocity.^2 + z_velocity.^2) / u_inf^2;
   end
end
