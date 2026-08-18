function [airfoil_x, airfoil_z, geom] = panelgen(airfoil, N, alpha, opts)
   % builds panel coords for a NACA 4 digit aerofoil + one wake panel
   arguments
       airfoil
       N
       alpha
       opts.WakeLength = 0.9
       opts.ClosedTE = true
   end

   if ~all(isstrprop(airfoil, 'digit'))
       error('panelgen:badCode', 'Aerofoil code must be four digits.');
   end

   % NACA digits
   m = str2double(airfoil(1)) / 100;
   p = str2double(airfoil(2)) / 10;
   t = str2double(airfoil(3:4)) / 100;

   if m > 0 && p == 0
       error('panelgen:badCamber', 'camber position cant be zero if there is camber');
   end

   alpha_rad = alpha * pi / 180;

   % cosine spacing, bunches panels near the leading and trailing edges
   idx = (0:N).';
   x_coords = (1 - 0.5 * (1 - cos(2*pi*idx/N))).';

   % find leading edge, first point where x starts increasing again
   leading_edge = find(diff(x_coords) > 0, 1) + 1;
   if isempty(leading_edge)
       leading_edge = 1;
   end

   % mean camber line
   y_camber = zeros(1, N+1);
   theta    = zeros(1, N+1);
   if m > 0
       fore = x_coords < p;
       aft  = ~fore;

       y_camber(fore) = (m / p^2) * (2*p*x_coords(fore) - x_coords(fore).^2);
       theta(fore)    = atan(2*m / p^2 * (p - x_coords(fore)));

       y_camber(aft) = (m / (1-p)^2) * (1 - 2*p + 2*p*x_coords(aft) - x_coords(aft).^2);
       theta(aft)    = atan(2*m / (1-p)^2 * (p - x_coords(aft)));
   end

   % thickness distribution
   c4 = 0.1036;
   if ~opts.ClosedTE
       c4 = 0.1015;
   end
   y_t = 5*t * (0.2969*sqrt(x_coords) - 0.1260*x_coords - 0.3516*x_coords.^2 ...
       + 0.2843*x_coords.^3 - c4*x_coords.^4);

   % combine camber + thickness, upper surface before the leading edge,
   % lower surface after
   isUpper = false(1, N+1);
   isUpper(1:leading_edge-1) = true;

   airfoil_x = zeros(1, N+1);
   airfoil_z = zeros(1, N+1);
   airfoil_x(isUpper)  = x_coords(isUpper)  + y_t(isUpper)  .* sin(theta(isUpper));
   airfoil_z(isUpper)  = y_camber(isUpper)  - y_t(isUpper)  .* cos(theta(isUpper));
   airfoil_x(~isUpper) = x_coords(~isUpper) - y_t(~isUpper) .* sin(theta(~isUpper));
   airfoil_z(~isUpper) = y_camber(~isUpper) + y_t(~isUpper) .* cos(theta(~isUpper));

   % force TE closed
   te_x = mean([airfoil_x(1), airfoil_x(end)]);
   te_z = mean([airfoil_z(1), airfoil_z(end)]);
   airfoil_x([1, end]) = te_x;
   airfoil_z([1, end]) = te_z;

   % wake panel trails off the TE in freestream direction
   airfoil_x(end+1) = te_x + opts.WakeLength * cos(alpha_rad);
   airfoil_z(end+1) = te_z + opts.WakeLength * sin(alpha_rad);

   if nargout > 2
       geom = panel_geometry(airfoil_x, airfoil_z, N);
   end
end
