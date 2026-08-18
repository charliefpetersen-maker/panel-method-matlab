function [srf, coef] = surface_pressure(mu, airfoil_x, airfoil_z, u_inf, alpha, N, opts)
   % gets the surface pressure distribution from the doublet strengths
   arguments
       mu
       airfoil_x
       airfoil_z
       u_inf
       alpha
       N
       opts.OffsetFactor = 1e-4  
       opts.MomentRef = 0.25     % quarter chord
   end

   alpha_rad = alpha * pi / 180;
   geom = panel_geometry(airfoil_x, airfoil_z, N);

   delta = opts.OffsetFactor * geom.len;
   px = geom.midx + delta .* geom.nx;
   pz = geom.midz + delta .* geom.nz;

   u = repmat(u_inf * cos(alpha_rad), N, 1);
   w = repmat(u_inf * sin(alpha_rad), N, 1);

   for k = 1:N+1
       p1 = [airfoil_x(k),   airfoil_z(k)];
       p2 = [airfoil_x(k+1), airfoil_z(k+1)];
       [ud, vd] = cdoublet_vec(px, pz, p1, p2);
       u = u + mu(k) * ud;
       w = w + mu(k) * vd;
   end

   Vt = u .* geom.tx + w .* geom.tz;

   srf.geom = geom;
   srf.x    = geom.midx;
   srf.z    = geom.midz;
   srf.Vt   = Vt;
   srf.Cp   = 1 - (Vt / u_inf).^2;
   srf.offsetFactor = opts.OffsetFactor;

   [~, leIdx] = min(geom.midx);
   srf.lower = (1:N).' <= leIdx;
   srf.upper = ~srf.lower;

   if nargout > 1
       coef = local_aero_coefficients(mu, srf, u_inf, alpha, N, opts.MomentRef);
   end
end

function coef = local_aero_coefficients(mu, srf, u_inf, alpha, N, momentRef)
   % works out Cl, Cd, Cm from mu and Cp distribution

   alpha_rad = alpha * pi / 180;
   geom = srf.geom;
   Cp   = srf.Cp;

   % lift from the circulation - kutta-joukowski
   coef.Cl_kutta = -2 * mu(N+1) / u_inf;

   dFx = -Cp .* geom.nx .* geom.len;
   dFz = -Cp .* geom.nz .* geom.len;

   Fx = sum(dFx);
   Fz = sum(dFz);

   % Fx, Fz are the aerofoil-axis force coeffs
   coef.Cn = Fz;
   coef.Ca = Fx;

   coef.Cl_press = -Fx * sin(alpha_rad) + Fz * cos(alpha_rad);
   coef.Cd_press =  Fx * cos(alpha_rad) + Fz * sin(alpha_rad);

   xr = geom.midx - momentRef;
   zr = geom.midz;
   coef.Cm_quarter = sum(-xr .* dFz + zr .* dFx);

   coef.Cl_error = abs(coef.Cl_kutta - coef.Cl_press);
   if abs(coef.Cl_kutta) > 1e-6
       coef.Cl_error_rel = coef.Cl_error / abs(coef.Cl_kutta);
   else
       coef.Cl_error_rel = NaN;
   end

   coef.Cl = coef.Cl_kutta;
end
