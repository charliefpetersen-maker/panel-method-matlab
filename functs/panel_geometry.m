function geom = panel_geometry(airfoil_x, airfoil_z, N)
   % works out panel midpoints, lengths, angles, tangents and normals
   % panel 1 = TE on lower surface, panel N = TE on upper surface, N+1 = wake

   dx = diff(airfoil_x(1:N+1)).';
   dz = diff(airfoil_z(1:N+1)).';

   geom.len  = hypot(dx, dz);
   geom.beta = atan2(dz, dx);

   geom.midx = (airfoil_x(1:N).' + airfoil_x(2:N+1).') / 2;
   geom.midz = (airfoil_z(1:N).' + airfoil_z(2:N+1).') / 2;

   geom.tx = cos(geom.beta);
   geom.tz = sin(geom.beta);

   % outward normal, points away from the body given the TE->LE->TE ordering
   geom.nx = -sin(geom.beta);
   geom.nz =  cos(geom.beta);

   cumLen = [0; cumsum(geom.len)];
   geom.s  = cumLen(1:N) + geom.len/2;
end
