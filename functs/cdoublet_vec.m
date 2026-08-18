function [u, v] = cdoublet_vec(px, pz, p1, p2)
   % velocity induced at (px, pz) by a constant-strength doublet panel
   % running from p1 to p2. equivalent to a vortex sheet, standard panel
   % method building block

   alpha = -atan2(p2(2) - p1(2), p2(1) - p1(1));
   ca = cos(alpha);
   sa = sin(alpha);

   dx1 = px - p1(1);   dz1 = pz - p1(2);
   dx2 = px - p2(1);   dz2 = pz - p2(2);

   d1 = ca .* dx1 - sa .* dz1;
   d2 = ca .* dx2 - sa .* dz2;
   dz = sa .* dx1 + ca .* dz1;

   r1 = d1.^2 + dz.^2;
   r2 = d2.^2 + dz.^2;

   % avoid dividing by zero at panel end point
   r1(r1 < eps) = eps;
   r2(r2 < eps) = eps;

   u1 = -dz .* (1 ./ r1 - 1 ./ r2) / (2*pi);
   v1 = (d1 ./ r1 - d2 ./ r2) / (2*pi);

   % point sits on the panel's own line, general formula is singular there
   onPanel = abs(dz) < 1e-6;
   if any(onPanel(:))
       dd1 = d1(onPanel);
       dd2 = d2(onPanel);
       dd1(abs(dd1) < eps) = eps;
       dd2(abs(dd2) < eps) = eps;
       u1(onPanel) = 0;
       v1(onPanel) = (1 ./ dd1 - 1 ./ dd2) / (2*pi);
   end

   u =  ca .* u1 + sa .* v1;
   v = -sa .* u1 + ca .* v1;
end
