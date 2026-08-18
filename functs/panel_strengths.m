function [mu, info] = panel_strengths(airfoil_x, airfoil_z, u_inf, alpha, N, opts)
   % sets up and solves for doublet strength at each panel
   arguments
       airfoil_x
       airfoil_z
       u_inf
       alpha
       N
       opts.FixGauge = true
   end

   if numel(airfoil_x) < N + 2
       error('panel_strengths:tooFewPoints', 'not enough coordinates for N panels');
   end

   alpha = alpha * pi / 180;

   geom = panel_geometry(airfoil_x, airfoil_z, N);
   beta = geom.beta;

   A = zeros(N+1, N+1);
   B = zeros(N+1, 1);
   B(1:N) = -u_inf * sin(alpha - beta);

   for j = 1:N+1
       p1 = [airfoil_x(j),   airfoil_z(j)];
       p2 = [airfoil_x(j+1), airfoil_z(j+1)];

       [u, v] = cdoublet_vec(geom.midx, geom.midz, p1, p2);

       A(1:N, j) = v .* cos(beta) - u .* sin(beta);
   end

   % kutta condition
   A(N+1, :)   = 0;
   A(N+1, 1)   =  1;
   A(N+1, N)   = -1;
   A(N+1, N+1) =  1;
   B(N+1)      =  0;

   if opts.FixGauge
       gaugeRow      = zeros(1, N+1);
       gaugeRow(1:N) = 1;
       Asolve = [A; gaugeRow];
       Bsolve = [B; 0];
   else
       Asolve = A;
       Bsolve = B;
   end

   mu = Asolve \ Bsolve;

   if nargout > 1
       info.A     = A;
       info.B     = B;
       info.geom  = geom;
       info.gauge = opts.FixGauge;
       info.cond     = cond(Asolve);
       info.residual = norm(A*mu - B);

       vn = (A(1:N, :) * mu - B(1:N)) / u_inf;
       info.maxVn = max(abs(vn));

       nullVec = [ones(N,1); 0];
       info.nullDefect = norm(A * nullVec) / max(norm(A(:)), eps);

       if ~opts.FixGauge
           warning('panel_strengths:noGauge', ...
               'gauge disabled, matrix has a null space, cond = %.2e', info.cond);
       elseif info.cond > 1e12
           warning('panel_strengths:illConditioned', ...
               'condition number is %.2e, check for bad panels', info.cond);
       end
   end
end
