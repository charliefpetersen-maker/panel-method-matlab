function pass = verify_solver()
   % runs a bunch of checks on the solver,
   fprintf('\n=== Panel method verification ===\n\n');

   results = cell(0,3);

   % test 1
   rng(0);
   p1 = [0.1, -0.05];
   p2 = [0.35, 0.02];
   px = randn(20,1)*0.5 + 0.5;
   pz = randn(20,1)*0.3;
   [uv, vv] = cdoublet_vec(px, pz, p1, p2);
   us = zeros(20,1); vs = zeros(20,1);
   for i = 1:20
       [us(i), vs(i)] = cdoublet_scalar([px(i), pz(i)], p1, p2);
   end
   e1 = max(max(abs(uv - us)), max(abs(vv - vs)));
   results(end+1,:) = check('Vectorised influence matches scalar', e1, 1e-10, '%.2e');

   % test 2
   sol = solve_case('0012', 200, 0, 10);
   e2 = abs(sol.coef.Cl);
   results(end+1,:) = check('Symmetric aerofoil at 0 deg gives zero lift', e2, 1e-6, '%.2e');

   % test 3
   sp = solve_case('0012', 200,  5, 10);
   sm = solve_case('0012', 200, -5, 10);
   e3 = abs(sp.coef.Cl + sm.coef.Cl);
   results(end+1,:) = check('Lift antisymmetric in alpha', e3, 1e-6, '%.2e');

   % test 4 - a real 12% thick aerofoil sits above 2*pi by roughly 0.77*t
   a = [-2, 0, 2, 4];
   Cl = zeros(size(a));
   for i = 1:numel(a)
       s0 = solve_case('0012', 400, a(i), 10);
       Cl(i) = s0.coef.Cl;
   end
   pf = polyfit(a*pi/180, Cl, 1);
   e4 = abs(pf(1) - 2*pi) / (2*pi);
   results(end+1,:) = check('Lift slope within 10% of 2*pi per radian', e4, 0.10, '%.3f');

   % test 5
   sol = solve_case('2412', 200, 6, 10);
   e5 = sol.info.maxVn;
   results(end+1,:) = check('Residual normal velocity below 0.5% of U_inf', e5, 5e-3, '%.2e');

   % test 6
   N = sol.N;
   e6 = abs(sol.mu(1) - sol.mu(N) + sol.mu(N+1)) / sol.u_inf;
   results(end+1,:) = check('Kutta condition satisfied', e6, 1e-6, '%.2e');

   % test 7
   e7 = sol.coef.Cl_error_rel;
   results(end+1,:) = check('Circulation and pressure lift agree within 8%', e7, 0.08, '%.4f');

   % test 8
   e8 = abs(sol.coef.Cd_press);
   results(end+1,:) = check('Inviscid drag is negligible', e8, 0.02, '%.4f');

   % test 9
   [x, z] = panelgen('2412', 200, 0);
   e9 = hypot(x(1) - x(N+1), z(1) - z(N+1));
   results(end+1,:) = check('Trailing edge is closed', e9, 1e-12, '%.2e');

   % test 10
   wState = warning('off', 'panel_strengths:noGauge');
   [~, infoFree] = panel_strengths(x, z, 10, 5, N, 'FixGauge', false);
   warning(wState);
   [~, infoFixed] = panel_strengths(x, z, 10, 5, N, 'FixGauge', true);
   fprintf('   cond: unconstrained %.2e -> gauge-fixed %.2e\n', infoFree.cond, infoFixed.cond);
   e10 = infoFixed.cond;
   results(end+1,:) = check('Gauge-fixed system is well conditioned', e10, 1e10, '%.2e');

   fprintf('\n%-52s %12s %8s\n', 'Test', 'Value', 'Result');
   fprintf('%s\n', repmat('-', 1, 74));
   pass = true;
   for i = 1:size(results, 1)
       fprintf('%-52s %12s %8s\n', results{i,1}, results{i,2}, results{i,3});
       pass = pass && strcmp(results{i,3}, 'PASS');
   end
   fprintf('%s\n', repmat('-', 1, 74));

   if pass
       fprintf('\nAll %d checks passed.\n\n', size(results,1));
   else
       fprintf('\nOne or more checks FAILED.\n\n');
   end
end

function row = check(name, value, tol, fmt)
   if value <= tol
       verdict = 'PASS';
   else
       verdict = 'FAIL';
   end
   row = {name, sprintf(fmt, value), verdict};
end

function [u,v] = cdoublet_scalar(p,p1,p2)
   % scalar-only reference implementation, kept separate from cdoublet_vec
   % purely so test 1 can check the vectorised version against something
   % independently written

   alfa=-atan2(p2(2)-p1(2),p2(1)-p1(1));
   T=[cos(alfa),-sin(alfa);sin(alfa),cos(alfa)];
   A=T*(p-p1)';
   B=T*(p-p2)';
   d1=A(1);
   d2=B(1);
   dz=A(2);

   if abs(dz)<0.000001
       u1=0;
       v1=(d1/(d1^2)-d2/(d2^2))/(2*pi);
   else
       u1=-dz*(1/(d1^2+dz^2)-1/(d2^2+dz^2))/(2*pi);
       v1=(d1/(d1^2+dz^2)-d2/(d2^2+dz^2))/(2*pi);
   end

   temp=[cos(alfa),sin(alfa);-sin(alfa),cos(alfa)]*[u1;v1];
   u=temp(1);
   v=temp(2);
end
