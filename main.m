clear; clc; close all;

% path setup works regardless of current folder
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'functs'));
dataDir = fullfile(thisDir, 'data');

% using defaults for now
airfoil = '2412';
N       = 200;
alpha   = 5;
u_inf   = 10;

% uncomment to change inputs
% airfoil = ask_airfoil();
% N       = ask_number('Number of panels', 200, 8, 5000);
% alpha   = ask_number('Angle of attack in degrees', 5, -30, 30);
% u_inf   = ask_number('Freestream velocity', 10, 1e-6, 1e6);

% 1 = single case, 2 = convergence studies, 3 = verification, 4 = batch
% choice = ask_number('Select an option (1-4)', 1, 0, 4);
choice = 1;

switch choice

    case 1 % single case

        fprintf('Solving...\n');
        sol = solve_case(airfoil, N, alpha, u_inf);

        fprintf('\nCl (from circulation)   = %.5f\n', sol.coef.Cl_kutta);
        fprintf('Cl (from pressure)      = %.5f\n', sol.coef.Cl_press);
        fprintf('agreement between the two = %.3f%%\n', 100*sol.coef.Cl_error_rel);
        fprintf('Cm about quarter chord   = %.5f\n', sol.coef.Cm_quarter);
        fprintf('Cd (should be ~0)        = %.2e\n', sol.coef.Cd_press);
        fprintf('min surface Cp           = %.5f\n', min(sol.surf.Cp));
        fprintf('condition number of A    = %.2e\n', sol.info.cond);
        fprintf('solve time               = %.3f s\n\n', sol.runtime);

        plot_results(sol);


    case 2 % convergence studies

        study_convergence('Airfoil', airfoil, 'XfoilFile', fullfile(dataDir, 'Xfoil.txt'));


    case 3 % automated checks

        verify_solver();

    case 4 % batch of test cases

        cases = { ...
            '0012', 200,  0, 10; ...
            '0012', 200,  8, 10; ...
            '2412', 200,  5, 10; ...
            '2412', 200, 10, 15; ...
            '4424', 200,  0, 10};

        fprintf('\nSolving %d cases...\n', size(cases,1));

        nCases = size(cases, 1);
        labels = cell(1, nCases);
        Cl_batch  = zeros(1, nCases);
        Cm_batch  = zeros(1, nCases);
        minCp_batch = zeros(1, nCases);

        for i = 1:nCases
            sol = solve_case(cases{i,1}, cases{i,2}, cases{i,3}, cases{i,4});

            labels{i} = sprintf('NACA%s, \\alpha=%.0f\\circ', sol.airfoil, sol.alpha);
            Cl_batch(i)    = sol.coef.Cl_kutta;
            Cm_batch(i)    = sol.coef.Cm_quarter;
            minCp_batch(i) = min(sol.surf.Cp);

            plot_results(sol);
        end

        plot_batch_summary(labels, Cl_batch, Cm_batch, minCp_batch);

    case 0
        fprintf('Exiting.\n');
end


%% comparison chart (case 4)

function plot_batch_summary(labels, Cl_batch, Cm_batch, minCp_batch)
   figure('Name', 'Batch comparison', 'Color', 'w');
   tiledlayout(1, 3, 'TileSpacing', 'compact');

   nexttile;
   bar(Cl_batch, 'FaceColor', [0.10 0.35 0.75]);
   set(gca, 'XTickLabel', labels, 'TickLabelInterpreter', 'tex');
   ylabel('C_L'); title('Lift'); grid on; box on

   nexttile;
   bar(Cm_batch, 'FaceColor', [0.85 0.15 0.15]);
   set(gca, 'XTickLabel', labels, 'TickLabelInterpreter', 'tex');
   ylabel('C_m (quarter chord)'); title('Pitching moment'); grid on; box on

   nexttile;
   bar(minCp_batch, 'FaceColor', [0.20 0.60 0.30]);
   set(gca, 'XTickLabel', labels, 'TickLabelInterpreter', 'tex');
   ylabel('min C_p'); title('Peak suction'); grid on; box on

   sgtitle('Batch comparison across cases');
end


%% input helpers (uncomment line 14)

function code = ask_airfoil()
   while true
       code = input('Enter the NACA airfoil (e.g., 2412): ', 's');
       if isempty(code)
           code = '2412';
       end
       code = strtrim(code);
       if length(code) == 4 && all(isstrprop(code, 'digit'))
           m = str2double(code(1));
           p = str2double(code(2));
           if m > 0 && p == 0
               fprintf('Error: camber specified but camber position is zero.\n');
               continue;
           end
           return;
       end
       fprintf('Error: input must be a 4-character string of digits.\n');
   end
end

function val = ask_number(prompt, default, lo, hi)
   while true
       raw = input(sprintf('%s [%g]: ', prompt, default), 's');
       if isempty(raw)
           val = default;
           return;
       end
       val = str2double(raw);
       if isnan(val)
           fprintf('Error: not a number.\n');
       elseif val < lo || val > hi
           fprintf('Error: enter a value between %g and %g.\n', lo, hi);
       else
           return;
       end
   end
end

function tf = ask_yesno(prompt, default)
   if default
       d = 'y';
   else
       d = 'n';
   end
   raw = lower(strtrim(input(sprintf('%s? (y/n) [%s]: ', prompt, d), 's')));
   if isempty(raw)
       raw = d;
   end
   tf = startsWith(raw, 'y');
end
