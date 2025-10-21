%% ROCKET RECOVERY SYSTEM SNATCH FORCE CALCULATOR
% This script calculates snatch forces on recovery hardware during
% parachute deployment and determines safety factors for all components.
% Designed for high-power rocketry (L1/L2/L3) and competition rockets.
%
% Author: [Your Team Name]
% Date: October 2025

clear; clc; close all;

%% ========================================================================
%  USER INPUT SECTION - MODIFY THESE VALUES FOR YOUR ROCKET
%  ========================================================================

% --- ROCKET PARAMETERS ---
rocket_mass_lb = 9.42;           % Rocket mass at landing (without motors), lb
apogee_ft = 5757;                 % Expected apogee altitude, ft

% --- PARACHUTE PARAMETERS ---
drogue_diameter_in = 24;          % Drogue parachute diameter, inches
drogue_Cd = 1.55;                 % Drogue drag coefficient (typical: 1.5-2.2)

main_diameter_in = 60;            % Main parachute diameter, inches
main_Cd = 2.20;                   % Main drag coefficient (typical: 1.5-2.5)
main_deploy_alt_ft = 900;         % Main parachute deployment altitude, ft

% --- DEPLOYMENT ALTITUDE STRATEGY ---
% Choose: 'absolute' or 'percentage'
deployment_strategy = 'absolute';  % 'absolute': use main_deploy_alt_ft as-is
                                   % 'percentage': deploy at % of apogee

target_deploy_percentage = 15;     % If using 'percentage' strategy (e.g., 15%)
                                   % Typical range: 10-20% for competition rockets

% --- TARGET DEPLOYMENT VELOCITY (Alternative Strategy) ---
use_target_velocity = false;       % Set to true to use velocity-based deployment
target_deploy_velocity_fps = 30;   % Target deployment velocity, ft/s

% --- RECOVERY HARDWARE RATINGS (lb) ---
eyebolt_rating_lb = 1300;          % Steel eyebolt capacity
quicklink_rating_lb = 1400;       % Quick link capacity
shockcord_rating_lb = 1500;        % Shock cord breaking strength
swivel_rating_lb = 1300;           % Swivel capacity (if used)

% Component names for plotting
hardware_names = {'Eyebolt', 'Shock Cord', 'Swivel', 'Quick Link'};
hardware_ratings = [eyebolt_rating_lb, shockcord_rating_lb, ...
                    swivel_rating_lb, quicklink_rating_lb];

% --- SHOCK CORD MATERIAL PROPERTIES ---
% Material options: 'kevlar', 'nylon', 'mixed'
shockcord_material = 'kevlar';    % Choose: 'kevlar', 'nylon', or 'mixed'

total_shockcord_length_ft = 15;   % Total shock cord length, ft
nylon_section_length_ft = 0;      % Length of nylon section (if mixed), ft

% Material properties (elongation % and shock absorption efficiency)
% These are empirical values from rocketry testing
kevlar_elongation = 0.04;         % 4% max elongation
kevlar_absorption = 0.95;         % Minimal absorption (95% of force transmitted)

nylon_elongation = 0.25;          % 25% max elongation  
nylon_absorption = 0.65;          % Moderate absorption (65% of force transmitted)

% --- OPENING SHOCK FACTORS ---
% K = 1.5: Ideal slow opening (with deployment bag)
% K = 2.0-2.5: Normal deployment
% K = 3.5-4.0: Aggressive opening
% K = 5.0: Worst case (rapid/violent opening)
K_base_typical = 2.5;             % Base typical opening shock factor
K_base_worst = 5.0;               % Base worst case opening shock factor

% --- SAFETY REQUIREMENTS ---
min_safety_factor = 3.0;          % Minimum recommended safety factor
warning_safety_factor = 2.0;      % Warning threshold

%% ========================================================================
%  CONSTANTS AND CALCULATIONS
%  ========================================================================

% Physical constants
g = 32.174;                       % Gravitational acceleration, ft/s²
rho = 0.002377;                   % Air density at sea level, slugs/ft³

% Convert mass to slugs
rocket_mass_slugs = rocket_mass_lb / g;

% Calculate parachute areas
A_drogue = pi * (drogue_diameter_in / 12 / 2)^2;  % ft²
A_main = pi * (main_diameter_in / 12 / 2)^2;      % ft²

%% ========================================================================
%  TERMINAL VELOCITIES
%  ========================================================================

% Terminal velocity under drogue
Vt_drogue = sqrt((2 * rocket_mass_lb) / (rho * drogue_Cd * A_drogue));

% Terminal velocity under main
Vt_main = sqrt((2 * rocket_mass_lb) / (rho * main_Cd * A_main));

%% ========================================================================
%  DEPLOYMENT VELOCITY ESTIMATION
%  ========================================================================

% Determine actual deployment altitude based on strategy
switch lower(deployment_strategy)
    case 'absolute'
        actual_deploy_alt = main_deploy_alt_ft;
        strategy_description = sprintf('Absolute: %d ft AGL', actual_deploy_alt);
        
    case 'percentage'
        actual_deploy_alt = apogee_ft * (target_deploy_percentage / 100);
        strategy_description = sprintf('%.0f%% of apogee = %d ft AGL', ...
                                       target_deploy_percentage, round(actual_deploy_alt));
    otherwise
        error('Invalid deployment strategy. Choose: absolute or percentage');
end

% Override with velocity-based if requested
if use_target_velocity
    % Calculate what altitude gives us target velocity
    % v = Vt * (1 - exp(-d/L)), solve for d
    % d = -L * ln(1 - v/Vt)
    velocity_fraction_target = target_deploy_velocity_fps / Vt_drogue;
    
    if velocity_fraction_target >= 1.0
        warning('Target velocity >= terminal velocity! Using 95% of terminal velocity.');
        velocity_fraction_target = 0.95;
    end
    
    drop_distance_for_target = -L_char * log(1 - velocity_fraction_target);
    actual_deploy_alt = apogee_ft - drop_distance_for_target;
    strategy_description = sprintf('Velocity-based: %.1f ft/s → %d ft AGL', ...
                                   target_deploy_velocity_fps, round(actual_deploy_alt));
    velocity_fraction = velocity_fraction_target;
else
    % Distance fallen from apogee when main deploys
    drop_distance = apogee_ft - actual_deploy_alt;
    
    % Estimate velocity at main deployment
    % Uses exponential approach to terminal velocity model
    L_char = 250;  % Characteristic length, ft (empirical for drogue descent)
    
    if drop_distance >= 500
        velocity_fraction = 0.95;  % Near terminal velocity after 500 ft
    else
        velocity_fraction = 1 - exp(-drop_distance / L_char);
    end
end

V_main_deploy = velocity_fraction * Vt_drogue;

% Calculate what percentage of apogee this represents
deploy_percentage_actual = (actual_deploy_alt / apogee_ft) * 100;

%% ========================================================================
%  SHOCK CORD MATERIAL ANALYSIS
%  ========================================================================

% Calculate effective shock absorption factor based on material
switch lower(shockcord_material)
    case 'kevlar'
        shock_absorption_factor = kevlar_absorption;
        cord_description = 'All Kevlar (minimal stretch)';
        
    case 'nylon'
        shock_absorption_factor = nylon_absorption;
        cord_description = 'All Nylon (high stretch)';
        
    case 'mixed'
        % Calculate weighted average based on cord lengths
        if nylon_section_length_ft > total_shockcord_length_ft
            warning('Nylon section longer than total cord! Using total length.');
            nylon_section_length_ft = total_shockcord_length_ft;
        end
        
        kevlar_fraction = (total_shockcord_length_ft - nylon_section_length_ft) / total_shockcord_length_ft;
        nylon_fraction = nylon_section_length_ft / total_shockcord_length_ft;
        
        % Weighted average of absorption factors
        % The nylon section has disproportionate effect due to its stretch
        shock_absorption_factor = kevlar_fraction * kevlar_absorption + ...
                                  nylon_fraction * nylon_absorption;
        
        cord_description = sprintf('Mixed: %.1f ft Kevlar + %.1f ft Nylon', ...
                                   total_shockcord_length_ft - nylon_section_length_ft, ...
                                   nylon_section_length_ft);
    otherwise
        error('Invalid shock cord material. Choose: kevlar, nylon, or mixed');
end

% Apply shock absorption to opening shock factors
K_typical = K_base_typical * shock_absorption_factor;
K_worst = K_base_worst * shock_absorption_factor;

% Calculate force reduction percentage
force_reduction_pct = (1 - shock_absorption_factor) * 100;

%% ========================================================================
%  SNATCH FORCE CALCULATIONS
%  ========================================================================

% Steady-state drag force at deployment velocity
F_drag_main = 0.5 * rho * V_main_deploy^2 * main_Cd * A_main;

% Snatch forces
F_snatch_typical = K_typical * F_drag_main;
F_snatch_worst = K_worst * F_drag_main;

% Safety factors for each component
SF_typical = hardware_ratings / F_snatch_typical;
SF_worst = hardware_ratings / F_snatch_worst;

%% ========================================================================
%  DISPLAY RESULTS
%  ========================================================================

fprintf('\n========================================\n');
fprintf('ROCKET RECOVERY SNATCH FORCE ANALYSIS\n');
fprintf('========================================\n\n');

fprintf('ROCKET PARAMETERS:\n');
fprintf('  Mass: %.2f lb\n', rocket_mass_lb);
fprintf('  Apogee: %d ft\n', apogee_ft);
fprintf('\n');

fprintf('DEPLOYMENT STRATEGY:\n');
fprintf('  Strategy: %s\n', strategy_description);
fprintf('  Main deployment altitude: %d ft (%.1f%% of apogee)\n', ...
        round(actual_deploy_alt), deploy_percentage_actual);
fprintf('  Distance fallen: %d ft\n', round(apogee_ft - actual_deploy_alt));
fprintf('\n');

fprintf('PARACHUTE SPECIFICATIONS:\n');
fprintf('  Drogue: %d" diameter, Cd = %.2f, Area = %.2f ft²\n', ...
        drogue_diameter_in, drogue_Cd, A_drogue);
fprintf('  Main:   %d" diameter, Cd = %.2f, Area = %.2f ft²\n', ...
        main_diameter_in, main_Cd, A_main);
fprintf('\n');

fprintf('SHOCK CORD CONFIGURATION:\n');
fprintf('  Material: %s\n', cord_description);
fprintf('  Total length: %.1f ft\n', total_shockcord_length_ft);
fprintf('  Shock absorption factor: %.2f (%.1f%% force reduction)\n', ...
        shock_absorption_factor, force_reduction_pct);
fprintf('  Effective K (typical): %.2f (base: %.2f)\n', K_typical, K_base_typical);
fprintf('  Effective K (worst):   %.2f (base: %.2f)\n', K_worst, K_base_worst);
fprintf('\n');

fprintf('DESCENT VELOCITIES:\n');
fprintf('  Terminal velocity (drogue): %.1f ft/s (%.1f mph)\n', ...
        Vt_drogue, Vt_drogue * 0.681818);
fprintf('  Terminal velocity (main):   %.1f ft/s (%.1f mph)\n', ...
        Vt_main, Vt_main * 0.681818);
fprintf('  Velocity at main deploy:    %.1f ft/s (%.1f mph)\n', ...
        V_main_deploy, V_main_deploy * 0.681818);
fprintf('  (%.0f%% of terminal velocity)\n', velocity_fraction * 100);
fprintf('\n');

fprintf('SNATCH FORCES:\n');
fprintf('  Steady-state drag: %.1f lb\n', F_drag_main);
fprintf('  Typical deployment (K=%.1f): %.1f lb\n', K_typical, F_snatch_typical);
fprintf('  Worst case (K=%.1f): %.1f lb\n', K_worst, F_snatch_worst);
fprintf('\n');

fprintf('========================================\n');
fprintf('HARDWARE SAFETY FACTOR ANALYSIS\n');
fprintf('========================================\n\n');

fprintf('%-15s | Rating | Typical SF | Worst SF | Status\n', 'Component');
fprintf('--------------------------------------------------------------\n');

for i = 1:length(hardware_names)
    % Determine status
    if SF_worst(i) >= min_safety_factor
        status = 'OK';
        symbol = char(10003);  % Check mark
    elseif SF_worst(i) >= warning_safety_factor
        status = 'WARNING';
        symbol = '!';
    else
        status = 'DANGER';
        symbol = 'X';
    end
    
    fprintf('%-15s | %5d  |   %5.2f    |  %5.2f   | %s %s\n', ...
            hardware_names{i}, hardware_ratings(i), ...
            SF_typical(i), SF_worst(i), symbol, status);
end

fprintf('\n');
fprintf('Legend:\n');
fprintf('  %c OK      - Safety factor >= %.1f:1 (even in worst case)\n', ...
        char(10003), min_safety_factor);
fprintf('  ! WARNING - Safety factor %.1f:1 to %.1f:1 (marginal)\n', ...
        warning_safety_factor, min_safety_factor);
fprintf('  X DANGER  - Safety factor < %.1f:1 (upgrade recommended)\n', ...
        warning_safety_factor);
fprintf('\n');

%% ========================================================================
%  VISUALIZATION
%  ========================================================================

% Create figure with multiple subplots
figure('Position', [100, 100, 1400, 900], 'Color', 'w');

% Subplot 1: Safety factors bar chart
subplot(2, 3, 1);
x = 1:length(hardware_names);
bar_data = [SF_typical; SF_worst]';

b = bar(x, bar_data, 'grouped');
b(1).FaceColor = [0.2 0.6 0.8];  % Blue for typical
b(2).FaceColor = [0.8 0.2 0.2];  % Red for worst case

hold on;
yline(min_safety_factor, 'g--', 'LineWidth', 2, ...
      'Label', sprintf('Minimum SF = %.1f', min_safety_factor));
yline(warning_safety_factor, 'y--', 'LineWidth', 2, ...
      'Label', sprintf('Warning = %.1f', warning_safety_factor));

set(gca, 'XTick', x, 'XTickLabel', hardware_names);
ylabel('Safety Factor');
title('Hardware Safety Factors');
legend('Typical (K=2.5)', 'Worst Case (K=5.0)', 'Location', 'northeast');
grid on;
xtickangle(45);

% Subplot 2: Force comparison
subplot(2, 3, 2);
force_data = [F_snatch_typical, F_snatch_worst];
bar_labels = {'Typical', 'Worst Case'};
b2 = bar(force_data);
b2.FaceColor = 'flat';
b2.CData(1,:) = [0.2 0.6 0.8];
b2.CData(2,:) = [0.8 0.2 0.2];

set(gca, 'XTickLabel', bar_labels);
ylabel('Force (lb)');
title('Expected Snatch Forces');
grid on;

% Add value labels on bars
for i = 1:length(force_data)
    text(i, force_data(i), sprintf('%.1f lb', force_data(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% Subplot 3: Safety status visualization
subplot(2, 3, 3);
status_matrix = zeros(length(hardware_names), 2);

for i = 1:length(hardware_names)
    % Typical case status
    if SF_typical(i) >= min_safety_factor
        status_matrix(i, 1) = 3;  % Green
    elseif SF_typical(i) >= warning_safety_factor
        status_matrix(i, 1) = 2;  % Yellow
    else
        status_matrix(i, 1) = 1;  % Red
    end
    
    % Worst case status
    if SF_worst(i) >= min_safety_factor
        status_matrix(i, 2) = 3;
    elseif SF_worst(i) >= warning_safety_factor
        status_matrix(i, 2) = 2;
    else
        status_matrix(i, 2) = 1;
    end
end

imagesc(status_matrix);
colormap([0.8 0.2 0.2; 1 0.8 0; 0.2 0.8 0.2]);  % Red, Yellow, Green
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Typical', 'Worst Case'});
set(gca, 'YTick', 1:length(hardware_names), 'YTickLabel', hardware_names);
title('Component Safety Status');
colorbar('Ticks', [1 2 3], 'TickLabels', {'Danger', 'Warning', 'OK'});

% Subplot 4: Apogee scaling analysis - deployment altitude vs apogee
subplot(2, 3, 4);

% Calculate required deployment altitude for different apogees
% to maintain current deployment velocity
apogee_range = 1000:500:10000;
required_deploy_alt = zeros(size(apogee_range));
deploy_pct = zeros(size(apogee_range));

for i = 1:length(apogee_range)
    % Calculate drop distance needed to reach current deployment velocity
    drop_needed = -L_char * log(1 - velocity_fraction);
    required_deploy_alt(i) = apogee_range(i) - drop_needed;
    deploy_pct(i) = (required_deploy_alt(i) / apogee_range(i)) * 100;
end

yyaxis left
plot(apogee_range, required_deploy_alt, 'b-', 'LineWidth', 2);
ylabel('Deployment Altitude (ft)', 'Color', 'b');
hold on;

% Mark current configuration
plot(apogee_ft, actual_deploy_alt, 'ro', 'MarkerSize', 10, ...
     'MarkerFaceColor', 'r');
text(apogee_ft, actual_deploy_alt, '  Current', 'FontSize', 9);

yyaxis right
plot(apogee_range, deploy_pct, 'r--', 'LineWidth', 1.5);
ylabel('Deployment % of Apogee', 'Color', 'r');
ylim([0 100]);

xlabel('Apogee (ft)');
title(sprintf('Maintain %.1f ft/s Deployment Velocity', V_main_deploy));
grid on;
legend('Deployment Altitude', 'Percentage', 'Location', 'northwest');

% Add annotation box
annotation_text = sprintf(['For consistent forces,\n' ...
                          'deploy at ~%.0f%% of apogee'], ...
                          mean(deploy_pct));
dim = [.32 .45 .1 .1];
annotation('textbox', dim, 'String', annotation_text, ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
           'EdgeColor', 'black', 'FontSize', 8);

% Subplot 5: Shock cord material comparison
subplot(2, 3, 5);

% Calculate forces for all three material types
materials = {'Kevlar', 'Mixed\n(50/50)', 'All Nylon'};
absorption_factors = [kevlar_absorption, ...
                      (kevlar_absorption + nylon_absorption) / 2, ...
                      nylon_absorption];
                      
worst_forces = zeros(1, 3);
for i = 1:3
    worst_forces(i) = K_base_worst * absorption_factors(i) * F_drag_main;
end

b3 = bar(worst_forces);
b3.FaceColor = 'flat';
b3.CData(1,:) = [0.7 0.7 0.7];  % Gray for Kevlar
b3.CData(2,:) = [0.8 0.6 0.2];  % Orange for Mixed
b3.CData(3,:) = [0.2 0.6 0.8];  % Blue for Nylon

set(gca, 'XTickLabel', materials);
ylabel('Worst-Case Force (lb)');
title('Effect of Shock Cord Material');
grid on;

% Add current selection indicator
hold on;
if strcmpi(shockcord_material, 'kevlar')
    current_idx = 1;
elseif strcmpi(shockcord_material, 'nylon')
    current_idx = 3;
else
    current_idx = 2;
end
plot(current_idx, worst_forces(current_idx), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
text(current_idx, worst_forces(current_idx) + max(worst_forces)*0.05, ...
     'Current', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'r');

% Add value labels
for i = 1:3
    text(i, worst_forces(i), sprintf('%.0f lb', worst_forces(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% Subplot 6: Nylon length vs safety factor (for mixed cord)
subplot(2, 3, 6);
nylon_lengths = 0:0.5:total_shockcord_length_ft;
SF_curve = zeros(1, length(nylon_lengths));

for i = 1:length(nylon_lengths)
    kevlar_frac = (total_shockcord_length_ft - nylon_lengths(i)) / total_shockcord_length_ft;
    nylon_frac = nylon_lengths(i) / total_shockcord_length_ft;
    abs_factor = kevlar_frac * kevlar_absorption + nylon_frac * nylon_absorption;
    force = K_base_worst * abs_factor * F_drag_main;
    SF_curve(i) = min(hardware_ratings) / force;
end

plot(nylon_lengths, SF_curve, 'b-', 'LineWidth', 2);
hold on;

% Show current configuration
if strcmpi(shockcord_material, 'mixed')
    current_SF = min(SF_worst);
    plot(nylon_section_length_ft, current_SF, 'ro', 'MarkerSize', 10, ...
         'MarkerFaceColor', 'r');
end

% Add reference lines
yline(min_safety_factor, 'g--', 'LineWidth', 1.5, ...
      'Label', sprintf('Target SF = %.1f', min_safety_factor));
yline(warning_safety_factor, 'y--', 'LineWidth', 1.5, ...
      'Label', sprintf('Minimum = %.1f', warning_safety_factor));

xlabel('Nylon Section Length (ft)');
ylabel('Safety Factor (Weakest Component)');
title('Nylon Length Effect on Safety Factor');
grid on;
legend('SF vs Nylon Length', 'Location', 'southeast');

% Add text annotation
xlim([0, total_shockcord_length_ft]);
text_str = sprintf('Total cord: %.1f ft', total_shockcord_length_ft);
text(total_shockcord_length_ft * 0.5, max(SF_curve) * 0.9, text_str, ...
     'HorizontalAlignment', 'center', 'FontSize', 9, ...
     'BackgroundColor', 'w', 'EdgeColor', 'k');

% Add overall title
sgtitle('Rocket Recovery System Analysis', 'FontSize', 16, 'FontWeight', 'bold');

%% ========================================================================
%  RECOMMENDATIONS
%  ========================================================================

fprintf('\n========================================\n');
fprintf('APOGEE SCALING ANALYSIS\n');
fprintf('========================================\n\n');

fprintf('IMPORTANT: Deployment altitude must scale with apogee!\n\n');

% Test what happens at higher apogees with current absolute altitude
test_apogees = [apogee_ft, 3000, 5500, 10000];
fprintf('Impact of increasing apogee (with %s deployment):\n\n', ...
        lower(deployment_strategy));

fprintf('Apogee    Deploy Alt   Distance   Velocity   Worst Force   Min SF\n');
fprintf('------    ----------   Fallen     at Deploy  -----------   ------\n');

for test_apo = test_apogees
    % Calculate for current strategy
    if strcmpi(deployment_strategy, 'percentage')
        test_deploy = test_apo * (target_deploy_percentage / 100);
    else
        test_deploy = main_deploy_alt_ft;  % Absolute stays same
    end
    
    test_drop = test_apo - test_deploy;
    
    if test_drop >= 500
        test_v_frac = 0.95;
    else
        test_v_frac = 1 - exp(-test_drop / L_char);
    end
    
    test_v = test_v_frac * Vt_drogue;
    test_F_drag = 0.5 * rho * test_v^2 * main_Cd * A_main;
    test_F_worst = K_worst * test_F_drag;
    test_min_SF = min(hardware_ratings) / test_F_worst;
    
    % Flag if safety factor too low
    flag = '';
    if test_min_SF < warning_safety_factor
        flag = '  *** DANGER ***';
    elseif test_min_SF < min_safety_factor
        flag = '  ! Warning';
    end
    
    fprintf('%5d ft   %6d ft   %5d ft   %5.1f ft/s   %6.1f lb     %.2f%s\n', ...
            round(test_apo), round(test_deploy), round(test_drop), ...
            test_v, test_F_worst, test_min_SF, flag);
end

fprintf('\n');

if strcmpi(deployment_strategy, 'absolute')
    fprintf('⚠️  WARNING: You are using ABSOLUTE deployment altitude!\n');
    fprintf('    At higher apogees, forces will INCREASE significantly.\n');
    fprintf('    Recommendation: Switch to PERCENTAGE strategy\n');
    fprintf('    Set: deployment_strategy = ''percentage''\n');
    fprintf('         target_deploy_percentage = %.0f%%\n', deploy_percentage_actual);
elseif strcmpi(deployment_strategy, 'percentage')
    fprintf('✓ GOOD: You are using PERCENTAGE deployment strategy.\n');
    fprintf('   Forces will remain consistent across different apogees.\n');
end

fprintf('\n========================================\n');
fprintf('RECOMMENDATIONS\n');
fprintf('========================================\n\n');

% Find weakest component
[min_SF_worst, idx_min] = min(SF_worst);

if min_SF_worst < warning_safety_factor
    fprintf('*** CRITICAL: %s has insufficient safety factor (%.2f:1) ***\n', ...
            hardware_names{idx_min}, min_SF_worst);
    fprintf('\nSuggested actions (in order of effectiveness):\n');
    fprintf('  1. Increase main deployment altitude (currently %d ft, %.1f%% of apogee)\n', ...
            round(actual_deploy_alt), deploy_percentage_actual);
    fprintf('  2. Upgrade %s to higher rating\n', hardware_names{idx_min});
    fprintf('  3. Use deployment bag to reduce opening shock factor\n');
    
    % Calculate how much nylon needed to achieve target SF
    if strcmpi(shockcord_material, 'kevlar')
        fprintf('\n  *** SHOCK CORD MATERIAL RECOMMENDATION ***\n');
        fprintf('  Current: All Kevlar (no shock absorption)\n');
        
        target_SF = min_safety_factor;
        target_force = hardware_ratings(idx_min) / target_SF;
        target_absorption = target_force / (K_base_worst * F_drag_main);
        
        if target_absorption > nylon_absorption && target_absorption < kevlar_absorption
            % Calculate required nylon length
            % target_absorption = kevlar_frac * kevlar_abs + nylon_frac * nylon_abs
            % target_absorption = (1 - nylon_frac) * kevlar_abs + nylon_frac * nylon_abs
            % Solve for nylon_frac
            nylon_frac_needed = (kevlar_absorption - target_absorption) / ...
                                (kevlar_absorption - nylon_absorption);
            nylon_length_needed = nylon_frac_needed * total_shockcord_length_ft;
            
            fprintf('\n  → Add %.1f ft of NYLON to your shock cord\n', nylon_length_needed);
            fprintf('    (%.0f%% of total %0.1f ft shock cord length)\n', ...
                    nylon_frac_needed * 100, total_shockcord_length_ft);
            fprintf('    This will reduce worst-case force to %.1f lb\n', target_force);
            fprintf('    New safety factor: %.2f:1\n', target_SF);
            
            % Show configuration suggestion
            kevlar_length_needed = total_shockcord_length_ft - nylon_length_needed;
            fprintf('\n  Suggested configuration:\n');
            fprintf('    - Main parachute end: %.1f ft Nylon (energy absorber)\n', ...
                    nylon_length_needed * 0.7);
            fprintf('    - Middle section: %.1f ft Kevlar (durable)\n', kevlar_length_needed);
            fprintf('    - Rocket attachment: %.1f ft Nylon (shock absorber)\n', ...
                    nylon_length_needed * 0.3);
            
        elseif target_absorption <= nylon_absorption
            fprintf('\n  → Switch to ALL NYLON shock cord\n');
            fprintf('    This may not be sufficient alone.\n');
            fprintf('    Also consider deployment bag or higher deployment altitude.\n');
        else
            fprintf('\n  → Adding nylon alone will not achieve required safety factor\n');
            fprintf('    Must also increase deployment altitude or upgrade hardware\n');
        end
    elseif strcmpi(shockcord_material, 'mixed')
        fprintf('\n  *** INCREASE NYLON SECTION ***\n');
        fprintf('  Current nylon length: %.1f ft\n', nylon_section_length_ft);
        
        % Calculate additional nylon needed
        target_SF = min_safety_factor;
        target_force = hardware_ratings(idx_min) / target_SF;
        target_absorption = target_force / (K_base_worst * F_drag_main);
        
        if target_absorption > nylon_absorption && target_absorption < kevlar_absorption
            nylon_frac_needed = (kevlar_absorption - target_absorption) / ...
                                (kevlar_absorption - nylon_absorption);
            nylon_length_needed = nylon_frac_needed * total_shockcord_length_ft;
            additional_nylon = nylon_length_needed - nylon_section_length_ft;
            
            if additional_nylon > 0
                fprintf('  → Add %.1f ft MORE nylon (total: %.1f ft)\n', ...
                        additional_nylon, nylon_length_needed);
                fprintf('    New configuration: %.1f ft Kevlar + %.1f ft Nylon\n', ...
                        total_shockcord_length_ft - nylon_length_needed, nylon_length_needed);
            else
                fprintf('  → Current nylon section is sufficient\n');
                fprintf('    Consider other options like deployment bag\n');
            end
        end
    end
    
    % Calculate recommended altitude for SF = 3.0
    target_force = hardware_ratings(idx_min) / min_safety_factor;
    target_drag = target_force / K_worst;
    target_velocity = sqrt((2 * target_drag) / (rho * main_Cd * A_main));
    target_v_fraction = target_velocity / Vt_drogue;
    
    if target_v_fraction < 1.0 && target_v_fraction > 0
        target_drop = -L_char * log(1 - target_v_fraction);
        recommended_alt = apogee_ft - target_drop;
        fprintf('\n  Alternative: Deploy at %.0f ft or higher\n', recommended_alt);
    end
    
elseif min_SF_worst < min_safety_factor
    fprintf('! WARNING: %s has marginal safety factor (%.2f:1)\n', ...
            hardware_names{idx_min}, min_SF_worst);
    fprintf('\nConsider:\n');
    fprintf('  1. Using deployment bag to control opening rate\n');
    fprintf('  2. Slightly increasing deployment altitude\n');
    
    if strcmpi(shockcord_material, 'kevlar')
        fprintf('  3. Add 2-3 ft of nylon to shock cord for energy absorption\n');
    end
    
    if strcmpi(deployment_strategy, 'absolute')
        fprintf('  4. Switch to PERCENTAGE strategy for multi-motor flexibility\n');
    end
else
    fprintf('✓ All components have adequate safety factors!\n');
    fprintf('\nCurrent configuration is suitable for flight.\n');
    
    % Show benefit of current material choice
    if ~strcmpi(shockcord_material, 'kevlar')
        fprintf('\nYour shock cord material is helping!\n');
        fprintf('  Force reduction from shock absorption: %.1f%%\n', force_reduction_pct);
        fprintf('  Without shock absorption, worst-case force would be: %.1f lb\n', ...
                K_base_worst * F_drag_main);
    end
end

fprintf('\n----------------------------------------\n');
fprintf('COMPETITION ROCKET RECOMMENDATIONS:\n');
fprintf('----------------------------------------\n');
fprintf('For rockets flying multiple motors at different altitudes:\n\n');
fprintf('  Strategy 1: PERCENTAGE deployment (Recommended)\n');
fprintf('    - Set deployment_strategy = ''percentage''\n');
fprintf('    - Set target_deploy_percentage = %.0f%%\n', deploy_percentage_actual);
fprintf('    - Forces remain consistent across all flights\n');
fprintf('    - No reprogramming needed for different motors\n\n');

fprintf('  Strategy 2: VELOCITY-based deployment\n');
fprintf('    - Set use_target_velocity = true\n');
fprintf('    - Set target_deploy_velocity_fps = %.1f\n', V_main_deploy);
fprintf('    - Deployment altitude auto-adjusts for each apogee\n');
fprintf('    - Requires altimeter with this capability\n\n');

fprintf('  Strategy 3: ABSOLUTE altitude (Not recommended for multi-motor)\n');
fprintf('    - Only suitable if flying same motor repeatedly\n');
fprintf('    - Must recalculate for each new motor/apogee combination\n');

fprintf('\n========================================\n');
fprintf('END OF ANALYSIS\n');
fprintf('========================================\n');