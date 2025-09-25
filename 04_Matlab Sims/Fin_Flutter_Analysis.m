% Rocket fin flutter analysis v3 (hardened)
% Based on John Bennett / Apogee "Peak of Flight" #615 (12/19/2023)
% Validated to ~Mach 1.5 per Bennett; use margin (>= 1.5) and higher-fidelity tools above that.
% Limitations: no tapered thickness; no explicit tip-to-tip modeling (see note on raising G).
% By Zach Pedrick (orig), corrections & hardening by JL/ChatGPT (2025-09-25)

close all; clear; clc;

%% -------------------------- Material Shear Modulus [psi] --------------------------
Balsa_Wood     =  33359;
Birch_Plywood  =  89000;
G10_Fiberglass = 600000;
Carbon_Fiber   = 600000;   % placeholder; depends on layup; see Bennett #615
Aluminum_6061  = 3800000;
Titanium_6M4V  = 6200000;
Steel_4130     =12000000;
PETG           =   9000;
PLA            = 348000;

%% ------------------------------- USER INPUTS --------------------------------------
format = 'RASAero';     % 'RASAero' or 'Openrocket'

% Launch-site elevation [ft]
GL = 2500;              % Spaceport America ~4600 ft for IREC

% Fin geometry (inches) and material
G  = G10_Fiberglass;    % Shear modulus [psi]
Rc = 16;                % Root chord [in]
Tc =  4;                % Tip chord [in]
m  =  6;              % Semi-span / height [in]
Sl = 10;                % Sweep length LE->LE along body/chord [in]
t  =  0.3125;             % Thickness [in]

% Optional: override area if planform is non-trapezoid (set [] to use trapezoid)
Fin_A_override = [];    % [in^2] e.g., from CAD; leave [] for trapezoid area

% Optional: approximate tip-to-tip reinforcement by boosting G (rule-of-thumb)
% G = 2*G;  % <- uncomment if you have aggressive tip-to-tip skins (per Bennett guidance)

%% ---------------------------- Select sim CSV file ---------------------------------
[file,location] = uigetfile({'*.csv;*.CSV'},'Select a Simulation File');
if isequal(file,0), error('User canceled file selection.'); end
fullpath = fullfile(location,file);
fprintf('User selected: %s\n', fullpath);

%% ------------------------- Robust CSV import (RASAero/OR) -------------------------
% Try header-aware read first; fall back to numeric columns if needed.
T = [];
try
    T = readtable(fullpath, 'FileType','text', 'TextType','string');
catch
    % Some older MATLABs struggle with BOM/encoding; try readmatrix later
end

vel = []; alt = [];
fmt = lower(strtrim(format));
if ~isempty(T)
    switch fmt
        case 'rasaero'
            % Common RASAero headers (adjust if your CSV differs)
            candV = ["Total Velocity (ft/sec)","Total Velocity (ft/s)","Total Velocity (fps)"];
            candH = ["Altitude AGL (ft)","Altitude Above Ground (ft)","Altitude (AGL ft)"];
        case 'openrocket'
            candV = ["Total velocity (ft/s)","Velocity (ft/s)","Total velocity [ft/s]"];
            candH = ["Altitude above ground (ft)","Altitude (above ground) [ft]","Altitude [ft]"];
        otherwise
            error("Format must be 'RASAero' or 'Openrocket'.");
    end

    % helper to pick the first matching header
    pickcol = @(cands) T.(cands(find(ismember(cands, string(T.Properties.VariableNames)),1)));
    try
        vel = pickcol(candV);
        alt = pickcol(candH);
    catch
        % Will fall back below
    end
end

if isempty(vel) || isempty(alt)
    % Fallback to numeric column indices used in the original script
    M = readmatrix(fullpath);
    if strcmpi(format,'RASAero')
        vel = M(:,18);         % Total velocity [ft/s]
        alt = M(:,23);         % Altitude AGL [ft]
    elseif strcmpi(format,'Openrocket')
        M   = M(~all(isnan(M),2),:);
        vel = M(:,5);          % Total velocity [ft/s]
        alt = M(:,2);          % Altitude AGL [ft]
    end
end

% Column-vector hygiene
vel = vel(:); alt = alt(:);

% Trim to apogee (remove descent)
[~, apogee_index] = max(alt);
idx = 1:apogee_index;
vel = vel(idx);
alt = alt(idx);

if any(~isfinite(vel)) || any(~isfinite(alt))
    error('Non-finite values found in imported data. Check your CSV export.');
end

%% ------------------------------- Sanity checks ------------------------------------
assert(all([Rc Tc m Sl t] > 0),'Fin geometry must be positive.');
assert(G > 0,'Shear modulus G must be positive.');
if strcmpi(format,'Openrocket')
    warnstate = 'off'; %#ok<NASGU> % OpenRocket files must be exported in Imperial units to match this script
end

%% ------------------------- Atmosphere (ISA) vs altitude ---------------------------
% ISA expects meters; alt+GL converts AGL->MSL; then convert outputs to ft/s and psi.
[~, a, P, ~] = atmosisa( (alt + GL) * 0.3048 ); % requires Aerospace Toolbox
a = a * 3.28084;   % m/s -> ft/s
P = P / 6894.757;  % Pa  -> psi
k  = 1.4;
P0 = 14.696;       % sea-level static pressure [psi]

%% ------------------------- Fin centroid & derived terms ---------------------------
% **** Correct centroid formula uses Sl (sweep), not m (span). ****
Cx = ((2*Tc*Sl) + (Tc^2) + (Sl*Rc) + (Tc*Rc) + (Rc^2)) / (3*(Tc+Rc));  % [in]
E  = (Cx/Rc) - 0.25;                     % epsilon: CG location behind 1/4-chord
DN = (24*E*k*P0) / pi;                   % denominator constant

% Planform area (trapezoid) unless overridden
if isempty(Fin_A_override)
    Fin_A = 0.5 * m * (Tc + Rc);         % [in^2]
else
    Fin_A = Fin_A_override;
end

thic_ratio = t / Rc;                     % t/c (using root chord per Bennett method)
lambda     = Tc / Rc;                    % taper ratio
A          = (m^2) / Fin_A;              % aspect ratio

%% ---------------------------- Flutter computation --------------------------------
Vf = zeros(size(alt));
for i = 1:numel(alt)
    Vf(i) = a(i) * sqrt( G / ( ((DN*A^3)/(thic_ratio^3*(A+2))) * ((lambda+1)/2) * (P(i)/P0) ) );
end

% Pull values of interest
[vel_max, vel_max_index] = max(vel);
alt_vmax   = alt(vel_max_index);
Vf_atvmax  = Vf(vel_max_index);
SF         = Vf_atvmax / vel_max;

% Local Mach numbers at that altitude
M_max     = vel_max   / a(vel_max_index);
M_flutter = Vf_atvmax / a(vel_max_index);

%% --------------------------------- Reporting --------------------------------------
fprintf('Maximum Rocket Velocity occurs at %.2f ft AGL\n', alt_vmax);
fprintf('Maximum Rocket Velocity: %.2f ft/s (Mach %.2f)\n', vel_max, M_max);
fprintf('Fin Flutter Velocity:    %.2f ft/s (Mach %.2f)\n', Vf_atvmax, M_flutter);
fprintf('Fin Flutter Safety Factor (Vf/Vmax): %.3f\n', SF);
if SF < 1.5
    warning('Flutter SF = %.3f (<1.5). Increase thickness/G, reduce span, or add tip-to-tip.', SF);
end

%% ---------------------------------- Plotting --------------------------------------
f = figure('Color','w');
ax = axes(f); hold(ax,'on'); grid(ax,'on');
title(ax,'Fin Flutter Analysis');
xlabel(ax,'Altitude [ft]'); ylabel(ax,'Velocity [ft/s]');
xlim(ax, [0, alt(end)]);

lp1 = plot(ax, alt, Vf,  'DisplayName','Flutter Velocity', 'LineWidth',2);
lp2 = plot(ax, alt, vel, 'DisplayName','Rocket Velocity',  'LineWidth',2);
xl  = plot(ax, [alt_vmax, alt_vmax], [vel_max, Vf_atvmax], '--k', 'DisplayName', sprintf('Safety Factor %.3f',SF), 'LineWidth',2);
legend(ax,'Location','northeast');

% Annotation arrows in normalized figure coordinates (no tightPosition dependency)
scale     = axis(ax);
pos       = get(ax,'Position');          % [x y w h] in normalized figure units
altN      = (alt_vmax - scale(1)) / (scale(2)-scale(1)) * pos(3) + pos(1);
VfN       = (Vf_atvmax - scale(3)) / (scale(4)-scale(3)) * pos(4) + pos(2);
VelN      = (vel_max   - scale(3)) / (scale(4)-scale(3)) * pos(4) + pos(2);
xAnn      = [altN + 0.06, altN];         % small right offset so arrow points left

annotation('textarrow', xAnn, [VfN  VfN],  'String', sprintf('%.2f ft/s',Vf_atvmax));
annotation('textarrow', xAnn, [VelN VelN], 'String', sprintf('%.2f ft/s',vel_max));
