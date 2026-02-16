clear; close all; clc;

%% ---------------- USER SETTINGS ----------------
csvPath = "OpenRocket flight data.csv";

% Desired x-axis max for ALL plots
xMaxPlot = 230;

% y-limits (leave [] for auto)
altYLim   = [];
velYLim   = [];
accYLim_G = [];

% Event detection tuning (works well for OpenRocket-style traces)
ignitionAccThresh_G = 0.30;
burnoutAccThresh_G  = 0.05;
burnoutHold_s       = 0.30;
mainSearchDelay_s   = 5.0;

% OPTIONAL overrides (sec). Leave NaN to auto-detect.
tIgnition_override = NaN;
tBurnout_override  = NaN;
tApogee_override   = NaN;
tMain_override     = NaN;

% Plot styling
lw = 2.0;

% Curve colors (requested)
cAlt = [0.0000 0.4470 0.7410]; % blue
cVel = [0.8500 0.0000 0.0000]; % red
cAcc = [1.0000 0.6000 0.0000]; % orange

% Event styling (+ ground hit)
eventNames  = ["Motor ignition","Motor burnout","Apogee","Main parachute deployment","Ground hit"];
eventColors = [...
    0.55 0.55 0.55;  % ignition (gray)
    0.30 0.30 0.30;  % burnout (dark gray)
    0.25 0.45 0.85;  % apogee (blue-ish)
    0.55 0.25 0.75;  % main (purple-ish)
    0.15 0.65 0.20]; % ground hit (green)
eventLineStyle = '--';
eventLW = 1.2;
labelFontSize = 9;

% Label vertical placement: fraction of y-span down from top
labelInsetFrac = 0.23;   % your preferred placement (0.23 => 77% of plot height)

% Export
doExport = true;
exportDPI = 300;
outAlt = "OpenRocket_Altitude.png";
outVel = "OpenRocket_Velocity.png";
outAcc = "OpenRocket_Acceleration.png";

%% ---------------- READ CSV (skip '#') ----------------
% OpenRocket exports often include comment lines starting with '#'
opts = detectImportOptions(csvPath, "NumHeaderLines", 0);
opts.VariableNamesLine = 0;
opts.DataLines = [1 Inf];
opts.CommentStyle = '#';
opts.Delimiter = ',';
opts.ExtraColumnsRule = 'ignore';
opts.EmptyLineRule = 'read';
opts = setvartype(opts, 'double');

T = readtable(csvPath, opts);

% Fallback: sometimes everything ends up in one column
if width(T) == 1
    raw = string(T{:,1});
    parts = split(raw, ",");
    M = str2double(parts);
    goodRow = all(isfinite(M), 2);
    M = M(goodRow, :);

    t    = M(:,1);
    alt  = M(:,2);   % already ft in your export
    vel  = M(:,3);   % m/s
    accG = M(:,4);   % already G in your export
else
    t    = T{:,1};
    alt  = T{:,2};
    vel  = T{:,3};
    accG = T{:,4};
end

good = isfinite(t) & isfinite(alt) & isfinite(vel) & isfinite(accG);
t = t(good); alt = alt(good); vel = vel(good); accG = accG(good);

% Labels (match OpenRocket tooltips / your report)
altLabel = "Altitude (ft)";
velLabel = "Vertical velocity (m/s)";
accLabel = "Vertical acceleration (G)";

%% ---------------- EVENTS (consistent times across all plots) ----------------
[tIgn, tBo, tApo, tMain] = detectEvents_OpenRocketUnits( ...
    t, alt, vel, accG, ...
    ignitionAccThresh_G, burnoutAccThresh_G, burnoutHold_s, mainSearchDelay_s);

% Apply overrides if provided
if isfinite(tIgnition_override), tIgn = tIgnition_override; end
if isfinite(tBurnout_override),  tBo  = tBurnout_override;  end
if isfinite(tApogee_override),   tApo = tApogee_override;   end
if isfinite(tMain_override),     tMain = tMain_override;    end

tGround = t(end);
eventTimes = [tIgn, tBo, tApo, tMain, tGround];

%% ---------------- FIGURE 1: ALTITUDE ----------------
fig1 = figure('Name','Altitude vs Time','Color','w');
plot(t, alt, 'LineWidth', lw, 'Color', cAlt); grid on;
xlabel('Time (s)'); ylabel(altLabel); title('Altitude vs Time');
xlim([0 xMaxPlot]);
if ~isempty(altYLim), ylim(altYLim); end

% Disable scientific notation on altitude axis
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%,.0f');  % comment out if you don't want commas

addEventMarkers_TopAnchored(gca, eventTimes, eventNames, eventColors, eventLineStyle, eventLW, labelFontSize, labelInsetFrac);

%% ---------------- FIGURE 2: VELOCITY ----------------
fig2 = figure('Name','Vertical velocity vs Time','Color','w');
plot(t, vel, 'LineWidth', lw, 'Color', cVel); grid on;
xlabel('Time (s)'); ylabel(velLabel); title('Vertical velocity vs Time');
xlim([0 xMaxPlot]);
if ~isempty(velYLim), ylim(velYLim); end

addEventMarkers_TopAnchored(gca, eventTimes, eventNames, eventColors, eventLineStyle, eventLW, labelFontSize, labelInsetFrac);

%% ---------------- FIGURE 3: ACCELERATION ----------------
fig3 = figure('Name','Vertical acceleration vs Time','Color','w');
plot(t, accG, 'LineWidth', lw, 'Color', cAcc); grid on;
xlabel('Time (s)'); ylabel(accLabel); title('Vertical acceleration vs Time');
xlim([0 xMaxPlot]);
if ~isempty(accYLim_G), ylim(accYLim_G); end

addEventMarkers_TopAnchored(gca, eventTimes, eventNames, eventColors, eventLineStyle, eventLW, labelFontSize, labelInsetFrac);

%% ---------------- EXPORT ----------------
if doExport
    exportgraphics(fig1, outAlt, 'Resolution', exportDPI);
    exportgraphics(fig2, outVel, 'Resolution', exportDPI);
    exportgraphics(fig3, outAcc, 'Resolution', exportDPI);
end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function [tIgn, tBo, tApo, tMain] = detectEvents_OpenRocketUnits(t, alt_ft, vel_ms, acc_G, ...
    ignThreshG, boThreshG, boHold_s, mainDelay_s)

    % Apogee: altitude peak (single definition)
    [~, iApo] = max(alt_ft);
    tApo = t(iApo);

    dt = median(diff(t));
    if ~isfinite(dt) || dt <= 0
        dt = 0.01;
    end

    % Ignition: first time accel exceeds threshold (G)
    iIgn = find(acc_G > ignThreshG, 1, 'first');
    if isempty(iIgn), iIgn = 1; end
    tIgn = t(iIgn);

    % Burnout: sustained accel below threshold (G)
    w = max(3, round(0.10 / dt));
    accSm = movmean(acc_G, w);
    holdN = max(1, round(boHold_s / dt));

    iBo = NaN;
    for k = iIgn:length(t)-holdN
        if all(accSm(k:k+holdN-1) < boThreshG)
            iBo = k;
            break;
        end
    end
    if ~isfinite(iBo)
        [~, iBo] = max(vel_ms); % fallback
    end
    tBo = t(iBo);

    % Main deploy: largest +acc spike after apogee + delay (G)
    i0 = find(t >= (tApo + mainDelay_s), 1, 'first');
    if isempty(i0), i0 = iApo; end
    iEnd = max(i0, length(t)-5);
    if iEnd <= i0, iEnd = length(t); end

    [~, iRel] = max(acc_G(i0:iEnd));
    iMain = i0 + iRel - 1;
    tMain = t(iMain);
end

function addEventMarkers_TopAnchored(ax, eventTimes, eventNames, eventColors, ls, lw, fs, insetFrac)
    % Event labels anchored to top edge of axes (data coords) and hang downward.

    hold(ax,'on'); drawnow;

    xl = xlim(ax);
    yl = ylim(ax);

    yTop = yl(2);
    yInset = insetFrac * (yl(2)-yl(1));
    yText = yTop - yInset;

    % Small x offset so label isn't exactly on the line
    xOff = 0.005 * (xl(2)-xl(1));

    for i = 1:numel(eventTimes)
        tE = double(eventTimes(i));
        if ~isfinite(tE), continue; end

        xline(ax, tE, ls, 'Color', eventColors(i,:), 'LineWidth', lw);

        text(ax, tE + xOff, yText, eventNames(i), ...
            'Units','data', ...
            'Rotation', 90, ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','top', ...
            'Color', eventColors(i,:), ...
            'FontSize', fs, ...
            'Clipping','on');
    end
end