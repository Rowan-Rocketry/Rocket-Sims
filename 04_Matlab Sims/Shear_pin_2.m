%% IREC L3: Shear Pins + Vent Holes Sizing Helper (OpenRocket CSV)
% James — this script is designed to be a *repeatable sizing assistant*:
% - Reads an OpenRocket-like CSV (time history)
% - Estimates peak acceleration and peak dynamic pressure (max-Q-ish)
% - Recommends starting shear-pin counts (4-40 nylon) for BOTH joints
% - Recommends vent-hole count + diameter for:
%     (A) Avionics bay (altimeter bay) vents
%     (B) Drogue bay "anti-pressure-lock" vents
%     (C) Main bay "anti-pressure-lock" vents
%
% IMPORTANT:
% - Nylon screw shear strength varies *a lot*. If you can, replace the default
%   pin shear force with a measured value from your own screws/test jig.
% - Vent sizing is a first-pass heuristic; you should validate with ground tests
%   and altimeter data quality checks.
% - This script does NOT compute black powder masses.

clear; clc;

%% ---------------- USER INPUTS ----------------
% Geometry
geom.body_ID_in            = 6.000;      % [in]
geom.body_OD_in            = 6.170;      % [in]
geom.coupler_OD_in         = 5.998;      % [in]
geom.coupler_length_in     = 16.0;       % [in] physical length
geom.coupler_engagement_in = 6.0;        % [in] engaged into each joint

% Mass breakdown (weights in lb; these are your provided values)
% Section convention (typical dual-deploy with avionics between bays):
% - M1: above avionics bay (forward section above the *upper* joint)
% - M2: avionics bay section itself (between joints)
% - M3: below avionics bay (aft section below the *lower* joint)
mass.W_M1_lb = 14.666;  % [lb] above avionics bay coupler (not including av-bay)
mass.W_M2_lb = 3.898;   % [lb] avionics bay coupler/section
mass.W_M3_loaded_lb = 36.336; % [lb] below avionics bay with loaded motor (not used for pin sizing here)
mass.W_M3_burnout_lb = 24.516; % [lb] below avionics bay w/ empty casing (not used for pin sizing here)

% Bay volumes (EMPTY free internal volume)
bays.V_drogue_in3 = 67.9314;  % [in^3]
bays.V_main_in3   = 122.25;   % [in^3]
bays.V_av_in3     = 91.665;   % [in^3]
bays.L_av_in      = 16.0;     % [in]

% Launch conditions (used for density estimate)
atm.launch_alt_ft     = 2920;    % [ft]
atm.launch_temp_F     = 85;      % [F]
atm.launch_press_mbar = 911.03;  % [mbar]
atm.use_ISA_density   = false;   % true = ignore provided T/P and use ISA

% Safety / design settings
design.SF_separation        = 1.75;   % safety factor against drag-separation
design.pin_type             = "4-40 nylon";
design.pin_shear_lbf        = 25;     % [lbf/pin] DEFAULT PLACEHOLDER
% >>> Replace pin_shear_lbf with YOUR measured value. Many 4-40 nylons
%     fall somewhere ~15–35 lbf in single shear depending on screw & fixture.

% Vent-hole heuristic settings
vents.av_prefer_4_holes     = true;   % 4 holes reduce AoA bias
vents.bay_antiLock_num      = 2;      % drogue/main bay vents (anti-pressure-lock) default count
vents.bay_antiLock_d_in     = 0.125;  % [in] default bay vent diameter (1/8")
% You can override vents.bay_antiLock_num and vents.bay_antiLock_d_in if needed.

% CSV input (OpenRocket-like)
csvCfg.enableCSV = true;
csvCfg.filepath  = "FLIGHT DATA.csv";  % <-- change path/name as needed
% If auto-detection fails, set these to exact column names from your CSV:
csvCfg.colTime_s     = "";   % e.g. "Time (s)"
csvCfg.colAlt_m      = "";   % e.g. "Altitude (m)" or "Altitude"
csvCfg.colAlt_ft     = "";   % e.g. "Altitude (ft)"
csvCfg.colVel_mps    = "";   % e.g. "Vertical velocity (m/s)" or "Velocity (m/s)"
csvCfg.colAccel_mps2 = "";   % e.g. "Acceleration (m/s^2)"

%% ---------------- RUN ----------------
results = struct();

% Read flight data (optional)
if csvCfg.enableCSV
    flight = readOpenRocketCSV(csvCfg);
else
    flight = struct();
end

% Compute flight-derived loads if available
if isfield(flight,"t_s") && ~isempty(flight.t_s)
    loads = computeLoadsFromFlight(flight, atm, geom);
else
    loads = struct();
    loads.a_peak_g = 9.16;  % fallback (from your earlier screenshot)
    loads.q_peak_psf = NaN;
    loads.q_peak_kPa = NaN;
    loads.info = "No CSV used; using fallback peak accel only.";
end

% Shear pins: upper & lower joints
pins = recommendShearPins(mass, loads, design);

% Vent holes: avionics bay vents
avVents = recommendAvBayVents(bays, loads, vents);

% Vent holes: drogue & main bay anti-pressure-lock vents
drogueVents = recommendRecoveryBayVents("Drogue", bays.V_drogue_in3, vents);
mainVents   = recommendRecoveryBayVents("Main",   bays.V_main_in3,   vents);

% Print summary
printSummary(loads, pins, avVents, drogueVents, mainVents, design);

%% ---------------- FUNCTIONS ----------------

function flight = readOpenRocketCSV(cfg)
    flight = struct("t_s",[],"alt_m",[],"vel_mps",[],"accel_mps2",[],"raw",[]);
    if ~isfile(cfg.filepath)
        warning("CSV file not found: %s", cfg.filepath);
        return;
    end

    T = readtable(cfg.filepath, 'PreserveVariableNames', true);
    flight.raw = T;

    % Helper to find a column by fuzzy matching
    names = string(T.Properties.VariableNames);

    % Time
    t = tryGetCol(T, names, cfg.colTime_s, ["time","t"]);

    % Altitude: prefer meters; fall back to feet then convert
    alt_m = tryGetCol(T, names, cfg.colAlt_m, ["altitude (m)","altitude_m","altitude"]);
    if isempty(alt_m)
        alt_ft = tryGetCol(T, names, cfg.colAlt_ft, ["altitude (ft)","altitude_ft"]);
        if ~isempty(alt_ft)
            alt_m = alt_ft * 0.3048;
        end
    end

    % Velocity (vertical)
    vel = tryGetCol(T, names, cfg.colVel_mps, ...
        ["vertical velocity","velocity (m/s)","vel (m/s)","speed (m/s)","vz"]);

    % Acceleration (if available)
    accel = tryGetCol(T, names, cfg.colAccel_mps2, ...
        ["acceleration","accel (m/s^2)","az","a (m/s^2)"]);

    % Basic validation
    if isempty(t) || isempty(alt_m) || isempty(vel)
        warning("Could not auto-detect required columns. Set csvCfg.col* to exact headers.");
        disp("Detected headers:");
        disp(names');
        return;
    end

    flight.t_s       = t;
    flight.alt_m     = alt_m;
    flight.vel_mps   = vel;
    flight.accel_mps2 = accel; % may be empty
end

function col = tryGetCol(T, names, exactName, patterns)
    col = [];
    if strlength(string(exactName)) > 0
        if any(names == string(exactName))
            col = T.(exactName);
            col = col(:);
            return;
        else
            warning("Exact column name not found: %s", exactName);
        end
    end

    % Fuzzy search: first match containing any pattern
    lowerNames = lower(names);
    idx = [];
    for p = patterns
        p = lower(string(p));
        hit = find(contains(lowerNames, p), 1, 'first');
        if ~isempty(hit)
            idx = hit; break;
        end
    end

    if ~isempty(idx)
        col = T.(names(idx));
        col = col(:);
    end
end

function loads = computeLoadsFromFlight(flight, atm, geom)
    loads = struct();

    g0 = 9.80665; % m/s^2

    % Acceleration
    if ~isempty(flight.accel_mps2)
        a_mps2 = flight.accel_mps2;
    else
        % Derive accel from velocity
        dt = gradient(flight.t_s);
        a_mps2 = gradient(flight.vel_mps) ./ max(dt, 1e-6);
    end
    loads.a_peak_mps2 = max(a_mps2);
    loads.a_peak_g    = loads.a_peak_mps2 / g0;

    % Density model
    if atm.use_ISA_density
        rho = isaDensityApprox(flight.alt_m);
        loads.rho_model = "ISA approx vs altitude";
    else
        % Use provided ground T/P, and apply a scale-height approximation for pressure vs altitude
        rho = densityFromTPScaleHeight(flight.alt_m, atm.launch_alt_ft, atm.launch_temp_F, atm.launch_press_mbar);
        loads.rho_model = "Ground T/P with scale-height pressure model";
    end

    % Dynamic pressure q = 0.5 * rho * V^2
    q_Pa = 0.5 .* rho .* (flight.vel_mps.^2);
    [q_peak_Pa, idx] = max(q_Pa);
    loads.q_peak_Pa  = q_peak_Pa;
    loads.q_peak_kPa = q_peak_Pa/1000;

    % Convert to psf for intuition
    loads.q_peak_psf = q_peak_Pa * 0.020885434273; % 1 Pa = 0.0208854 psf

    loads.q_peak_time_s = flight.t_s(idx);
    loads.q_peak_alt_m  = flight.alt_m(idx);
    loads.q_peak_vel_mps = flight.vel_mps(idx);

    % Cross-sectional area for reference
    D_od_m = geom.body_OD_in * 0.0254;
    loads.area_ref_m2 = pi*(D_od_m/2)^2;
end

function rho = densityFromTPScaleHeight(alt_m, launch_alt_ft, launch_temp_F, launch_press_mbar)
    % Approximate: assume isothermal layer at ground temp and exponential pressure decay.
    % This is good enough for q-estimation and vent heuristics.

    R = 287.05; % J/kg-K
    T = (launch_temp_F - 32)*5/9 + 273.15; % K
    p0 = launch_press_mbar * 100; % Pa
    z0 = launch_alt_ft * 0.3048; % m
    z = alt_m;

    g0 = 9.80665;
    H = R*T/g0; % scale height [m] for isothermal atmosphere
    p = p0 .* exp(-(z - z0)./H);
    rho = p./(R*T);
end

function rho = isaDensityApprox(alt_m)
    % Simple ISA troposphere approx up to 11 km
    % Not perfect but fine for q estimates.
    T0 = 288.15; p0 = 101325; g = 9.80665; R = 287.05; L = 0.0065;
    z = max(0, alt_m);
    T = T0 - L*z;
    p = p0 .* (T./T0).^(g/(R*L));
    rho = p./(R*T);
end

function pins = recommendShearPins(mass, loads, design)
    % Inertial separation loads at peak acceleration (conservative)
    % Upper joint sees mass above it: M1
    % Lower joint sees mass above it: M1 + M2
    a_ftps2 = loads.a_peak_g * 32.174;

    F_upper_lbf = mass.W_M1_lb * (a_ftps2/32.174);           % W*(a/g)
    F_lower_lbf = (mass.W_M1_lb + mass.W_M2_lb) * (a_ftps2/32.174);

    % Add a conservative aero margin factor if q is available (optional)
    % NOTE: We do NOT convert q directly into joint-opening force; we apply a small multiplier.
    if isfield(loads,'q_peak_psf') && ~isnan(loads.q_peak_psf)
        aeroFactor = 1.15; % mild margin bump for high-Q rockets
    else
        aeroFactor = 1.00;
    end

    F_upper_req = F_upper_lbf * design.SF_separation * aeroFactor;
    F_lower_req = F_lower_lbf * design.SF_separation * aeroFactor;

    pins.F_upper_inertial_lbf = F_upper_lbf;
    pins.F_lower_inertial_lbf = F_lower_lbf;
    pins.F_upper_req_lbf = F_upper_req;
    pins.F_lower_req_lbf = F_lower_req;

    pins.pin_shear_lbf = design.pin_shear_lbf;

    pins.N_upper = max(1, ceil(F_upper_req / design.pin_shear_lbf));
    pins.N_lower = max(1, ceil(F_lower_req / design.pin_shear_lbf));

    % Practical constraints for symmetry/repeatability
    % Prefer 3 over 2 on larger rockets; cap at a reasonable number.
    pins.N_upper = enforcePracticalPinCount(pins.N_upper);
    pins.N_lower = enforcePracticalPinCount(pins.N_lower);
end

function N = enforcePracticalPinCount(N)
    % Encourage 3 pins minimum once you get above 2, for symmetry.
    if N == 2
        N = 3;
    end
    % Cap to avoid "charge escalation by pin escalation"
    if N > 6
        N = 6;
    end
end

function av = recommendAvBayVents(bays, loads, vents)
    % Heuristic based on bay volume and high-speed profile:
    % - 4 holes if preferred, otherwise 3
    % - Diameter chosen from a coarse lookup using volume and "fast ascent" hint.

    V = bays.V_av_in3;

    % Fast-ascent hint: if peak velocity > ~200 m/s or qpeak > ~800 psf, treat as fast
    fast = false;
    if isfield(loads,'q_peak_psf') && ~isnan(loads.q_peak_psf) && loads.q_peak_psf > 800
        fast = true;
    end
    if isfield(loads,'q_peak_vel_mps') && ~isempty(loads.q_peak_vel_mps) && loads.q_peak_vel_mps > 200
        fast = true;
    end

    if vents.av_prefer_4_holes
        av.N = 4;
    else
        av.N = 3;
    end

    % Diameter selection (starting points, 6" class typical ranges)
    % For ~90 in^3 at high speed: 4 × 3/16" is a good baseline.
    if V < 60
        d = 0.15625; % 5/32
    elseif V < 120
        d = 0.1875;  % 3/16
    else
        d = 0.25;    % 1/4
    end

    % If fast ascent, bias one step larger, but not beyond 1/4
    if fast
        if d < 0.1875
            d = 0.1875;
        elseif d < 0.25 && V > 120
            d = 0.25;
        end
    end

    av.d_in = d;
    av.fast_ascent = fast;
    av.V_in3 = V;
    av.L_in = bays.L_av_in;
end

function bay = recommendRecoveryBayVents(name, V_in3, vents)
    % These are NOT altimeter vents—these are "anti-pressure-lock" vents.
    % You generally want modest venting so the bay equalizes and doesn't fight separation.

    bay.name = name;
    bay.V_in3 = V_in3;

    % Simple heuristic: small bays can use 1 hole; larger bays 2 holes.
    if V_in3 < 80
        bay.N = max(1, round(vents.bay_antiLock_num/2));
    else
        bay.N = vents.bay_antiLock_num;
    end

    bay.d_in = vents.bay_antiLock_d_in;
end

function printSummary(loads, pins, avVents, drogueVents, mainVents, design)
    fprintf("\n==================== RESULTS ====================\n");

    fprintf("\n--- Flight-derived loads ---\n");
    fprintf("Peak acceleration: %.2f g\n", loads.a_peak_g);

    if isfield(loads,'q_peak_psf') && ~isnan(loads.q_peak_psf)
        fprintf("Peak dynamic pressure (q): %.0f psf  (%.1f kPa)\n", loads.q_peak_psf, loads.q_peak_kPa);
        fprintf("  at t=%.2f s, alt=%.0f m, v=%.1f m/s\n", loads.q_peak_time_s, loads.q_peak_alt_m, loads.q_peak_vel_mps);
        fprintf("Density model: %s\n", loads.rho_model);
    else
        fprintf("Peak dynamic pressure (q): (not computed — no CSV)\n");
    end

    fprintf("\n--- Shear pins (%s) ---\n", design.pin_type);
    fprintf("Assumed shear strength per pin: %.1f lbf  (UPDATE THIS WITH YOUR TEST DATA)\n", design.pin_shear_lbf);
    fprintf("Safety factor: %.2f\n", design.SF_separation);

    fprintf("Upper joint inertial load estimate: %.0f lbf\n", pins.F_upper_inertial_lbf);
    fprintf("Lower joint inertial load estimate: %.0f lbf\n", pins.F_lower_inertial_lbf);
    fprintf("Upper joint required (SF, margin): %.0f lbf -> Recommend %d pins\n", pins.F_upper_req_lbf, pins.N_upper);
    fprintf("Lower joint required (SF, margin): %.0f lbf -> Recommend %d pins\n", pins.F_lower_req_lbf, pins.N_lower);

    fprintf("\n--- Vent holes ---\n");
    fprintf("Avionics bay vents (altimeter bay): %d holes @ %.3f in diameter\n", avVents.N, avVents.d_in);
    fprintf("  (Bay V=%.1f in^3, L=%.1f in, fast-ascent=%d)\n", avVents.V_in3, avVents.L_in, avVents.fast_ascent);

    fprintf("%s bay anti-pressure-lock vents: %d holes @ %.3f in diameter\n", drogueVents.name, drogueVents.N, drogueVents.d_in);
    fprintf("%s bay anti-pressure-lock vents: %d holes @ %.3f in diameter\n", mainVents.name,   mainVents.N,   mainVents.d_in);

    fprintf("\nNotes:\n");
    fprintf("- Pin counts are a conservative starting point. Validate with full-stack ground tests.\n");
    fprintf("- Av-bay vent recommendation is a starting pattern; place holes on a clean cylindrical section.\n");
    fprintf("- This script intentionally does not compute black powder quantities.\n");
    fprintf("=================================================\n\n");
end
