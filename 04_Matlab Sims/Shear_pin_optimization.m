% Shear Pin Calculator
% Calculates the minimium number and type of shear pins needed
% based the force from drag separation and internal pressure in recovery bays
% Internal pressure is based on research from Off We Go Rockectry's static port hole document
% combined with a mass flow balance
% Drag separation force in caclucated using an equation from MIT
% Uses RASAero II sim data format
% By Zach Pedrick 4/5/2025

close all
clear
clc

%-------------------------------------------------------------------------------------------------------------------------------------------
% User Inputs
% Factor of Safety
FS = 1.5;

% Max recovery G's
G = 20;

% Has two parachutes?
chutes2 = 'no';

% Mass (section division is based on tradional recovery setup with the main by the nosecone and the drogue by the motor)
M1 = 27.09;      % [lbm] Section 1 is forward the main recovery bay separation point
M2 = 6.40;      % [lbm] Section 2 is between the main and drogue recovery bay sepe\aration points
M3 = 31;      % [lbm] Section 3 is aft of the drogue recovery bay separation point

% Coefficient of Drag
CD1 = 0.081;     % CD above drogue recovery bay separation point
CD2 = 0.292;     % CD below drouge recovery bay separation point

% Cross sectional area of the rocket
A1 = 29.9;     % [in] cross sectional area above the recovery bay separation point
A2 = 34.775;     % [in] cross sectional area below the recovery bay separation point

% Recovery Bay Dimensions
D_body = 6.0;        % [in] Inner diameter of the body tube
L_drogue = 11.3215;      % [in] length of drogue recovery bay
L_main = 20.375;        % [in] length of main recovery bay

% Vent Holes
num_vent_drogue = 2;        % number of vent holes in drogue recovery bay
num_vent_main = 2;          % number of vent holes in main recovery bay
D_vent_drogue = 0.125;     % [in] diameter of the vent holes
D_vent_main = 0.125;       % [in] diameter of the vent holes

% Ground conditions
GL = 2920;                  % Altitude of launch site [ft]
temp_ground = 90;           % Temperature at ground level at launch site [f]

file_path = 'C:\Users\james\OneDrive\Documents\GitHub\Rocket-Sims\02_RASAero Sims\Flight Test.CSV';
%-------------------------------------------------------------------------------------------------------------------------------------------

% Miscellaneous variables
R_air = 287;                % [J/(kg*K)] gas constant of air
K_air = 1.4;                % specific heat ratio of air
p_air = (2/(K_air+1))^(K_air/(K_air-1));    % critical pressure ratio for air
screw_name = ["2-56", "M2.5", "M3"];   % list of screws types
screw_shear = [26, 37, 55];        % [lbf] list of shear forces for in order of screw types

% Converting units into metric for calculations
M1 = M1*0.4536;
M2 = M2*0.4536;
M3 = M3*0.4536;
A1 = A1*0.00064516;
A2 = A2*0.00064516;
D_body = D_body*0.0254;
L_drogue = L_drogue*0.0254;
L_main = L_main*0.0254;
D_vent_drogue = D_vent_drogue*0.0254;
D_vent_main = D_vent_main*0.0254;

% Calculating volumes and areas
V_drogue = (pi/4)*(L_drogue*D_body^2);                         % [m^3] Volume of the drogue recovery bay
V_main = (pi/4)*(L_main*D_body^2);                             % [m^3] Volume of the drogue recovery bay
A_vent_drogue = num_vent_drogue*(pi/4)*D_vent_drogue^2;        % [m^2] Total area of vent holes
A_vent_main = num_vent_main*(pi/4)*D_vent_main^2;              % [m^2] Total area of vent holes
A_bulkhead = (pi/4)*(D_body^2);                                % [m^2] Area of the bulkhead

% Importanting data from flight sim
sim_data = readmatrix(file_path);
time = sim_data(:,1);
vel = sim_data(:,18);
vel = vel*0.3048;
AGL = sim_data(:,23);

% Important flight data points
apogee = max(AGL);
apogee_index = find(AGL==apogee,1);

% remove unused data from flight sim arrays. All data after apogee.
AGL(apogee_index+1:length(AGL)) = [];
time(apogee_index+1:length(time)) = [];
vel(apogee_index+1:length(vel)) = [];

% Creating atmospheric data for each time step
altitude = (AGL+GL)/3.281;              % turns altitude above ground level to altitude above sea level [m]
[T, ~, P, rho] = atmosisa(altitude);    % creates temperature, ambient pressure, and density data for each time step of the flight sim.

% Creating arrays for data storage
Mdot_vent = zeros([apogee_index 2]);        % creates a matrix of zeros to store mass flow rate through the vent holes column 1: drouge bay, column 2: main bay, column 3: avionics bay, column 4: motor bay
P_bay_ascent = zeros([apogee_index 2]);     % creates a matrix of zeros to store the absolute pressure [pa] of the recovery bays column 1: drouge bay, column 2: main bay, column 3: avionics bay, column 4: motor bay

% Ascent Simulation
P_bay_ascent(1,:) = P(1);
for i = 1:apogee_index-1  
    Mdot_vent(i,1) = -0.62*A_vent_drogue*sqrt(2*(P_bay_ascent(i,1)^2)/(R_air*T(i))*(K_air/(K_air-1))*((P(i)/P_bay_ascent(i,1))^(2/K_air)-(P(i)/P_bay_ascent(i,1))^((K_air+1)/K_air)));
    Mdot_vent(i,2) = -0.62*A_vent_main*sqrt(2*(P_bay_ascent(i,2)^2)/(R_air*T(i))*(K_air/(K_air-1))*((P(i)/P_bay_ascent(i,2))^(2/K_air)-(P(i)/P_bay_ascent(i,2))^((K_air+1)/K_air)));
    

    P_bay_ascent(i+1,1) = Mdot_vent(i,1)*(time(i+1)-time(i))/V_drogue*R_air*T(i+1)+P_bay_ascent(i,1);
    P_bay_ascent(i+1,2) = Mdot_vent(i,2)*(time(i+1)-time(i))/V_main*R_air*T(i+1)+P_bay_ascent(i,2);
end

Mdot_vent = real(Mdot_vent);
P_bay_ascent = real(P_bay_ascent);


% creating gage pressure data
P_bay_ascent_gauge = (P_bay_ascent-P)/6895;

% Drag Separation Force
V_max = max(vel);
vel_index = find(vel==V_max,1);
P_drogue_ascent = P_bay_ascent_gauge(vel_index:length(P_bay_ascent_gauge),1)*6895;
F_rho = rho(vel_index:length(rho));
F_vel = vel(vel_index:length(vel));

D1 = A1*CD1.*F_rho.*(F_vel.^2)/2;
D2 = A2*CD2.*F_rho.*(F_vel.^2)/2;
N = -(M3*D1 - (M1+M2)*D2)/(M1+M2+M3);

% Combined Separation force from drag and internal pressure
F_sep = (N + P_drogue_ascent*A_bulkhead)*0.2248;
F_drogue = max(F_sep) * FS;
F_main = (M1/0.4536)*G + P_bay_ascent_gauge(length(P_bay_ascent_gauge(:,2)),2)*A_bulkhead*1550;

% Selecting shear pin
number_pins_drogue = 200; % Intializing the number of shear pins for optimization
number_pins_main = 200;   % Intializing the number of shear pins for optimization
if F_drogue <= 0
    shear_pin_drogue = "No Shear Pins Required";
    screw_index_drogue = 1;
    number_pins_drogue = 0;
elseif F_drogue <= 2*screw_shear(1)
    shear_pin_drogue = screw_name(1);
    screw_index_drogue = 1;
    number_pins_drogue = 2;
else
    for i = 1:length(screw_shear)
        temp_number = ceil(F_drogue/screw_shear(i));
        if temp_number ~= 1 && temp_number < number_pins_drogue
            screw_index_drogue = i;
            number_pins_drogue = temp_number;
        end
    end
    shear_pin_drogue = screw_name(screw_index_drogue);
end


    if F_main <= 0
        shear_pin_main = "No Shear Pins Required";
        screw_index_main = 1;
        number_pins_main = 0;
    elseif F_main <= 2*screw_shear(1)
        shear_pin_main = screw_name(1);
        number_pins_main = 2;
    else
        for i = 1:length(screw_shear)
            temp_number = ceil(F_main/screw_shear(i));
            if temp_number ~= 1 && temp_number < number_pins_main
                screw_index_main = i;
                number_pins_main = temp_number;
            end
        end
        shear_pin_main = screw_name(screw_index_main);
    end




% Ejection charge pressure
F_shear_drogue = number_pins_drogue*screw_shear(screw_index_drogue);
EC_pressure_drogue = F_shear_drogue*FS/(A_bulkhead*1550);
F_shear_main = number_pins_main*screw_shear(screw_index_main);
EC_pressure_main = F_shear_main*FS/(A_bulkhead*1550);

% Ejection Charge size
EC_drogue_grams = (EC_pressure_drogue*(V_drogue*61024)*453.593)/(265.92*3307);
EC_main_grams = (EC_pressure_main*(V_main*61024)*453.593)/(265.92*3307);

% Print Shear pins type and number to console
fprintf('----------Parachute Bay----------\n')
fprintf('Max Separation Force %.2f lbs\n',F_drogue)
fprintf('%d %s Required\n',number_pins_drogue,shear_pin_drogue)
fprintf('Shear Force %d lbs\n',F_shear_drogue)
fprintf('Required Ejection Charge Pressure %.2f psi\n',EC_pressure_drogue)
fprintf('Ejection Charge %.2f grams\n',EC_drogue_grams)

fprintf('-----------Payload-----------\n')
fprintf('Max Separation Force %.2f lbs\n',F_main)
fprintf('%d %s Required\n',number_pins_main,shear_pin_main)
fprintf('Shear Force %d lbs\n',F_shear_main)
fprintf('Required Ejection Charge Pressure %.2f psi\n',EC_pressure_main)
fprintf('Ejection Charge %.2f grams\n',EC_main_grams)

% Plotting the data
drogue_ylim = ceil(max(P_bay_ascent_gauge(:,1)));
%main_ylim = ceil(max(P_bay_ascent_gauge(:,2)));
ascent_xlim = time(apogee_index);
f1 = figure;        % ascent data
t1 = tiledlayout(2,1);

nexttile    % Drogue bay during ascent
hold on
grid on
plot(time,P_bay_ascent_gauge(:,1),'LineWidth',2);
title('Drogue Recovery Bay')
ylim([0 drogue_ylim])
xlim([0 ascent_xlim])
ax=gca;
ax.FontSize = 15;

nexttile    % main bay during ascent
hold on
grid on
plot(time,P_bay_ascent_gauge(:,2),'LineWidth',2);
title('Main Recovery Bay')
ylim([0 drogue_ylim])
xlim([0 ascent_xlim])
ax=gca;
ax.FontSize = 15;

title(t1,'\fontsize{15}Rocket Recovery Bay Gauge Pressure vs. Time During Ascent')
xlabel(t1,'\fontsize{15}Time [sec]')
ylabel(t1,'\fontsize{15}Gauge Pressure [psi]')
