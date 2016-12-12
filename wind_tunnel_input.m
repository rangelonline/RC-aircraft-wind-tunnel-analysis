%{
This input deck streamlines the wind tunnel data collection and
post-processing associated with wind tunnel analysis. While this script
post-processes data collected on an RC aircraft, the input deck outline can
easily be modified or expanded to examine the aerodynamic characteristics 
of a wide variety of aircraft and structures. 

The present data was collected at the Oran W. Nicks Low Speed Wind Tunnel 
in College Station, Texas as part of an aerospace engineering senior 
design capstone project at Texas A&M University in College Station, TX. 

wind tunnel:            Oran W. Nicks LSWT
                        1775 George Bush Dr W.
                        College Station, TX 77845
test date:              November 15, 2016
company:                ATLAS
aircraft:               AG-16 Dragonfly (1:6 scale RC version)

CALLED FUNCTIONS:
1)                      data_files(): 
2)                      read_WT_data(filename, filetype, sref, cref)
3)                      moment_transfer(mystruct, field, point1, point2, norm)
4)                      aerodynamic_center(struct, xb, cref, cutoff)

By John Rangel 12/12/16
%}

close all
clc
clear

%--------------------------------------------------------------------------
%                           geometries
%--------------------------------------------------------------------------

% aircraft geometry:
Sref = 4.125; %reference area [ft^2]
cref = 0.75; %reference chord [ft]
b = 5.5; %wingspan [ft]
AR = (b^2)/Sref;

% setup geometry:
xb = 17.625/12.0; % nose to WT balance [ft]
xcg = 19.5/12.0; % nose to CG [ft]

% deg to rad:
deg2rad = pi/180.0;

%--------------------------------------------------------------------------
%                           read wind tunnel data
%--------------------------------------------------------------------------

% specify run files:
WT_data_files = data_files();

% aerodynamic and longitudinal behavior testing:
dynamic_pressure_sweep = read_WT_data(WT_data_files.Qsweep, 'stat', Sref, cref);
alpha_sweep = read_WT_data(WT_data_files.alpha_sweep, 'stat', Sref, cref);

% directional stability testing:
beta_sweep = read_WT_data(WT_data_files.beta_sweep, 'stat', Sref, cref);

% control surface effectiveness testing:
takeoff_elevator_sweep = read_WT_data(WT_data_files.takeoff_elevator_sweep, 'stat', Sref, cref);
cruise_aileron_sweep = read_WT_data(WT_data_files.cruise_aileron_sweep, 'stat', Sref, cref);
cruise_elevator_sweep = read_WT_data(WT_data_files.cruise_elevator_sweep, 'stat', Sref, cref);
cruise_rudder_sweep = read_WT_data(WT_data_files.cruise_rudder_sweep, 'stat', Sref, cref);

% engine 
takeoff_Emax_elevator = read_WT_data(WT_data_files.takeoff_Emax_elevator, 'stat', Sref, cref);
engine_sweep_with_elevator = read_WT_data(WT_data_files.cruise_Esweep_ElevatorSweep,'stat', Sref, cref);
engine_sweep_with_rudder = read_WT_data(WT_data_files.cruise_Esweep_RudderSweep, 'stat', Sref, cref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
%                        lift, drag analysis
%--------------------------------------------------------------------------
% NOTE: Check the CL vs alpha curve to spot where trend is no longer
%       linear. This is the alpha where seperation begins. All data past 
%       this point should not be used in calculations. 

% separation cutoff:
separation_points = 4;
cutoff = length(alpha_sweep.CL)-separation_points;

% convert deg to rad:
%alpha_sweep.alpha = alpha_sweep.alpha*deg2rad;

% calculate lift curve slope:
CL_vs_alpha_coefficients = polyfit(alpha_sweep.alpha(1:cutoff), alpha_sweep.CL(1:cutoff), 1);
linear_lift_curve = polyval(CL_vs_alpha_coefficients, alpha_sweep.alpha);
linear_lift_curve_slope = CL_vs_alpha_coefficients(1);

% lift characteristics:
zero_lift_alpha_poly = polyfit(alpha_sweep.CL, alpha_sweep.alpha, 1);
zero_lift_alpha = polyval(zero_lift_alpha_poly, 0.0);
CL_max = max(alpha_sweep.CL);
critical_alpha = max(alpha_sweep.alpha);

% print lift characteristics to console:
fprintf(['LIFT CHARACTERISTCS:\n',...
    'CL vs alpha slope: %.3f per deg \n', ...
    'CL max: %.3f \n', ...
    'zero-lift alpha: %.3f deg \n', ...
    'critical alpha: %.3f deg \n'], linear_lift_curve_slope, CL_max, ...
                                    zero_lift_alpha, critical_alpha);

%--------------------------------------------------------------------------
%                        drag polar, L/D
%--------------------------------------------------------------------------

% drag polar calculations:
% NOTE: assume quadratic drag polar
drag_polar_coefficients = polyfit(alpha_sweep.CL(1:cutoff), alpha_sweep.CD(1:cutoff), 2);
CL_range = linspace(min(alpha_sweep.CL), max(alpha_sweep.CL), 40);
quadratic_drag_polar = polyval(drag_polar_coefficients, CL_range);
CD_zero_lift = polyval(drag_polar_coefficients, 0.0);
[CD_min, I] = min(quadratic_drag_polar);
CL_min_drag = CL_range(I);

% calculate Oswald efficiency factor using area minimization:
e_range = linspace(0.1, 0.99, 20);
area_diff = zeros(length(e_range), 1);
xmin = min(alpha_sweep.CL);
xmax = max(alpha_sweep.CL);
count = 1;
for e_i = e_range
    drag_polar_1 = @(x) drag_polar_coefficients(1)*x.^2 + ...
                    drag_polar_coefficients(2)*x + drag_polar_coefficients(3);
    drag_polar_2 = @(x) CD_min + (1/(pi*e_i*AR))*(x - CL_min_drag).^2;
    area_1 = integral(drag_polar_1, xmin, xmax);
    area_2 = integral(drag_polar_2, xmin, xmax);
    area_diff(count) = abs(area_1 - area_2);
    count = count+1;
end
[~, e_index] = min(area_diff);
e_factor = e_range(e_index);
drag_polar = @(x, e) CD_min + (1/(pi*e*AR))*(x - CL_min_drag).^2;

% calculate L/D max:
alpha_sweep.LoD = alpha_sweep.CL./alpha_sweep.CD;
[LoD_max, LoD_max_index] = max(alpha_sweep.LoD);
LoD_max_alpha = alpha_sweep.alpha(LoD_max_index);

% print drag polar, LoD characteristics to console:
fprintf(['\nDRAG POLAR, LoD CHARACTERISTICS: \n', ...
        'CD zero-lift: %.3f \n', ...
        'CDmin: %.3f \n', ...
        'CL min-drag: %.3f \n',...
        'e: %.3f \n', ...
        'LoD max: %.3f \n',...
        'LoD max alpha: %.3f deg \n\n'], CD_zero_lift, CD_min, CL_min_drag, e_factor, ...
                                LoD_max, LoD_max_alpha);

%--------------------------------------------------------------------------
%                       longitudinal stability
%--------------------------------------------------------------------------

% move measured moments from balance center to CG:
alpha_sweep.Cm_cg = moment_transfer(alpha_sweep, 'Cm', xb, xcg, cref);

% Cm_alpha:
Cm_vs_alpha_coefficients = polyfit(alpha_sweep.alpha, alpha_sweep.Cm_cg, 1);
Cm_alpha = Cm_vs_alpha_coefficients(1);
Cm_vs_alpha_fit = polyval(Cm_vs_alpha_coefficients, alpha_sweep.alpha);

% aerodynamic center:
xac_bar = aerodynamic_center(alpha_sweep, xb, cref, cutoff);
xac = xac_bar*cref;

% static margin:
static_margin = 100.0*(xac - xcg)/cref;

% print results to console:
fprintf(['LONGITUDINAL STABILITY: \n', ...
        'Cm_alpha: %.3f per deg \n', ...
        'xac: %.3f ft \n', ...
        'static margin: %.1f percent \n\n'], Cm_alpha, xac, static_margin)

%--------------------------------------------------------------------------
%                       directional stability
%--------------------------------------------------------------------------

% move measured moments from balance center to aircraft CG:
beta_sweep.Cl_cg = beta_sweep.Cl;
beta_sweep.Cn_cg = moment_transfer(beta_sweep, 'Cn', xb, xcg, b);

% weathercock stability (Cn_beta):
Cn_vs_beta_coefficients = polyfit(beta_sweep.beta, beta_sweep.Cn_cg, 1);
Cn_beta = Cn_vs_beta_coefficients(1);
Cn_vs_beta_fit = polyval(Cn_vs_beta_coefficients, beta_sweep.beta);

% lateral stability (Cl_beta):
Cl_vs_beta_coefficients = polyfit(beta_sweep.beta, beta_sweep.Cl_cg, 1);
Cl_beta = Cl_vs_beta_coefficients(1);
Cl_vs_beta_fit = polyval(Cl_vs_beta_coefficients, beta_sweep.beta);

% print results to console:
fprintf(['DIRECTIONAL STABILITY: \n', ...
        'weathercock stability (Cn beta): %.3f per deg \n', ...
        'lateral stability (Cl beta): %.3f per deg \n\n'], Cn_beta, Cl_beta)

%--------------------------------------------------------------------------
%                     control surface effectiveness
%--------------------------------------------------------------------------

% move measured moments from balance center to aircraft CG:
cruise_elevator_sweep.Cm_cg = moment_transfer(cruise_elevator_sweep, 'Cm', xb, xcg, cref);
cruise_rudder_sweep.Cn_cg = moment_transfer(cruise_rudder_sweep, 'Cn', xb, xcg, b);
cruise_aileron_sweep.Cl_cg = cruise_aileron_sweep.Cl;

% add angle measurement deflections:
% NOTE: negative angles represent aileron, rudder left deflection
% NOTE: positive angles represent elevator down deflection
cruise_elevator_sweep.elevator_deflection = [0.0 20.0 10.0 0.0 -5.0 -10.0 0.0];
cruise_aileron_sweep.aileron_deflection = [0.0 -16.5 -8.25 0.0 8.5 17.0 0.0];
cruise_rudder_sweep.rudder_deflection = [0.0 -11.2 -5.6 0.0 6.0 12.0 0.0];

% elevator effectiveness:
% NOTE: negative Cm --> counterclockwise rotation
Cm_vs_ElevatorDeflection_coefficients = polyfit(cruise_elevator_sweep.elevator_deflection, ...
                            cruise_elevator_sweep.Cm_cg, 1);
elevator_effectiveness = Cm_vs_ElevatorDeflection_coefficients(1);
Cm_vs_ElevatorDeflection_fit = polyval(Cm_vs_ElevatorDeflection_coefficients, ...
                                cruise_elevator_sweep.elevator_deflection);

% rudder effectiveness:
Cn_vs_RudderDeflection_coefficients = polyfit(cruise_rudder_sweep.rudder_deflection, ...
                                        cruise_rudder_sweep.Cn_cg, 1);
rudder_effectiveness = Cn_vs_RudderDeflection_coefficients(1);
Cn_vs_RudderDeflection_fit = polyval(Cn_vs_RudderDeflection_coefficients, cruise_rudder_sweep.rudder_deflection);

% aileron effectiveness:
Cl_vs_AileronDeflection_coefficients = polyfit(cruise_aileron_sweep.aileron_deflection, ...
                                        cruise_aileron_sweep.Cl_cg, 1);
aileron_effectiveness = Cl_vs_AileronDeflection_coefficients(1);
Cl_vs_AileronDeflection_fit = polyval(Cl_vs_AileronDeflection_coefficients, cruise_aileron_sweep.aileron_deflection);

% print results to console:
fprintf(['CONTROL SURFACE EFFECTIVENESS: \n', ...
        'aileron effectiveness: %.3f per deg \n', ...
        'elevator effectiveness: %.3f per deg \n', ...
        'rudder effectiveness: %.3f per deg \n\n'], ...
        aileron_effectiveness, elevator_effectiveness, rudder_effectiveness)

%--------------------------------------------------------------------------
%                          propellor influence
%--------------------------------------------------------------------------

% move measured moments from balance center to aircraft CG:
engine_sweep_with_elevator.Cm_cg = moment_transfer(engine_sweep_with_elevator, 'Cm', xb, xcg, cref);
engine_sweep_with_rudder.Cn_cg = moment_transfer(engine_sweep_with_rudder, 'Cn', xb, xcg, b);

% throttle percentage:
throttle_percentage = [0.0 50.0 100.0];

% organize elevator data:
MaxElevator_no_throttle = max(cruise_elevator_sweep.Cm_cg);
MaxElevator_half_throttle = engine_sweep_with_elevator.Cm_cg(2);
MaxElevator_full_throttle = engine_sweep_with_elevator.Cm_cg(3);
MaxElevator_throttle_sweep = [MaxElevator_no_throttle MaxElevator_half_throttle MaxElevator_full_throttle];
MaxElevator_Cm_loss = MaxElevator_throttle_sweep(1) - MaxElevator_throttle_sweep(3);

MinElevator_no_throttle = min(cruise_elevator_sweep.Cm_cg);
MinElevator_half_throttle = engine_sweep_with_elevator.Cm_cg(4);
MinElevator_full_throttle = engine_sweep_with_elevator.Cm_cg(5);
MinElevator_throttle_sweep = [MinElevator_no_throttle MinElevator_half_throttle MinElevator_full_throttle];
MinElevator_Cm_loss = MinElevator_throttle_sweep(3) - MinElevator_throttle_sweep(1);

% organize rudder data:
MaxRudder_no_throttle = max(cruise_rudder_sweep.Cn_cg);
MaxRudder_half_throttle = engine_sweep_with_rudder.Cn_cg(4);
MaxRudder_full_throttle = engine_sweep_with_rudder.Cn_cg(5);
MaxRudder_throttle_sweep = [MaxRudder_no_throttle MaxRudder_half_throttle MaxRudder_full_throttle];
MaxRudder_Cn_loss = MaxRudder_throttle_sweep(1) - MaxRudder_throttle_sweep(3);

MinRudder_no_throttle = min(cruise_rudder_sweep.Cn_cg);
MinRudder_half_throttle = engine_sweep_with_rudder.Cn_cg(2);
MinRudder_full_throttle = engine_sweep_with_rudder.Cn_cg(3);
MinRudder_throttle_sweep = [MinRudder_no_throttle MinRudder_half_throttle MinRudder_full_throttle];
MinRudder_Cn_loss = MinRudder_throttle_sweep(1) - MinRudder_throttle_sweep(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
%                       plot lift, drag analysis
%--------------------------------------------------------------------------

% plot CL vs alpha:
figure()
subplot(1,2,1)
hold on
scatter(alpha_sweep.alpha, alpha_sweep.CL)
plot(alpha_sweep.alpha, linear_lift_curve)
grid on
title('Aircraft CL vs alpha curve')
xlabel('angle of attack [deg]')
ylabel('CL')
legend_text = sprintf('Aircraft lift curve slope: %.3f per deg', linear_lift_curve_slope);
legend('WT data', legend_text)

% plot CD vs alpha:
subplot(1,2,2)
scatter(alpha_sweep.alpha, alpha_sweep.CD)
grid on
title('CD vs alpha curve')
xlabel('angle of attack [deg]')
ylabel('CD')
legend('WT data')

%--------------------------------------------------------------------------
%                       plot drag polars
%--------------------------------------------------------------------------

% plot of converging e through area minimization technique:
offset = 3;
e1 = e_range(e_index+offset);
e2 = e_range(e_index-offset);
figure()
hold on
scatter(alpha_sweep.CL, alpha_sweep.CD)
fplot(@(x) drag_polar(x, e2), [min(alpha_sweep.CL) max(alpha_sweep.CL)], '-.r')
fplot(@(x) drag_polar(x, e_factor), [min(alpha_sweep.CL) max(alpha_sweep.CL)], 'b')
fplot(@(x) drag_polar(x, e1), [min(alpha_sweep.CL) max(alpha_sweep.CL)], '--r')
e1_legend = sprintf('e=%.3f', e1);
e2_legend = sprintf('e=%.3f', e2);
e_factor_legend = sprintf('e=%.3f', e_factor);
grid on
legend('WT test data', e2_legend, e_factor_legend, e1_legend);
xlabel('CL')
ylabel('CD')
title('Drag polar: e estimation')

% plot L/D
figure()
subplot(1,2,1)
scatter(alpha_sweep.alpha, alpha_sweep.LoD)
title('lift over drag')
xlabel('angle of attack [deg]')
ylabel('L/D')
legend('WT data')
grid on

% plot drag polar:
subplot(1,2,2)
hold on
scatter(alpha_sweep.CL, alpha_sweep.CD)
plot(CL_range, quadratic_drag_polar)
grid on
title('Drag polar')
xlabel('CL')
ylabel('CD')
legend('WT data', 'drag polar: quadratic fit')

%--------------------------------------------------------------------------
%                       plot stability
%--------------------------------------------------------------------------

% create subplot figure:
figure()

% plot longitudinal stability (Cm alpha):
subplot(1,3,1)
hold on
scatter(alpha_sweep.alpha, alpha_sweep.Cm_cg)
plot(alpha_sweep.alpha, Cm_vs_alpha_fit)
grid on
title('Longitudinal stability')
xlabel('alpha [deg]')
ylabel('Cm about aircraft CG')
Cm_alpha_legend = sprintf('Cm alpha: %.3f per deg', Cm_alpha);
legend('WT data', Cm_alpha_legend)

% plot weathercock stability (Cn beta):
subplot(1,3,2)
hold on
scatter(beta_sweep.beta, beta_sweep.Cn_cg)
plot(beta_sweep.beta, Cn_vs_beta_fit)
grid on
title('Weathercock stability')
xlabel('beta [deg]')
ylabel('Cn about aircraft CG')
Cn_beta_legend = sprintf('Cn beta: %.3f per deg', Cn_beta);
legend('WT data', Cn_beta_legend)

% plot lateral stability (Cl beta):
subplot(1,3,3)
hold on
scatter(beta_sweep.beta, beta_sweep.Cl_cg)
plot(beta_sweep.beta, Cl_vs_beta_fit)
grid on
title('Lateral stability')
xlabel('beta [deg]')
ylabel('Cl about CG')
Cl_beta_legend = sprintf('Cl beta: %.3f per deg', Cl_beta);
legend('WT data', Cl_beta_legend)

%--------------------------------------------------------------------------
%                    plot control surface effectiveness
%--------------------------------------------------------------------------

% create subplot figure:
figure()

% plot elevator effectiveness:
subplot(1,3,1)
hold on
scatter(cruise_elevator_sweep.elevator_deflection, cruise_elevator_sweep.Cm_cg)
plot(cruise_elevator_sweep.elevator_deflection, Cm_vs_ElevatorDeflection_fit)
title('elevator effectiveness')
xlabel('elevator deflection [deg]')
ylabel('Cm about CG')
grid on
Cm_ElevatorDeflection_legend = sprintf('elevator effectiveness: %.3f per deg', elevator_effectiveness);
legend('WT test data', Cm_ElevatorDeflection_legend)

% plot rudder effectiveness:
subplot(1,3,2)
hold on
scatter(cruise_rudder_sweep.rudder_deflection, cruise_rudder_sweep.Cn_cg)
plot(cruise_rudder_sweep.rudder_deflection, Cn_vs_RudderDeflection_fit)
title('rudder effectiveness')
xlabel('rudder deflection [deg]')
ylabel('Cn about CG')
grid on
Cn_RudderDeflection_legend = sprintf('rudder effectiveness: %.3f per deg', rudder_effectiveness);
legend('WT test data', Cn_RudderDeflection_legend)

% plot aileron effectiveness:
subplot(1,3,3)
hold on
scatter(cruise_aileron_sweep.aileron_deflection, cruise_aileron_sweep.Cl_cg)
plot(cruise_aileron_sweep.aileron_deflection, Cl_vs_AileronDeflection_fit)
title('aileron effectiveness')
xlabel('aileron deflection [deg]')
ylabel('Cl about CG')
grid on
Cl_AileronDeflection_legend = sprintf('aileron effectiveness: %.3f per deg', aileron_effectiveness);
legend('WT test data', Cl_AileronDeflection_legend)

%--------------------------------------------------------------------------
%                       plot propeller influence
%--------------------------------------------------------------------------

% prop effects on elevator:
figure()
subplot(2,1,1)
plot(throttle_percentage, MaxElevator_throttle_sweep, '--ob')
grid on
title('prop effect on elevator at max deflection')
xlabel('engine throttle [%]')
ylabel('Cm max elevator deflection')
legend_text = sprintf('WT test data\nCm loss: %.3f', MaxElevator_Cm_loss);
legend(legend_text)

subplot(2,1,2)
plot(throttle_percentage, MinElevator_throttle_sweep, '--ob')
grid on
title('prop effect on elevator at min deflection')
xlabel('engine throttle [%]')
ylabel('Cm min elevator deflection')
legend_text = sprintf('WT test data\nCm loss: %.3f', MinElevator_Cm_loss);
legend(legend_text)

% prop effects on rudder:
figure()
subplot(2,1,1)
plot(throttle_percentage, MaxRudder_throttle_sweep, '--ob')
grid on
title('prop effect on rudder at max deflection')
xlabel('engine throttle [%]')
ylabel('Cn max rudder deflection')
legend_text = sprintf('WT test data\nCn loss: %.3f', MaxRudder_Cn_loss);
legend(legend_text)

subplot(2,1,2)
plot(throttle_percentage, MinRudder_throttle_sweep, '--ob')
grid on
title('prop effect on rudder at min deflection')
xlabel('engine throttle [%]')
ylabel('Cn min rudder deflection')
legend_text = sprintf('WT test data\nCn loss: %.3f', MinRudder_Cn_loss);
legend(legend_text)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               CLEAN UP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all open files:
fclose all;


