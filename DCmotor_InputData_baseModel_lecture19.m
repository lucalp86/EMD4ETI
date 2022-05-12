clc;
close all;
clear;

Data = readtable('MotorData.xlsx');
clc;
motorID = 1;

Power = Data{4, 5 + 3 * (motorID-1)};
Voltage = Data{5, 5 + 3 * (motorID-1)};
Speed_rpm = Data{6, 5 + 3 * (motorID-1)};
Speed_radsec = Data{7, 5 + 3 * (motorID-1)};

Torque = Data{8, 5 + 3 * (motorID-1)};
Current = Data{9, 5 + 3 * (motorID-1)};
kE = Data{11, 5 + 3 * (motorID-1)};
kT = Data{12, 5 + 3 * (motorID-1)};

Ra = Data{15, 5 + 3 * (motorID-1)};
La = Data{16, 5 + 3 * (motorID-1)};

Inertia = Data{19, 5 + 3 * (motorID-1)};
ViscousDamping = Data{20, 5 + 3 * (motorID-1)};

Inertia_Load = Inertia;

% Dynamic System
Ra_20 = Ra;
Jeq = Inertia + Inertia_Load;
k_d = 0; % ViscousDamping;
La = La * 10^-3;

% Time constant
tau_a = La / Ra_20;
disp(['Electrical Time constant=', num2str(tau_a),'sec'])
tau_m = Ra_20 * Jeq / (kT * kE); 
disp(['Mechanical Time constant=', num2str(tau_m),'sec'])

% inital condition of integrator
pos0 = 0;
w0 = 0;
i0 = 0;

% Solve the model
T_sim = 15;
T_step = 0.01;
Sol = sim('DCmotor_EMD4ETI22_Simulink_baseModel_lecture19', 'StartTime','0','StopTime', num2str(T_sim),'FixedStep', num2str(T_step));


figure;
subplot(2,2,1);
hold all;
plot(Sol.tout, Sol.DCmotor_In_Out.signals(1).values, 'b', 'LineWidth',2);
% scatter(t_on, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
% scatter(t_load, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
plot(Sol.tout, Voltage * ones(size(Sol.tout)), 'b--')
xlabel('Time [sec]');
ylabel('Armature Voltage [V]')
grid on;

subplot(2,2,2);
hold all;
plot(Sol.tout, Sol.DCmotor_In_Out.signals(3).values, 'Color',[0.3 0.75 0.9], 'LineWidth',2);
% scatter(t_on, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
% scatter(t_load, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
xlabel('Time [sec]');
ylabel('Armature Current [A]')
grid on;

TT = zeros(2,1);
subplot(2,2,3);
hold all;
TT(1) = plot(Sol.tout, Sol.DCmotor_In_Out.signals(2).values(:, 2), 'r', 'LineWidth',2);
TT(2) = plot(Sol.tout, Sol.DCmotor_In_Out.signals(2).values(:, 1), 'black', 'LineWidth',2, 'linestyle','--');
legend(TT, {'Load Torque', 'Electromagnetic Torque'});
% scatter(t_on, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
% scatter(t_load, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
xlabel('Time [sec]');
ylabel('Load Torque [Nm]')
grid on;

subplot(2,2,4);
hold all;
plot(Sol.tout, Sol.DCmotor_In_Out.signals(4).values, 'Color',[0.9 0.7 0.15], 'LineWidth',2);
plot(Sol.tout, Speed_rpm * ones(size(Sol.tout)), 'r--')
% scatter(t_on, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
% scatter(t_load, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
xlabel('Time [sec]');
ylabel('Mechanical Speed [rpm]')
grid on;



