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
k_d = ViscousDamping;
La = La * 10^-3;

% Time constant
tau_a = La / Ra_20;
disp(['Electrical Time constant=', num2str(tau_a),'sec'])
tau_m = Ra_20 * Jeq / (kT * kE); 
disp(['Mechanical Time constant=', num2str(tau_m),'sec'])

% ----- DC-DC buck converter
% V_DC = 40; % [V] DC link voltage
V_DC = 1.5 * Voltage; % [V] DC link voltage
%f_switching = 2; % [kHz] switching frequency
f_switching = Speed_radsec/(2*pi) * 30 / 1000; % [kHz] switching frequency
T_triang = 1/(1000 * f_switching); % [sec] Triangular waveform period

% ----- Power Supply
Rs = 0.1;
Cs = 0.0003;
V_S = V_DC;

% inital condition of integrator
pos0 = 0;
w0 = 0;
i0 = 0;
vC0 = V_DC;

% Solve the model
T_sim = 15;
T_step = T_triang/100;
Sol = sim('DCmotor_EMD4ETI22_Simulink_baseModel_lecture19_PEC_Supply', 'StartTime','0','StopTime', num2str(T_sim),'FixedStep', num2str(T_step));


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



hhhhhhhhhhhhhhhhhh

Rth = 0.1;
Cth = 10;
Tamb = 20;
alpha_cu = 0.005;

% ---- Data manipulation
kT = kE;
% kE = kT;
Jeq = 2 * J;
Ra_20 = Ra;
La = La * 10^-3;   % Armature Inductance [H]
Ia_max = Pn / (wn * pi / 30) / kT;

tau_a = La / Ra_20;
tau_m = Ra_20 * Jeq / (kT * kE); 



Vdc = 60;
f_triang = 5000;
T_triang = 1/f_triang;

V_applied = 24;
T_Load = 0;

t_on = 0.5;
V_sec= 24;
V_Lim = Vn;

t_load = 5;
Delta_Load = 0.5;
T_Lim = Tn;
T_sec = Delta_Load;

t_va = 8.5;
Delta_V = Vn;

pos0 = 0;


T_sim = 15; %7000; %250 * max(tau_a, tau_m);
% Sol = sim('DCmotor_InputData_lecture17_EMD4ETI22_Model_Tot', T_sim);
% Sol = sim('DCmotor_InputData_lecture17_EMD4ETI22_Model_Thermal_AdvMech', T_sim);
% Sol = sim('DCmotor_InputData_lecture17_EMD4ETI22_Model_Thermal_GenLoad', T_sim);
% Sol = sim('DCmotor_InputData_lecture17_EMD4ETI22_baseModel_Thermal_PE', T_sim);
% Sol = sim('DCmotor_InputData_lecture17_EMD4ETI22_baseModel_Thermal', T_sim);
Sol = sim('DCmotor_InputData_lecture19_EMD4ETI22_baseModel', T_sim);

lllllllllllllllllllllllll

figure;
hold all;
subplot(2,2,1);
hold all;
plot(Sol.tout, Sol.ScopeDataIN.signals(1).values, 'b', 'LineWidth',2);
scatter(t_on, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
scatter(t_load, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
plot(Sol.tout, Vn * ones(size(Sol.tout)), 'b--')
xlabel('Time [sec]');
ylabel('Armature Voltage [V]')

subplot(2,2,3);
hold all;
plot(Sol.tout, Sol.ScopeDataIN.signals(2).values, 'r', 'LineWidth',2);
scatter(t_on, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
scatter(t_load, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
xlabel('Time [sec]');
ylabel('Load Torque [Nm]')

subplot(2,2,2);
hold all;
plot(Sol.tout, Sol.ScopeDataOUT.signals(1).values, 'b', 'LineWidth',2);
scatter(t_on, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
scatter(t_load, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
xlabel('Time [sec]');
ylabel('Armature Current [A]')

subplot(2,2,4);
hold all;
plot(Sol.tout, Sol.ScopeDataOUT.signals(2).values, 'r', 'LineWidth',2);
plot(Sol.tout, wn * ones(size(Sol.tout)), 'r--')
scatter(t_on, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
scatter(t_load, 0, MarkerFaceColor = 'white', MarkerEdgeColor='black')
xlabel('Time [sec]');
ylabel('Mechanical Speed [rpm]')


figure;
subplot(2,1,1);
hold all;
A(1) = plot(Sol.Pel.time, Sol.Pel.signals.values, 'green', 'DisplayName', 'P electrical', 'LineWidth',2);
A(2) = plot(Sol.Pel.time, Sol.Pmec.signals.values, 'red', 'DisplayName', 'P mechanical', 'LineWidth',2);
A(3) = plot(Sol.Pel.time, Sol.PJoule.signals.values, 'blue', 'DisplayName', 'P Joule', 'LineWidth',2);
%A(4) = plot(Sol.Pel.time, Sol.Pel.signals.values - Sol.Pmec.signals.values, 'color', 'black', 'DisplayName', 'd/dt E', 'LineWidth',2);
legend(A);
xlabel('Time [sec]');
ylabel('Power [W]')

subplot(2,1,2);
hold all;
plot(Sol.Pel.time, Sol.Temp_arm.signals.values, 'black', 'LineWidth',2);
plot(Sol.Pel.time, Tamb * ones(size(Sol.tout)), 'black', 'LineWidth',2, LineStyle='--');
xlabel('Time [sec]');
ylabel('Temperature Armature [C]')



% 
% 
% 
% 
% 
% 
% 
% 
% % ----- Control Analysis
% w0 = 1 / sqrt(tau_a * tau_m);
% xi = 1 / (2 * w0 * tau_a);
% 
% % ----- Control Settings
% se = - xi * w0 * (1 + sqrt(1 - 1/xi^2));
% sm = - xi * w0 * (1 - sqrt(1 - 1/xi^2));
% 
% kP = 200;
% kI = kP * sm;
% 
% % ----- Input
% T_Load = 0.2;
% Vdc = 24;
% 
% Wm_ref = 1000;
% 
% wm_ref = Wm_ref * pi / 30;
% T_sim = 10 * max(tau_a, tau_m);
% Sol = sim('BlockDiagram_DCmotor_lecture20_control_converter', T_sim);
% 
% figure;
% hold all;
% plot(Sol.tout, Sol.yout.signals(1).values * 30 / pi, 'r');
% plot(Sol.tout, Wm_ref * ones(size(Sol.tout)), 'r--')
% xlabel('Time');
% ylabel('Speed [rpm]')
% 
% figure;
% hold all;
% plot(Sol.tout, Sol.yout.signals(2).values, 'b');
% xlabel('Time');
% ylabel('Armature Current')
% 
% figure;
% hold all;
% A(1) = plot(Sol.P_resp.time, Sol.P_resp.signals.values, 'r', 'DisplayName', 'P - response');
% A(2) = plot(Sol.I_resp.time, Sol.I_resp.signals.values, 'b', 'DisplayName', 'I - response');
% A(3) = plot(Sol.PI_resp.time, Sol.PI_resp.signals.values, 'color', [0, 0.7,0], 'DisplayName', 'PI - response');
% legend(A);
% xlabel('Time');
% ylabel('Controller Signal')

