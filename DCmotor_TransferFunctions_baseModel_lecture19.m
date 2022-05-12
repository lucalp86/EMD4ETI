clc;
close all;
clear;

Data = readtable('MotorData.xlsx');
clc;

motorID = 3;

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
k_d = 0; %ViscousDamping;
La = La * 10^-3;

A = [-Ra_20/La -kE/La;
    kT/Jeq -k_d/Jeq];

B = [1/La 0;
    0 -1/Jeq];

C = [1 0;
    0 1];

D = 0;

sys = ss(A,B,C,D);

figure;
subplot(2,2,1)
step(sys(1,1));
title('H_{ee}', 'bold');
grid on;
subplot(2,2,2)
step(sys(1,2))
title('H_{em}', 'bold');
grid on;
subplot(2,2,3)
step(sys(2,1))
title('H_{me}', 'bold');
grid on;
subplot(2,2,4)
step(sys(2,2))
title('H_{mm}', 'bold');
grid on;


figure;
subplot(2,2,1)
bode(sys(1,1));
title('H_{ee}', 'bold');
grid on;
subplot(2,2,2)
bode(sys(1,2))
title('H_{em}', 'bold');
grid on;
subplot(2,2,3)
bode(sys(2,1))
title('H_{me}', 'bold');
grid on;
subplot(2,2,4)
bode(sys(2,2))
title('H_{mm}', 'bold');
grid on;


figure;
subplot(2,2,1)
rlocus(sys(1,1));
title('H_{ee}', 'bold');
grid on;
subplot(2,2,2)
rlocus(sys(1,2))
title('H_{em}', 'bold');
grid on;
subplot(2,2,3)
rlocus(sys(2,1))
title('H_{me}', 'bold');
grid on;
subplot(2,2,4)
rlocus(sys(2,2))
title('H_{mm}', 'bold');
grid on;

% s = tf('s');
% Den = (s+Ra_20/La)*(s+k_d/Jeq) + kT * kE / (La * Jeq);

% ----- Hee
% get one of the transfer function
U = tf(sys(1,1))
% determine the numerator and denumerator
[num,den] = tfdata(U);
% compute the zeros (z), poles (p), and gain (k)
[z,p,k] = tf2zp(num{1},den{1})

% ----- Hem
% get one of the transfer function
U = tf(sys(1,2))
% determine the numerator and denumerator
[num,den] = tfdata(U);
% compute the zeros (z), poles (p), and gain (k)
[z,p,k] = tf2zp(num{1},den{1})

% ----- Hmm
% get one of the transfer function
U = tf(sys(2,2))
% determine the numerator and denumerator
[num,den] = tfdata(U);
% compute the zeros (z), poles (p), and gain (k)
[z,p,k] = tf2zp(num{1},den{1})

