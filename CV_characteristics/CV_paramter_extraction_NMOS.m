clear;
close all;
clc;

% Experminetal Data
data = csvread('./data/CV_run1.csv');
V = data(:,1);
C = data(:,2);
plot(V,C);
grid on;
xlabel("Voltage (V)");
ylabel("Capacitance (uF/cm^2)");

% Paramters
eps_ox = 3.9*8.85e-12;
KT_q = 25e-3;
q = 1.6e-19;
n_i = 1e10;

% Parameter extraction

% 1) Accumulation CV
C_max = C(1);
C_ox = C_max*1e-2; % in units of F/m^2

T_ox = eps_ox/C_ox; % in m.


% 2) Inversion CV

C_min = C(end)*1e-2;  % in units of F/m^2

Na1 = 1e17; 
Na2 = 1;

while (abs(Na1-Na2)> 1e-8)
    Na2 = Na1;
    phi = KT_q*log(Na2/n_i);
    Na1 = 4*phi/(q*eps_ox*((1/C_min) - (1/C_ox))^2);
end
 
N_a = Na1;   % /m^3

%3) Midgap CV
C_mg = 1/(1/C_ox + sqrt(2*eps_ox*phi/(q*N_a))/eps_ox );
C_mg = C_mg*1e2; % in units of uF/cm^2
V_mg = 0;
for i = 1:length(C)
    if (abs(C(i)-C_mg) <= 0.002)
        V_mg = V(i);
    end
end

%4) Flatband CV
C_fb = 1/(1/C_ox + sqrt(eps_ox*KT_q/(q*N_a))/eps_ox );
V_fb = 0;
for i = 1:length(C)
    if (abs(C(i)-C_fb) <= 0.002)
        V_fb = V(i);
    end
end

phi_m = V_fb - phi;

fprintf("T_ox = %i m \n", T_ox)
fprintf("Doping density = %i /m^3\n", N_a)
fprintf("Midgap voltage = %i V\n", V_mg)
fprintf("phi_m = %i \n", phi_m)















