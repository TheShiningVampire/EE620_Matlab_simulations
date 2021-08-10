clear;
clc;

% Constants
q = 1.6e-19;
eps_0 = 8.85e-12;
Vt = 26e-3;               %thermal voltage
k = 8.314/(6.023*10^23);
T = 300;

%% Parameters

% Substrate Parameters
k_si = 12;                  %realtive permitivitty of silicon substrate
ni = 1.5e10*1e6;              %Intrinsic carrier density ( /cm^3)    
Eg = 1.1*q;                 %Band gap of silicon
Ea=4.05*q;                  %electron affinity for substrate
L_sub = 1e-6;
eps_si = k_si*eps_0;        %permitivitty of silicon substrate
chi_si = 4.05*q;

k_ox = 4;
Eg_ox = 9*q;
eps_ox = k_ox*eps_0;
Ec_off = 3.1*q;  


%%
data = csvread('run1.csv');

vgs = data(:,1);
caps = data(:,2);

figure(1); plot(vgs,caps); hold on

caps= caps * 1e-6*1e4;

cox = caps(1);
cmin = caps(end);

t_ox = eps_ox/cox;

cdep = 1/( 1/cmin - 1/cox );
wdep = eps_si / cdep ; 


ni = 1.5e10*1e6;
% Initial guess
Na = 1e17*1e6;

for i = 1:10
    phi_b = Vt*log(Na/ni);
    Na = 2*eps_si*phi_b/q/wdep^2;
end

phi_b = Vt*log(Na/ni);
Eg = q*1.1;
Vfb = -Eg/(2*q) - phi_b;

W_mg = (2*eps_si*phi_b/q/Na)^0.5;
Cfb = 1/(1/cox + wi/eps_si);


