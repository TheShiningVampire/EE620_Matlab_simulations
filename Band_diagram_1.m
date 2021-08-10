% This Poisson Solver calculates electrostatic parameters and plots the
% Energy Band Diagram based on the following inputs
% a) Metal Workfunction
% b) Oxide Thickness
% c) Substrate doping
% d) psi_s (Surface band bending in the substrate)

clear;
close all;
clc;

%% Inputs
% Give the inputs here (if default values are not used)

phi_m = -0.2532e-19         % Metal_Workfunction in eV 
t_ox = 10e-9                % Oxide_Thickness in m 
N_sub = -1e17*1e6           % Substrate_doping in /m^3 
                            %-ve for p-substrate(NMOS), +ve for n-substrate(PMOS)

psi_s = -0.4085              % Surface_band_bending in V  %Option 2 

%% Parameters

% Constants
q = 1.6e-19;
eps_0 = 8.85e-12;
Vt = 26e-3;               %thermal voltage

% Substrate Parameters
k_si = 12;                  %realtive permitivitty of silicon substrate
ni = 1.5e10*1e6;            %Intrinsic carrier density ( /cm^3)    
Eg = 1.1*q;                 %Band gap of silicon
Ea=4.05*q;                  %electron affinity for substrate
L_sub = 1e-6;
eps_si = k_si*eps_0;        %permitivitty of silicon substrate

%calculate mobile carriers (n and p) for substrate
if N_sub< 0 
    psub = abs(N_sub); nsub = ni^2/abs(N_sub);
elseif N_sub>0
    psub = ni^2/N_sub; nsub = N_sub;
else
    psub = ni; nsub = ni;
end

% Oxide parameters
k_ox = 4;
Eg_ox = 9*q;                % band gap of oxide
eps_ox = k_ox*eps_0;        % permitivitty of SiO2
Ec_off = 3.1*q;             % conduction band offset

% Metal parameters
t_m = 10e-9;
phi_b = -sign(N_sub)*Vt*log(abs(N_sub)/ni);      %Ei-Ef for semiconductor

V_fb = sign(N_sub)*Eg/(2*q) - phi_b;             %flatband voltage

% Solver
dx = 10e-9;             % Meshing size
NRmax = 10;             % Maximum NR iterations
tol = 1e-6;             % Error tolerance in potential

%% Solver
%Mesh
x = 0:dx:L_sub;
N = length(x);

%Newton Raphson Matrices

A = -2*diag(ones(1,N) , 0)+ 1*diag(ones(1,N-1),1)+ 1*diag(ones(1,N-1),-1);
A(1,:) = 0; A(1,1) = 1;
A(N,:) = 0; A(N,N) = 1;

b = zeros(N,1);
b(1) = psi_s; b(N) = 0;

psi = zeros(N,1);               % Option 1 for defining the initial guess
%psi = A\b;                     % Option 2 for defining the initial guess

%Newton Raphson iterations
for i = 1:NRmax
    
    p = psub*exp(-psi/Vt);
    n = nsub*exp(psi/Vt);
    rho = q*(N_sub + p - n);
    b = -rho/eps_si * dx^2;
    b(1) = psi_s; b(N) = 0;
    f = A*psi - b;
    
    % Jacobian Matrix
    delp = -1/Vt * p;
    deln = 1/Vt * n;
    delrho = q*(delp-deln);
    delb = -delrho/eps_si * dx^2;
    delb(1) = 0; delb(N) = 0;
    J = A - diag(delb);
    dV = -J\f;
    
    if max(abs(dV))<tol            % If the change is lower than the tolerance then stop the iterations
        break;
    end
    
    psi = psi+dV;
end


%% Addition of Oxide

E_surface = -(psi(2) - psi(1))/dx;      % Electric field at surface = -d(psi)/dx
E_ox = E_surface*(eps_si/eps_ox);       % Equal displacement vectors across the interface
x_ox = -t_ox:dx:0;                      % Mesh in oxide region
psi_ox = psi(1) - E_ox*x_ox;
V_ox = E_ox*t_ox;                       % Voltage drop across the oxide

%% Band Diagram

% Here we will take the intrinsic level in the substrate to be reference

% Energy bands in substrate
Ei_band = -q*psi;
Ec_band = Ei_band + Eg/2;               % Conduction band
Ev_band = Ei_band - Eg/2;               % Valence band
Ef_band = (Ei_band(end) - q*phi_b)*ones(length(Ei_band));       % Fermi level

% Energy bands in Oxide
Ec_ox_band = -q*psi_ox + Ec_off+ Eg/2;
Ev_ox_band = Ec_ox_band - Eg_ox;

Vg = V_fb + V_ox + psi_s;               % Voltage drop
Em_band = (Ef_band(end) - q*Vg)*ones(size(Ec_ox_band));          % Because q*Vg = Em - Ef
x_m = -t_ox-t_m:dx:-t_ox;               % Meshing in metal layer

figure(1); plot([x_m x_ox x]*1e9 ,[Em_band' ; Ec_ox_band' ; Ec_band]/q , 'r',[x_m x_ox x]*1e9 , [Em_band'; Ev_ox_band'; Ev_band]/q ,'b',x*1e9,Ei_band/q,'k--',x*1e9,Ef_band/q,'k','LineWidth' , 2);
xlim([-20,150])
axes.LineWidth = 1; Axes.FontSize = 18; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('x (nm) -->'); ylabel('Energy (eV)');
leg = legend('E_c','E_v','E_i','E_f');
title("Energy Band Diagram for MOS capacitors")




    




    