% This Poisson Solver calculates electrostatic parameters and plots the
% Energy Band Diagram based on the following inputs
% a) Metal Workfunction
% b) Oxide Thickness
% c) Substrate doping
% d) psi_s (Surface band bending in the substrate)

clear;
close all;
clc;

% Constants
q = 1.6e-19;
eps_0 = 8.85e-12;
Vt = 26e-3;               %thermal voltage
k = 8.314/(6.023*10^23);
T = 300;
%% Inputs
% Give the inputs here (if default values are not used)

phi_m = 4.1*q;               % Metal_Workfunction in eV 
t_ox = 10e-9;                % Oxide_Thickness in m 
N_sub = -1e17*1e6;           % Substrate_doping in /m^3 
                             %-ve for p-substrate(NMOS), +ve for n-substrate(PMOS)

V_g = 2; 
%psi_s = -0.4085             % Surface_band_bending in V  %Option 2 

%While using option 2 comment out lines 64-71 and remove the comment from
%line 146

%% Parameters

% Substrate Parameters
k_si = 12;                  %realtive permitivitty of silicon substrate
ni = 1e10*1e6;              %Intrinsic carrier density ( /cm^3)    
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
c_ox=eps_ox/t_ox;

% Metal parameters
t_m = 10e-9;
phi_b = -sign(N_sub)*Vt*log(abs(N_sub)/ni);         %Ei-Ef for semiconductor
phi_sub=(Ea+Eg/2+phi_b*q);                          %work function for substrate
V_fb = (phi_m-phi_sub)/q;                           %flatband voltage
V_th = V_fb -sign(N_sub)*(4*q*eps_si*abs(N_sub)*abs(phi_b))^0.5/c_ox + 2*phi_b;

%%Calculate psi_s from V_g

% psi_init= -(V_g - V_fb)*sign(N_sub);
% pre=(2*eps_si*Vt*q*abs(N_sub))^0.5;
% 
% func = @(psi_s) -V_g + V_fb -(1/c_ox *-sign(V_g-V_fb)*pre*(exp(-psi_s/Vt)+psi_s/Vt-1+(ni^2/(abs(N_sub))^2)*(exp(psi_s/Vt)-psi_s/Vt-1))^0.5) + psi_s;
% x0 = psi_init;
% psi_s = fzero(func,x0);


% Sign of substrate charges in the current mode, since we have signum(psi_s) in formula of charge

if(V_g>=V_fb)
    signq = -1;       %Depletion for NMOS, Accumulation for PMOS
else                  % i.e. (V_g<V_fb)
    signq = +1;       %Accumulation for NMOS, Depletion for PMOS
end

%Same charge sign for depletion and inversion

if sign(N_sub) <= 0
    F = @(s) V_fb-(signq*((2*eps_si*k*T*abs(N_sub))^0.5*((exp(-q*s/(k*T))+q*s/(k*T)-1)+(ni^2/(abs(N_sub))^2)*(exp(+q*s/(k*T))-q*s/(k*T)-1))^0.5)/c_ox)+s-V_g;
    psi_s = fsolve(F,-signq*0.1); % Initial guess, sign taken due to fsolve
    
else                              % sign(N_sub) > 0
    F = @(s) V_fb-(signq*((2*eps_si*k*T*abs(N_sub))^0.5*((ni^2/(abs(N_sub))^2)*(exp(-q*s/(k*T))+q*s/(k*T)-1)+(exp(+q*s/(k*T))-q*s/(k*T)-1))^0.5)/c_ox)+s-V_g;
    psi_s = fsolve(F,-signq*0.1); % Initial guess, sign taken due to fsolve
end


%%

% Solver
dx = 1e-9;             % Meshing size
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
E_ox = -sign(N_sub)*E_surface*(eps_si/eps_ox);       % Equal displacement vectors across the interface
x_ox = linspace(-t_ox,0,10);                      % Mesh in oxide region
psi_ox = psi(1) - E_ox*x_ox;
V_ox = E_ox*t_ox;                       % Voltage drop across the oxide

%% Band Diagram

% Here we will take the intrinsic level in the substrate to be reference

% Energy bands in substrate
Ei_band = q*psi*sign(N_sub);
Ec_band = Ei_band + Eg/2;               % Conduction band
Ev_band = Ei_band - Eg/2;               % Valence band
Ef_band = (Ei_band(end)- q*phi_b)*ones(length(Ei_band));       % Fermi level

% Energy bands in Oxide
Ec_ox_band = -q*psi_ox + Ec_off+ Eg/2;
Ev_ox_band = Ec_ox_band - Eg_ox;

%V_g = V_fb + V_ox + psi_s;             % Voltage drop
Em_band = (Ef_band(end) - q*V_g)*ones(size(Ec_ox_band));          % Because q*Vg = Em - Ef
x_m = linspace(-t_ox-t_m,-t_ox,10);               % Meshing in metal layer

figure(1); plot([x_m x_ox x]*1e9 ,[Em_band' ; Ec_ox_band' ; Ec_band]/q , 'r',[x_m x_ox x]*1e9 , [Em_band'; Ev_ox_band'; Ev_band]/q ,'b',x*1e9,Ei_band/q,'k--',x*1e9,Ef_band/q,'k','LineWidth' , 2);
axes.LineWidth = 1; Axes.FontSize = 18; axes.FontWeight = 'bold'; axes.Box = 'on';
xlim([-20,200])
xlabel('x (nm) -->'); ylabel('Energy (eV)');
leg = legend('E_c','E_v','E_i','E_f');
title("Energy Band Diagram for MOS capacitors")

%% Plot of carrier density

if N_sub < 0
    maj = p;        % Majority carriers
    min = n;        % Minority carriers
else
    maj = n;        % Majority carriers
    min = p;        % Minority carriers
end

figure(2); plot(x*1e9 , maj , 'r','LineWidth' , 2);
xlim([0,200])
axes.LineWidth = 1; Axes.FontSize = 18; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('x (nm) -->'); ylabel('Carrier density (/m^3)');
title("Majority carrier density in substrate for MOS capacitors")

figure(3); plot(x*1e9 , min , 'b','LineWidth' , 2);
xlim([0,200])
axes.LineWidth = 1; Axes.FontSize = 18; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('x (nm) -->'); ylabel('Carrier density (/m^3)');
title("Minority carrier density in substrate for MOS capacitors")




    




    