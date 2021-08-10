%% Code to generate Band daigram with Oxide Charges


%%
clear;
clc;

%% Input

V_gs = 1;
Vright = 0;

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
N_sub = -3*1e17*1e6;           % Substrate_doping in /m^3 
                             %-ve for p-substrate(NMOS), +ve for n-substrate(PMOS)

k_ox = 4;
Eg_ox = 9*q;
eps_ox = k_ox*eps_0;
Ec_off = 3.1*q;                             
%psi_s = -0.4085             % Surface_band_bending in V  %Option 2 

%% Parameters

% Substrate Parameters
k_si = 12;                  %realtive permitivitty of silicon substrate
ni = 1.5e10*1e6;              %Intrinsic carrier density ( /cm^3)    
Eg = 1.1*q;                 %Band gap of silicon
Ea=4.05*q;                  %electron affinity for substrate
L_sub = 1e-6;
eps_si = k_si*eps_0;        %permitivitty of silicon substrate
chi_si = 4.05*q;

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


% Solver parameters
dx = 1e-9;              % Meshing size
NRmax = 10;             % Maximum NR iterations
tol = 1e-6;             % Error tolerance in potential


%Mesh
x = 0:dx:L_sub;
N = length(x);

%% Mesh

x_oxide = -t_ox:dx:0;
xsi = dx:dx:L_sub;
x = [x_oxide,xsi];
N = length(x);
index_oxide = x<=0;
index_si = x>0;
index_0 = find(x == 0);

%% Newton Raphson Matrices (Reference - Code presented by TA during live interaction)
ki = k_ox*index_oxide + k_si*index_si;
kip = 0.5*(ki + circshift(ki , -1));
kim = 0.5*(ki + circshift(ki , 1));
kipm = kip + kim;
A = diag(kipm, 0) - diag(kip(1:end-1) , 1) - diag(kim(2:end),-1);
A(1,:) = 0; A(1,1) = 1;
A(N,:) = 0; A(N,N) = 1;

b = zeros(N,1); 
psi = zeros(N,1);
%% Fixed Oxide charges
rho_ox = zeros(size(x_oxide)); 
size_rho_ox = size(rho_ox);
rho_ox(floor(size_rho_ox(2)/2)) = 2e12*1e4/dx*q; 
rho_oxfx = rho_ox;

%% Newton Raphson - Iterations

% Voltage Sweep
for Vg = V_gs 
    % Inputs at each voltage   
    V1 = Vg - V_fb;
    Qsub = [];
    
    % Bias point
    Vleft = V1;

    for i = 1:NRmax
        psi_si = psi(index_si);
        psi_s = psi(index_0);
        E1 = q*psi_s + Eg/2 - q*phi_b;
        E2 = Eg/2;
   
        p = psub*exp(-psi_si/Vt);
        n = nsub*exp(psi_si/Vt);
        delp = -p/Vt;
        deln = n/Vt;
        
        if (sign(N_sub) == 1)
            majc = n;
            delmajc = deln;
            minc = p;
            delminc = delp;
        elseif (sign(N_sub) == -1)
            majc = p;
            delmajc = delp;
            minc = n;
            delminc = deln;
        end

        rho_ox(end) = rho_oxfx(end);
        rho_si = q*(N_sub + sign(N_sub)*(minc - majc));
        
        rho  = [rho_ox' ; rho_si];
        b = rho*dx^2 /eps_0;
        b(1) = Vleft;
        b(N) = Vright;
        f = A*psi-b;
        
        % Jacobian Calculations
        delrho_ox = zeros(size(rho_ox));
        delrho_ox(end) = delrho_ox(end);
        
        delrho_si = q*(sign(N_sub) * (delminc - delmajc));
        
        delrho = [delrho_ox' ; delrho_si];      
        delb = delrho*dx^2/eps_0;
        delb(1) = 0;
        delb(N) = 0;
        J = A-diag(delb);
        dV = -J\f;
        
        if (max(abs(dV))<tol)
            break;
        end
        
        psi = psi+dV;
    end
    
    Q1 = sum(rho_si)*dx ;
    minc1 = minc;
    
   
end


    

%% Band Diagram

% Here we will take the intrinsic level in the substrate to be reference

% Energy bands in substrate

Ei_band = -q*psi(index_si);
Ec_band = Ei_band + Eg/2;               % Conduction band
Ev_band = Ei_band - Eg/2;               % Valence band
Ef_band = (Ei_band(end) - q*phi_b)*ones(length(Ei_band));       % Fermi level

% Energy bands in Oxide
Ec_ox_band = -q*psi(index_oxide) + Ec_off+ Eg/2 ;
Ev_ox_band = Ec_ox_band - Eg_ox;

% Energy band in metal
x_m = -t_ox-t_m:dx:-t_ox;               % Meshing in metal layer
Em_band = (Ef_band(end) - q*Vg)*ones(size(x_m));          % Because q*Vg = Em - Ef

figure(1); plot([x_m, x_oxide, xsi]*1e9 ,[Em_band' ; Ec_ox_band ; Ec_band]/q , 'r',[x_m x_oxide xsi]*1e9 , [Em_band'; Ev_ox_band; Ev_band]/q ,'b',xsi*1e9,Ei_band/q,'k--',xsi*1e9,Ef_band/q,'k','LineWidth' , 2);
xlim([-20,150])
axes.LineWidth = 1; Axes.FontSize = 18; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('x (nm) -->'); ylabel('Energy (eV)');
leg = legend('E_c','E_v','E_i','E_f');
title("Energy Band Diagram for MOS capacitors")




