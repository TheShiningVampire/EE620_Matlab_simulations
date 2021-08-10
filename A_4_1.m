%% Code to generate HFCV and LFCV curves for NMOS and PMOS for fixed oxide charge and interface trap densities
%% Inputs : Gate workfunction, oxide thickness, channel doping, sweep_type
%% sweep_type = 0 for LFCV
%%            = 1 for HFCV
%%            = 2 for fast sweep   



%%
clear;
clc;

%% Input

V_gs = (-1.9:50e-3:2.5);
Vwiggle = 1e-3;
Vright = 0;
sweep_type = 1; 
dit_dpsi = eps;



% Constants
q = 1.6e-19;
eps_0 = 8.85e-12;
Vt = 26e-3;               %thermal voltage
k = 8.314/(6.023*10^23);
T = 300;

%% Inputs
% Give the inputs here (if default values are not used)

phi_m = 3.9*q;               % Metal_Workfunction in eV 
t_ox = 2e-9;                % Oxide_Thickness in m 
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
dx = 1e-9;             % Meshing size
NRmax = 10;             % Maximum NR iterations
tol = 1e-6;             % Error tolerance in potential


%Mesh
x = 0:dx:L_sub;
N = length(x);

%% Mesh

x_oxide = -t_ox:dx:0;
x_si = dx:dx:L_sub;
x = [x_oxide,x_si];
N = length(x);
i_oxide = x<=0;
index_si = x>0;
index_0 = find(x == 0);

%% Newton Raphson Matrices (Reference - Code presented by TA during live interaction)
ki = k_ox*i_oxide + k_si*index_si;
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
rho_ox(end) = 2e12*1e4/dx*q; 
rho_oxfx = rho_ox;
%% Interface charges
dit  = 3e12  *1e4/q;
ditf = @(E) dit*ones(size(E));
%% Newton Raphson - Iterations
C = [];
psi_ss = [];
% Voltage Sweep
for Vg = V_gs  
    disp(Vg);
    % Inputs at each voltage   
    V1 = Vg - V_fb;
    V2 = V1 +Vwiggle;
    Qsub = [];
    
    % Bias point
    Vleft = V1;
    for i = 1:NRmax
        psi_si = psi(index_si);
        psi_s = psi(index_0);
        E1 = q*psi_s + Eg/2 - q*phi_b;
        E2 = Eg/2;
        rho_dit = q/dx*integral(ditf , E1, E2);
        rho_dit_dpsi = q/dx * integral(ditf , E1 + q*dit_dpsi,E2);
        delrho_dit = 1/dit_dpsi * (rho_dit_dpsi - rho_dit);
        
        p = psub*exp(-psi_si/Vt);
        n = nsub*exp(psi_si/Vt);
        delp = -p/Vt;
        deln = n/Vt;
        
        if (sweep_type == 2)
            minc = 0;
            delminc = 0;
            ditq = 0;
            ditdq = 0;
            if (sign(N_sub) == 1)
                majc = n;
                delmajc = deln;
            elseif (sign(N_sub) == -1)
                majc = p;
                delmajc = delp;
            end
        else
            ditq = rho_dit;
            ditdq = delrho_dit;
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
        end
        
        rho_ox(end) = rho_oxfx(end) + ditq;
        rho_si = q*(N_sub + sign(N_sub)*(minc - majc));
        
        rho  = [rho_ox' ; rho_si];
        b = rho*dx^2 /eps_0;
        b(1) = Vleft;
        b(N) = Vright;
        f = A*psi-b;
        
        % Jacobian Calculations
        delrho_ox = zeros(size(rho_ox));
        delrho_ox(end) = delrho_ox(end) + ditdq;
        
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
    
    Q1 = sum(rho_si)*dx + ditq*dx;
    minc1 = minc;
    ditq1 = ditq;
    
    %Wiggle point
    Vleft = V2;
    
    for i= 1:NRmax
        psi_si = psi(index_si);
        psi_s = psi(index_0);
        E1 = q*psi_s + Eg/2 - q*phi_b;
        E2 = Eg/2;
        rho_dit = q/dx*integral(ditf , E1, E2);
        rho_dit_dpsi = q/dx * integral(ditf , E1 + q*dit_dpsi,E2);
        delrho_dit = 1/dit_dpsi * (rho_dit_dpsi - rho_dit);
        
        p = psub*exp(-psi_si/Vt);
        n = nsub*exp(psi_si/Vt);
        delp = -p/Vt;
        deln = n/Vt;
        
        if (sweep_type == 2)
            minc = 0;
            delminc = 0;
            ditq = 0;
            ditdq = 0;
            if (sign(N_sub) == 1)
                majc = n;
                delmajc = deln;
            elseif (sign(N_sub) == -1)
                majc = p;
                delmajc = delp;
            end
        elseif sweep_type ==1
            minc = minc1;
            delminc = 0;
            ditq = ditq1;
            ditdq = 0;
            if (sign(N_sub) == 1)
                majc = n;
                delmajc = deln;
            elseif (sign(N_sub) == -1)
                majc = p;
                delmajc = delp;
            end
        else
            ditq = rho_dit;
            ditdq = delrho_dit;
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
        end
        
        rho_ox(end) = rho_oxfx(end) + ditq;
        rho_si = q*(N_sub + sign(N_sub)*(minc - majc));
        
        rho  = [rho_ox' ; rho_si];
        b = rho*dx^2 /eps_0;
        b(1) = Vleft;
        b(N) = Vright;
        f = A*psi-b;
        
        % Jacobian Calculations
        delrho_ox = zeros(size(rho_ox));
        delrho_ox(end) = delrho_ox(end) + ditdq;
        
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
    
    Q2 = sum(rho_si)*dx + ditq*dx;
    
    psi_ss = [psi_ss psi_s];
    Cap = (Q1 - Q2)/Vwiggle;
    C = [C Cap];   
end


    
figure(1);
plot(V_gs,C,'LineWidth',1.5);

xlabel('Voltage(in V)');
ylabel('Capacitance(in F/cm^2)');
grid on;
ax=gca;
ax.FontSize=10;
hold off;










