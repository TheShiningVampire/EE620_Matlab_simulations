clear;
close all;
clc;


% Constants
q = 1.6e-19;
eps_0 = 8.85e-12;
Vt = 26e-3;               %thermal voltage
k = 8.314/(6.023*10^23);
T = 300;
V_wiggle=10e-3;

%% Inputs
% Give the inputs here (if default values are not used)

phi_m = 4.1*q;               % Metal_Workfunction in eV 
t_ox = 10e-9;                % Oxide_Thickness in m 
N_sub = -1e17*1e6;           % Substrate_doping in /m^3 
                             %-ve for p-substrate(NMOS), +ve for n-substrate(PMOS)

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


% Solver parameters
dx = 1e-9;             % Meshing size
NRmax = 10;             % Maximum NR iterations
tol = 1e-6;             % Error tolerance in potential


%Mesh
x = 0:dx:L_sub;
N = length(x);



%%

% V_g = -2.5:0.1:2.5;
% C=zeros(length(V_g),1);
% 
% for i=1:length(V_g)
%     
%    psi = Poisson_solver(V_g(i) , psub , nsub);
%    p = psub*exp(-psi/Vt);
%    n = nsub*exp(psi/Vt);
%    rho_si_dc = q*(N_sub + p - n);
%    Q_dc = sum(rho_si_dc) * dx;
%    
%    psi = Poisson_solver(V_g(i) + V_wiggle , psub , nsub);
%    p = psub*exp(-psi/Vt);
%    n = nsub*exp(psi/Vt);
%    rho_si_dcwiggle = q*(N_sub + p - n);
%    Q_dcwiggle = sum(rho_si_dcwiggle) * dx;   
%    C(i) = -(Q_dcwiggle - Q_dc)/V_wiggle;   
% end
% 
% plot(V_g , C)
% 
% 

%%
V_g = -2.5:0.1:2.5;
C=zeros(length(V_g),1);
Q1=zeros(length(V_g),1);
Q2=zeros(length(V_g),1);
for i=1:length(V_g)
    
    if(V_g(i)>=V_fb)
        signq = -1;       %Depletion for NMOS, Accumulation for PMOS
    else                  % i.e. (V_g(i)<V_fb)
        signq = +1;       %Accumulation for NMOS, Depletion for PMOS
    end

    %Same charge sign for depletion and inversion

    if sign(N_sub) <= 0
        F = @(s) V_fb-(signq*((2*eps_si*k*T*abs(N_sub))^0.5*((exp(-q*s/(k*T))+q*s/(k*T)-1)+(ni^2/(abs(N_sub))^2)*(exp(+q*s/(k*T))-q*s/(k*T)-1))^0.5)/c_ox)+s-V_g(i);
        psi_s = fsolve(F,-signq*0.1); % Initial guess, sign taken due to fsolve

    else                              % sign(N_sub) > 0
        F = @(s) V_fb-(signq*((2*eps_si*k*T*abs(N_sub))^0.5*((ni^2/(abs(N_sub))^2)*(exp(-q*s/(k*T))+q*s/(k*T)-1)+(exp(+q*s/(k*T))-q*s/(k*T)-1))^0.5)/c_ox)+s-V_g(i);
        psi_s = fsolve(F,-signq*0.1); % Initial guess, sign taken due to fsolve
    end


    % Solver parameters
    dx = 1e-9;             % Meshing size
    NRmax = 10;             % Maximum NR iterations
    tol = 1e-6;             % Error tolerance in potential


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

   rho_si_dc = q*(N_sub + p - n);
   Q1(i) = sum(rho_si_dc) * dx;

 
   
   
   
   
   
   
   
       if(V_g(i) + V_wiggle>=V_fb)
        signq = -1;       %Depletion for NMOS, Accumulation for PMOS
    else                  % i.e. (V_g(i) + V_wiggle<V_fb)
        signq = +1;       %Accumulation for NMOS, Depletion for PMOS
    end

    %Same charge sign for depletion and inversion

    if sign(N_sub) <= 0
        F = @(s) V_fb-(signq*((2*eps_si*k*T*abs(N_sub))^0.5*((exp(-q*s/(k*T))+q*s/(k*T)-1)+(ni^2/(abs(N_sub))^2)*(exp(+q*s/(k*T))-q*s/(k*T)-1))^0.5)/c_ox)+s-V_g(i) - V_wiggle(i);
        psi_s = fsolve(F,-signq*0.1); % Initial guess, sign taken due to fsolve

    else                              % sign(N_sub) > 0
        F = @(s) V_fb-(signq*((2*eps_si*k*T*abs(N_sub))^0.5*((ni^2/(abs(N_sub))^2)*(exp(-q*s/(k*T))+q*s/(k*T)-1)+(exp(+q*s/(k*T))-q*s/(k*T)-1))^0.5)/c_ox)+s-V_g(i) - V_wiggle(i);
        psi_s = fsolve(F,-signq*0.1); % Initial guess, sign taken due to fsolve
    end


    % Solver parameters
    dx = 1e-9;             % Meshing size
    NRmax = 10;             % Maximum NR iterations
    tol = 1e-6;             % Error tolerance in potential


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
        deln = 0;
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
   rho_si_dcwiggle = q*(N_sub + p - n);
   Q2(i) = sum(rho_si_dcwiggle) * dx;
end
