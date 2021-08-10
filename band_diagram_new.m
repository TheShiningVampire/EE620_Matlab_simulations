clear;
close all;
clc;

%% Parameters

% Constants
q = 1.6e-19;
eps_0 = 8.85e-12;
VkBT = 26e-3;               %thermal voltage

% Substrate Parameters
Nsub = -1e17*1e6;           %-ve for p-substrate(NMOS), +ve for n-substrate(PMOS), input doping density here
k_si = 12;                  %realtive permitivitty of silicon substrate
ni = 1.5e10*1e6;            %Intrinsic carrier density ( /cm^3)    
Eg = 1.1*q;                 %Band gap of silicon
L = 1e-6;
eps_si = k_si*eps_0;        %permitivitty of silicon substrate

%calculate mobile carriers (n and p) for substrate
if Nsub< 0 
    psub = abs(Nsub); nsub = ni^2/abs(Nsub);
elseif Nsub>0
    psub = ni^2/Nsub; nsub = Nsub;
else
    psub = ni; nsub = ni;
end

% Oxide parameters
tox = 10e-9;                % Thickness
k_ox = 4;
Eg_ox = 9*q;                % band gap of oxide
eps_ox = k_ox*eps_0;        % permitivitty of SiO2
Ec_off = 3.1*q;             % conduction band offset

% Metal parameters
tm = 10e-9;
phi_b = -sign(Nsub)*VkBT*log(abs(Nsub)/ni);     %Ei-Ef for semiconductor
Vfb = sign(Nsub)*Eg/(2*q) - phi_b;              %flatband voltage

% Solver
dx = 10e-9;             % Meshing size
NRmax = 10;             % Maximum NR iterations
tol = 1e-6;             % Error tolerance in potential

%% Input 
psi_s = phi_b;


%% Solver
%Mesh
x = 0:dx:L;
N = length(x);

%Newton Raphson Matrices

A = -2*diag(ones(1,N) , 0)+ 1*diag(ones(1,N-1),1)+ 1*diag(ones(1,N-1),-1);
A(1,:) = 0; A(1,1) = 1;
A(N,:) = 0; A(N,N) = 1;

b = zeros(N,1);
b(1) = psi_s; b(N) = 0;

psi = zeros(N,1);
%psi = A\b;

%Newton Raphson iterations
for i = 1:NRmax
    
    p = psub*exp(-psi/VkBT);
    n = nsub*exp(psi/VkBT);
    rho = q*(Nsub + p - n);
    b = -rho/eps_si * dx^2;
    b(1) = psi_s; b(N) = 0;
    f = A*psi - b;
    
    % Jacobian Matrix
    delp = -1/VkBT * p;
    deln = 1/VkBT * n;
    delrho = q*(delp-deln);
    delb = -delrho/eps_si * dx^2;
    delb(1) = 0; delb(N) = 0;
    J = A - diag(delb);
    dV = -J\f;
    
    if max(abs(dV))<tol
        break;
    end
    
    psi = psi+dV;
end

    




    