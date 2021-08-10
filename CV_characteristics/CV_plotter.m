clear all

q = 1.6e-19;        % electron charge
eps_0 = 8.85e-12;   % Vacuum permitivity
VKBt = 26e-3;       % KBT/q at 300K

%% Semiconductor
Nsub = -1e17*1e6;               % Negative for p-substrate(nmos), Positive for n-substrate (pmos)
k_si = 12;                      % Dielectric Constant
ni = 1.5e10*1e6;                % Intrinsic Carrier concentration
Eg = 1.1 * q;                   % Bandgap of the semiconductor material
Lsub = 250e-9;                       % substrate thickness
eps_si = k_si*eps_0;            % permitivity
Ea_sub = 4.05*q;                % electron affinity for substrate
if Nsub < 0
    psub = abs(Nsub);  nsub = ni^2/abs(Nsub);  % can be made more accurate
elseif Nsub > 0
    psub = ni^2/Nsub; nsub = Nsub;
else
    psub = ni; nsub = ni;
end

%% Oxide
tox = 3e-9;                     % Thickness
k_ox = 4;                       % Dielectric Constant
Eg_ox = 9*q;                    % Bandgap
eps_ox = k_ox*eps_0;            % Permittivity
Ec_off = 3.1*q;                 % Conduction band offset
C_ox = eps_ox/tox;              % Oxide capacitance

%% Flatband / Metal related
phi_m = (4.05)*q;               %work function of metal
tm = 10e-9;                     % gate thickness (has no significance for metal gate)
phi_b = -sign(Nsub) * VKBt * log(abs(Nsub)/ni);                     % (Ei - Ef)/q of semiconductor
phi_sub = Ea_sub + sign(Nsub) * Eg/(2) + phi_b*q;                   %work function for substrate
Vfb =  (phi_m-phi_sub)/q;                                           % n+/nmos, p+/pmos, i/i

%% Solver
dx = 0.1e-9;                      % Meshing size
NRmax = 100;                     % Maximum NR iterations
tol = 1e-6;                     % Error tolerance in potential

%% Mesh
xox = -tox:dx:0;
xsi = dx:dx:Lsub;
x = [xox, xsi];                 % mesh, interface is origin
N = length(xsi) + 1;                  % number of points
iox = x <= 0;
isi = x > 0;
i0 = find(x == 0);

%% Newton Raphson - Matrices
A = (-2*diag(ones(1, N),0) + 1*diag(ones(1,N-1),1) + 1*diag(ones(1,N-1),-1));
A(1, :) = 0; A(1,1) = 1;
A(N, :) = 0; A(N,N) = 1;

b = zeros(N,1);
psi = zeros(N,1);

%% Inputs 
Vgs = (-5: 50e-3: 5);
Vwiggle = 10e-3;
Vright = 0;
sweep_type = 0;                 % 0 = LFCV, 1 = HFCV, 2 = Fast Sweep

Caps = [];
psi_ss = [];
caplast = C_ox;

%% Newton Raphson Iterations

Caps = [];
psi_ss = [];

% Voltage Sweep
for Vg = Vgs
    disp(Vg)
    
    % Inputs at each voltage
    V1 = Vg - Vfb;
    V2 = V1 + Vwiggle;
    Qsub = [];
    
    % Bias input------------------
    Vleft = V1;
    psi_s_now = Vleft;
    Fact = sqrt(2*eps_si*q*VKBt*abs(Nsub));
    init = psi_s_now;
    for i = 1:1000                  % USED SIGN OF Q_S to OBTAIN CONVERGENT SOLUTION
        init = psi_s_now;
        Q_s =  Fact * (exp(-psi_s_now/VKBt) + psi_s_now /VKBt-1 + ni^2 /(abs(Nsub))^2*(exp(psi_s_now/VKBt)-psi_s_now/VKBt-1))^0.5;
        if (Vleft >= 0)
           Q_s = -Q_s; 
        end
        f= -Vleft + Vfb - Q_s/C_ox + psi_s_now;
        df=-1/C_ox*(Fact)^2*0.5*(-1/VKBt*exp(-psi_s_now/VKBt)+1/VKBt+ni^2/abs(Nsub)^2*(1/VKBt*exp(psi_s_now/VKBt)-1/VKBt))/Q_s + 1;
        psi_s_now = psi_s_now - f/df;
        if(abs(psi_s_now - init) < 1e-6)
            i
            break;
        end
    end    
    disp(Q_s);
    disp(psi_s_now);
    for i = 1:NRmax
        % Generate next guess
        p = psub*exp(-psi/VKBt);        % psi(x), p(psi)
        n = nsub*exp(psi/VKBt);         
        rho = q*(Nsub + p - n);         % Nsub = Nd - Na
        b = -rho/eps_si*dx^2;           
        b(1) = psi_s_now; b(N) = 0;         % Dirichlet Boundary Condition
        f = A*psi - b;

        % Jacobian calculations
        delp = -1/VKBt*p;
        deln = 1/VKBt*n;
        delrho = q*(delp - deln);
        delb = -delrho/eps_si*dx^2;
        delb(1) = 0; delb(N) = 0;
        J = A - diag(delb);             % delb gets added to diagonal entries
        dV = -J\f;                      % Potential update across x

    if max(abs(dV)) < tol
        break;
    end

    psi = psi + dV;
    end
    Q1 =  sum(rho)*dx; 
    
    % Wiggle point
    Vleft = V2;
    psi_s_now = Vleft;
    Fact = sqrt(2*eps_si*q*VKBt*abs(Nsub));
    init = psi_s_now;
    for i = 1:1000                  % USED SIGN OF Q_S to OBTAIN CONVERGENT SOLUTION
        init = psi_s_now;
        Q_s =  Fact * (exp(-psi_s_now/VKBt) + psi_s_now /VKBt-1 + ni^2 /(abs(Nsub))^2*(exp(psi_s_now/VKBt)-psi_s_now/VKBt-1))^0.5;
        if (Vleft >= 0)
           Q_s = -Q_s; 
        end
        f= -Vleft + Vfb - Q_s/C_ox + psi_s_now;
        df=-1/C_ox*(Fact)^2*0.5*(-1/VKBt*exp(-psi_s_now/VKBt)+1/VKBt+ni^2/abs(Nsub)^2*(1/VKBt*exp(psi_s_now/VKBt)-1/VKBt))/Q_s + 1;
        psi_s_now = psi_s_now - f/df;
        if(abs(psi_s_now - init) < 1e-6)
            i
            break;
        end
    end    
    for i = 1:NRmax
        % Generate next guess
        p = psub*exp(-psi/VKBt);        % psi(x), p(psi)
%       n = nsub*exp(psi/VKBt);       % For HFCV use the n obtained for bias only
        rho = q*(Nsub + p - n);         % Nsub = Nd - Na
        b = -rho/eps_si*dx^2;           
        b(1) = psi_s_now; b(N) = 0;         % Dirichlet Boundary Condition
        f = A*psi - b;

        % Jacobian calculations
        delp = -1/VKBt*p;
        deln = 1/VKBt*n;
        delrho = q*(delp - deln);
        delb = -delrho/eps_si*dx^2;
        delb(1) = 0; delb(N) = 0;
        J = A - diag(delb);             % delb gets added to diagonal entries
        dV = -J\f;                      % Potential update across x

    if max(abs(dV)) < tol
        break;
    end

    psi = psi + dV;
    end
    Q2 =  sum(rho)*dx;
    Cap = 1/Vwiggle*(Q1 - Q2);
    if abs(Cap) >= 100*abs(caplast) || abs(Cap) <= 0.01 * abs(caplast)
        Cap = caplast;
    end
    Caps = [Caps Cap];
    caplast = Cap;
end                                     % end for V1 and V2

plot(Vgs, Caps)


































