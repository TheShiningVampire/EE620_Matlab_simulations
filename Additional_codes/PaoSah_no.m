

% Constants
q = 1.6e-19;
eps_0 = 8.85e-12;
Vt = 26e-3;               %thermal voltage
k = 8.314/(6.023*10^23);
T = 300;
kT = k*T;

%% Parameters

% Substrate Parameters
k_si = 12;                  %realtive permitivitty of silicon substrate
ni = 1.5e10*1e6;              %Intrinsic carrier density ( /cm^3)    
Eg = 1.1*q;                 %Band gap of silicon
Ea=4.05*q;                  %electron affinity for substrate
eps_si = k_si*eps_0;        %permitivitty of silicon substrate
chi_si = 4.05*q;

N_sub = -4*1e17*1e6;           % Substrate_doping in /m^3 
                               %-ve for p-substrate(NMOS), +ve for n-substrate(PMOS)

W = 1e-6;                   % Width
L = 1e-6;                   % Length
muf = 200*1e-4;             % effective mobility

Na = abs(N_sub);             % Doping value

% Oxide parameters
t_ox = 10e-9;
k_ox = 4;
Eg_ox = 9*q;                % band gap of oxide
eps_ox = k_ox*eps_0;        % permitivitty of SiO2
Ec_off = 3.1*q;             % conduction band offset
c_ox=eps_ox/t_ox;


% Gate parameters
phi_m = 4.1*q;                 % Metal_Workfunction in eV
t_m = 100e-9;
phi_b = -sign(N_sub)*Vt*log(abs(N_sub)/ni);         %Ei-Ef for semiconductor
phi_sub=(Ea+Eg/2+phi_b*q);                          %work function for substrate
V_fb = (phi_m-phi_sub)/q;                           %flatband voltage
V_th = V_fb -sign(N_sub)*(4*q*eps_si*abs(N_sub)*abs(phi_b))^0.5/c_ox + 2*phi_b;


%% Pao-Sah Model

Vg = 2;
Id = [];
for Vd = 0:0.1:2
    
    % Reference- Code presented during live session
    psi_s_min = -Vg-Vd-abs(V_fb);
    psi_s_max = Vg+Vd+abs(V_fb);
    dpsi_s = 10e-3;

    psi_svec = psi_s_min:dpsi_s:psi_s_max;

    IdVs = [];
    dV = 1e-3;
    del_psi = 10e-3;

    for V = 0:dV:Vd
        f1 = @(psi) ni^2/Na * exp(q*(psi - V)/kT);
        f2 = @(psi) (2*kT*Na/eps_si)^0.5 * (q*psi/kT + f1(psi)/Na).^0.5;
        f3 = @(psi_s) V_fb + psi_s +eps_si/c_ox * f2(psi_s);
        f1byf2 = @(psi) f1(psi)./f2(psi);

        Vgs = f3(psi_svec);
        psi_s = interp1(real(Vgs),real(psi_svec),Vg);

        inn_int = integral(f1byf2 , del_psi , psi_s);
        IdVs = [IdVs q*muf*W/L*inn_int];
    end
    
    Id = [Id sum(IdVs)*dV]
end
