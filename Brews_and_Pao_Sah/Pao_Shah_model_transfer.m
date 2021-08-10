% The following code is the Pao-Shah model for the following NMOS
% L = 1 um, Tox = 10 nm, n+ poly-Si gate, Na (choose) so that V_T -> 0.8V

% const eff. mobility values of = 200 cm^2 / Vs (e-) and 100 (holes)

% Constants 
q = 1.6e-19;
eps_0 = 8.85e-12;
kT = 26e-3*q;

% Semiconductor
Nsub = -4e17*1e6;       % -ve p-substrate -> nmos
k_si = 12;
ni = 1.5e10*1e6;
Eg = 1.1*q;
eps_si = k_si*eps_0;
chi_si = 4.05*q;
Na = abs(Nsub);
L = 1e-6;
W = 1e-6;
if Nsub < 0
    muf = 200e-4;
else 
    muf = 100e-4;
end

%Oxide
tox = 10e-9;
k_ox = 4;
eps_ox = k_ox*eps_0;
Cox = eps_ox/tox;

% Metal 
if Nsub < 0
    phi_m = chi_si/q;
else 
    phi_m = chi_si/q + Eg/q;
end

phi_b = -sign(Nsub) * kT/q * log(abs(Nsub)/ni);
phi_s = chi_si/q + Eg/(2*q) + phi_b;
Vfb = phi_m - phi_s;

%% Plot the Id - Vg characteristics for different values of Vd (transfer)
% Vd = 0.1 V (Lin)
% Vd = 1 V (Sat)

Id_vec = [];
Vg_vec = 0:50e-3:2;
Vd = 0.1;
dpsi_s = 1e-3;
dV = 1e-3;
for Vg = Vg_vec
    
    psi_s_min = -Vg-Vd-abs(Vfb);
    psi_s_max = Vg+Vd+abs(Vfb);
    
    psi_svec = psi_s_min:dpsi_s:psi_s_max;

    IdVs = []; 
    del_psi = 10e-3;

    for V = 0:dV:Vd
        f1 = @(psi) ni^2/Na * exp(q*(psi - V)/kT);
        f2 = @(psi) (2*kT*Na/eps_si)^0.5 * ( q*psi/kT + f1(psi)/Na ).^0.5;
        f3 = @(psi_s) Vfb + psi_s + eps_si/Cox*f2(psi_s);
        f1byf2 = @(psi) f1(psi)./f2(psi);

        Vgs = f3(psi_svec);
        psi_s = interp1(real(Vgs), real(psi_svec), Vg);
        inn_int = integral(f1byf2, del_psi, psi_s);
        IdVs = [IdVs q*muf*W/L*inn_int];
    end

    Id = sum(IdVs)*dV;
    
    Id_vec = [Id_vec Id];
end

figure(2);
yyaxis left;
semilogy(Vg_vec, Id_vec*1e6);
axes = gca;
axes.LineWidth = 1; axes.FontSize = 14; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('Vg'); ylabel('Id (\mu A/\mu m)');
lines = axes.Children;
set(lines, 'LineWidth', 2);
legend('Lin (0.1 V_D)');

yyaxis right; 
plot(Vg_vec, Id_vec*1e6);
axes = gca;
axes.LineWidth = 1; axes.FontSize = 14; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('Vg'); ylabel('Id (\mu A/\mu m)');
lines = axes.Children;
set(lines, 'LineWidth', 2);

hold on;

Id_vec = [];
Vg_vec = 0:50e-3:2;
Vd = 1;
dpsi_s = 1e-3;
dV = 1e-3;

for Vg = Vg_vec
    
    psi_s_min = -Vg-Vd-abs(Vfb);
    psi_s_max = Vg+Vd+abs(Vfb);
    
    psi_svec = psi_s_min:dpsi_s:psi_s_max;

    IdVs = []; 
    del_psi = 10e-3;

    for V = 0:dV:Vd
        f1 = @(psi) ni^2/Na * exp(q*(psi - V)/kT);
        f2 = @(psi) (2*kT*Na/eps_si)^0.5 * ( q*psi/kT + f1(psi)/Na ).^0.5;
        f3 = @(psi_s) Vfb + psi_s + eps_si/Cox*f2(psi_s);
        f1byf2 = @(psi) f1(psi)./f2(psi);

        Vgs = f3(psi_svec);
        psi_s = interp1(real(Vgs), real(psi_svec), Vg);
        inn_int = integral(f1byf2, del_psi, psi_s);
        IdVs = [IdVs q*muf*W/L*inn_int];
    end

    Id = sum(IdVs)*dV;
    
    Id_vec = [Id_vec Id];
end

figure(2);
yyaxis left;
semilogy(Vg_vec, Id_vec*1e6,'--');
axes = gca;
axes.LineWidth = 1; axes.FontSize = 14; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('Vg'); ylabel('Id (\mu A/\mu m)');
lines = axes.Children;
set(lines, 'LineWidth', 2);
legend('Sat (1 V_D)');

yyaxis right; 
plot(Vg_vec, Id_vec*1e6,'--');
axes = gca;
axes.LineWidth = 1; axes.FontSize = 14; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('Vg'); ylabel('Id (\mu A/\mu m)');
lines = axes.Children;
set(lines, 'LineWidth', 2);

hold off; 










