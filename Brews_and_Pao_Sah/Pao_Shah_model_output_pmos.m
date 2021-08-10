% Solving Pao-Shah model for the following NMOSFET
% L = 1 um, Tox = 10 nm, n+ poly-Si gate, Na (choose) so that V_T -> 0.8V
% range for the Vg and Vd sweeps is 0-5V

% Plot the Id - Vg characteristics for different values of Vd (transfer)
% Plot the Id - Vd characteristics for different values of Vg ( Vg > Vt )

% Compare it with brews and piece-wise model 

% Re-do above for a PMOSFET, with p+ poly Si gate
% use const eff. mobility values of = 200 cm^2 / Vs (e-) and 100 (holes)

% Constants 
q = 1.6e-19;
eps_0 = 8.85e-12;
kT = 26e-3*q;

% Semiconductor
Nsub = 4e17*1e6;       % -ve p-substrate -> nmos & +ve for pmos
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

%% Plot Id - Vd characteristics for different values of Vg (output)
Vd_vec = -sign(Nsub) * (0:50e-3:2);
dpsi_s = -sign(Nsub) * 1e-3;

for Vg = -sign(Nsub) * (0.5:0.5:2)
    Id_vec = [];
    for Vd = -sign(Nsub) * (0:50e-3:2)
        
    psi_s_min = -Vg-Vd+sign(Nsub)*abs(Vfb);
    psi_s_max = Vg+Vd-sign(Nsub)*abs(Vfb);
    
    psi_svec = psi_s_min:dpsi_s:psi_s_max;

    IdVs = []; 
    del_psi = -sign(Nsub) * 10e-3;

    for V = 0:dV:Vd
        f1 = @(psi) ni^2/Na * exp(-sign(Nsub)*q*(psi - V)/kT);
        f2 = @(psi) (2*kT*Na/eps_si)^0.5 * ( abs( -sign(Nsub)*q*psi/kT + f1(psi)/Na ) ).^0.5;
        f3 = @(psi_s) Vfb + psi_s - sign(Nsub)*eps_si/Cox*f2(psi_s);
        f1byf2 = @(psi) f1(psi)./f2(psi);

        Vgs = f3(psi_svec);
        psi_s = interp1(real(Vgs), real(psi_svec), Vg);
        inn_int = integral(f1byf2, del_psi, psi_s);
        IdVs = [IdVs q*muf*W/L*inn_int];
    end

    Id = sum(IdVs)*dV;
    
    Id_vec = [Id_vec Id]; 
    end
    figure(3);
    plot(Vd_vec, abs(Id_vec*1e6));
    axes = gca;
    axes.LineWidth = 1; axes.FontSize = 14; axes.FontWeight = 'bold'; axes.Box = 'on';
    xlabel('Vd'); ylabel('Id (\mu A/\mu m)');
    lines = axes.Children;
    set(lines, 'LineWidth', 2);
    hold on;
end
figure(3);
if Nsub < 0
    lgd = legend('0.5','1','1.5','2');
    title(lgd,'V_G')
else 
    lgd = legend('-0.5','-1','-1.5','-2');
    title(lgd,'V_G')
end    
hold off;
