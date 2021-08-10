function [psi] = Poisson_solver(V_g(i) + V_wiggle(i) , psub , nsub)
    %% Poisson Solver for a given Vg 
    % Takes in Vg as the inout and outputs psi in the substrate
    %%Calculate psi_s from V_g(i) + V_wiggle

    % psi_init= -(V_g(i) + V_wiggle - V_fb)*sign(N_sub);
    % pre=(2*eps_si*Vt*q*abs(N_sub))^0.5;
    % 
    % func = @(psi_s) -V_g(i) + V_wiggle + V_fb -(1/c_ox *-sign(V_g(i) + V_wiggle-V_fb)*pre*(exp(-psi_s/Vt)+psi_s/Vt-1+(ni^2/(abs(N_sub))^2)*(exp(psi_s/Vt)-psi_s/Vt-1))^0.5) + psi_s;
    % x0 = psi_init;
    % psi_s = fzero(func,x0);


    % Sign of substrate charges in the current mode, since we have signum(psi_s) in formula of charge
    
    %Declaring the global variables
    
    global V_fb;
    global N_sub
    global eps_si;
    global k;
    global T;
    global q;
    global ni;
    global c_ox;
    global L_sub;
    global N;
    global Vt;
    global x;
    
    if(V_g(i) + V_wiggle>=V_fb)
        signq = -1;       %Depletion for NMOS, Accumulation for PMOS
    else                  % i.e. (V_g(i) + V_wiggle<V_fb)
        signq = +1;       %Accumulation for NMOS, Depletion for PMOS
    end

    %Same charge sign for depletion and inversion

    if sign(N_sub) <= 0
        F = @(s) V_fb-(signq*((2*eps_si*k*T*abs(N_sub))^0.5*((exp(-q*s/(k*T))+q*s/(k*T)-1)+(ni^2/(abs(N_sub))^2)*(exp(+q*s/(k*T))-q*s/(k*T)-1))^0.5)/c_ox)+s-V_g(i) + V_wiggle(i);
        psi_s = fsolve(F,-signq*0.1); % Initial guess, sign taken due to fsolve

    else                              % sign(N_sub) > 0
        F = @(s) V_fb-(signq*((2*eps_si*k*T*abs(N_sub))^0.5*((ni^2/(abs(N_sub))^2)*(exp(-q*s/(k*T))+q*s/(k*T)-1)+(exp(+q*s/(k*T))-q*s/(k*T)-1))^0.5)/c_ox)+s-V_g(i) + V_wiggle(i);
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
end