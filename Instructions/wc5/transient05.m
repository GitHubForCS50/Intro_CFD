%% Solves: Unsteady, compressible convection-diffusion problems.
% Description:
% This program solves unsteady convection-diffusion problems
% using the transient simple algorithm described in ch. 8.7.1 in "Computational
% Fluid Dynamics" by H.K. Versteeg and W. Malalasekera. Symbols and
% variables follow exactly the notations in this reference, and all
% equations cited are from this reference unless mentioned otherwise.

% Converted from C to Matlab by YTANG
% References: 1. Computational Fluid Dynamics, H.K. Versteeg and W. Malalasekera, Longman Group Ltd, 1995

clear
close all
clc
%% declare all variables and contants
% variables
global x x_u y y_v u v pc T rho mu Gamma b SMAX SAVG aP aE aW aN aS eps k...
    u_old v_old pc_old T_old Dt eps_old k_old uplus yplus yplus1 yplus2 f f_old time
% constants
global NPI NPJ XMAX YMAX LARGE U_IN SMALL Cmu sigmak sigmaeps C1eps C2eps kappa ERough Ti Sc Sct

NPI        = 80;        % number of grid cells in x-direction [-]
NPJ        = 20;        % number of grid cells in y-direction [-]
XMAX       = 0.1;       % width of the domain [m]
YMAX       = 0.025;     % height of the domain [m]
MAX_ITER   = 10;       % maximum number of outer iterations [-]
U_ITER     = 1;         % number of Newton iterations for u equation [-]
V_ITER     = 1;         % number of Newton iterations for v equation [-]
PC_ITER    = 30;        % number of Newton iterations for pc equation [-]
T_ITER     = 1;         % number of Newton iterations for T equation [-]
EPS_ITER   = 1;         % number of Newton iterations for Eps equation [-]
K_ITER     = 1;         % number of Newton iterations for K equation [-]
F_ITER     = 1;         % number of Newton iterations for f equation [-]
SMAXneeded = 1E-8;      % maximum accepted error in mass balance [kg/s]
SAVGneeded = 1E-9;      % maximum accepted average error in mass balance [kg/s]
LARGE      = 1E30;      % arbitrary very large value [-]
SMALL      = 1E-30;     % arbitrary very small value [-]
P_ATM      = 101000.;   % athmospheric pressure [Pa]
U_IN       = 0.2;       % in flow velocity [m/s]

s          = 4.0;       % Stoichiometric ratio of chemical reaction [-]
mox_0      = 0.233;     % Mass fraction oxygen in inlet air stream [-]
mfu_1      = 1.0;       % Mass fraction fuel in inlet fuelstream [-]
fst        = mox_0/(s*mfu_1+mox_0); % Stoichiometric mixture fraction
Sc         = 0.1181;    % Schmidt number
Sct        = 0.5;       % turbulent Schmidt number

Cmu        = 0.09;
sigmak     = 1.;
sigmaeps   = 1.3;
C1eps      = 1.44;
C2eps      = 1.92;
kappa      = 0.4187;
ERough     = 9.793;
Ti         = 0.04;

Dt         = 0.001;
TOTAL_TIME = 1.;

%% start main function here
init(); % initialization
bound(); % apply boundary conditions
timestep = 0;

for time = Dt:Dt:TOTAL_TIME
    iter = 0;
    timestep = timestep+1;
    
    % outer iteration loop
    while iter < MAX_ITER && SMAX > SMAXneeded && SAVG > SAVGneeded
        
        derivatives();
        ucoeff();
        for iter_u = 1:U_ITER
            u = solve(u, b, aE, aW, aN, aS, aP);
        end
        
        vcoeff();
        for iter_v = 1:V_ITER
            v = solve(v, b, aE, aW, aN, aS, aP);
        end
        
        bound();
        
        pccoeff();
        for iter_pc = 1:PC_ITER
            pc = solve(pc, b, aE, aW, aN, aS, aP);
        end
        
        velcorr(); % Correct pressure and velocity
        
        kcoeff();
        for iter_k = 1:K_ITER
            k = solve(k, b, aE, aW, aN, aS, aP);
        end
        
        epscoeff();
        for iter_eps = 1:EPS_ITER
            eps = solve(eps, b, aE, aW, aN, aS, aP);
        end
        
        fcoeff();
        for iter_f = 1:F_ITER
            f = solve(f, b, aE, aW, aN, aS, aP);
        end


        Tcoeff();
        for iter_T = 1:T_ITER
            T = solve(T, b, aE, aW, aN, aS, aP);
        end
%         
        viscosity();
        bound();
        
        % begin:storeresults()=============================================
        % Store data at current time level in arrays for "old" data
        % To newly calculated variables are stored in the arrays ...
        % for old variables, which can be used in the next timestep       
        u_old(3:NPI+1,2:NPJ+1)   = u(3:NPI+1,2:NPJ+1);       
        v_old(2:NPI+1,3:NPJ+1)   = v(2:NPI+1,3:NPJ+1);        
        pc_old(2:NPI+1,2:NPJ+1)  = pc(2:NPI+1,2:NPJ+1);        
        T_old(2:NPI+1,2:NPJ+1)   = T(2:NPI+1,2:NPJ+1);        
        eps_old(2:NPI+1,2:NPJ+1) = eps(2:NPI+1,2:NPJ+1);
        k_old(2:NPI+1,2:NPJ+1)   = k(2:NPI+1,2:NPJ+1);
        f_old(2:NPI+1,2:NPJ+1)   = f(2:NPI+1,2:NPJ+1);        
        % end: storeresults()==============================================
        
        % increase iteration number
        iter = iter +1;        
    end % end of while loop (outer interation)
    
    % begin: printConv(time,iter)=========================================
    % print convergence to the screen
    if time == Dt
        fprintf ('Iter\t Time\t u\t v\t T\t SMAX\t SAVG\n');
    end    
    fprintf ("%4d %10.3e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\n",iter,...
        time,u(ceil(3*(NPI+1)/10),ceil(2*(NPJ+1)/5)),v(ceil(3*(NPI+1)/10),ceil(2*(NPJ+1)/5)),...
        T(ceil(3*(NPI+1)/10),ceil(2*(NPJ+1)/5)), f(ceil(3*(NPI+1)/10),ceil(2*(NPJ+1)/5)), SMAX, SAVG);
    % end: printConv(time, iter)===========================================
    if timestep==1 || mod(timestep,10)==0
        createGIF(f);
    end
    
    % reset SMAX and SAVG
    SMAX = LARGE;
    SAVG = LARGE;   
end % end of calculation

%% determine mass fraction of three components according to SCRS model 
% mox = zeros(NPI+2,NPJ+2);
% mfu = zeros(NPI+2,NPJ+2);
% mpr = zeros(NPI+2,NPJ+2);
% mni = zeros(NPI+2,NPJ+2);
%   for I = 1:NPI+2
%                 for J = 1:NPJ+2
%                     if f(I,J) < 1 && f(I,J) >= fst
%                         mox(I,J) = 0;
%                         mfu(I,J) = (f(I,J)-fst)/(1-fst)*mfu_1;
%                     elseif f(I,J) > 0 && f(I,J) < fst
%                         mfu(I,J) = 0;
%                         mox(I,J) = ((fst - f(I,J))/fst)*mox_0;
%                     end
%                     
%                     mni(I,J) = mni_0*(1-f(I,J));
%                     mpr(I,J) = 1 - (mox(I,J) + mfu(I,J) + mni(I,J)); 
%                 end
%   end
%% begin: output()
% Print all results in output.txt
fp = fopen('output.txt','w');
for I = 1:NPI+1
    i = I;
    for J = 2:NPJ+1
        j = J;
        ugrid = 0.5*(u(i,J)+u(i+1,J)); % interpolated horizontal velocity
        vgrid = 0.5*(v(I,j)+v(I,j+1)); % interpolated vertical velocity
        fprintf(fp,'%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\n',...
            x(I), y(J), ugrid, vgrid, pc(I,J), T(I,J), rho(I,J), mu(I,J), Gamma(I,J), ...
            k(I,J), eps(I,J), uplus(I,J), yplus(I,J), yplus1(I,J), yplus2(I,J), f(I,J));
    end
    fprintf(fp, '\n');
end
fclose(fp);

% Plot vorticity in vort.txt
vort = fopen('vort.txt','w');
for I = 2:NPI+1
    i = I;
    for J = 2:NPJ+1
        j = J;
        vorticity = (u(i,J) - u(i,J-1)) / (y(J) - y(J-1)) - (v(I,j) - v(I-1,j)) / (x(I) - x(I-1));
        fprintf(vort, '%11.5e\t%11.5e\t%11.5e\n',x(I), y(J), vorticity);
    end
    fprintf(vort,'\n');
end
fclose(vort);

% Plot streamlines in str.txt
str = fopen('str.txt', 'w');
for I = 1:NPI+1
    i = I;
    for J = 1:NPJ+1
        j = J;
        stream = -0.5*(v(I+1,j)+v(I,j))*(x(I+1)-x(I))+0.5*(u(i,J+1)+u(i,J))*(y(J+1)-y(J));
        fprintf(str, '%11.5e\t%11.5e\t%11.5e\n',x(I), y(J), stream);
    end
    fprintf(str,'\n');
end
fclose(str);

% Plot horizontal velocity components in velu.txt
velu = fopen('velu.txt','w');
for I = 2:NPI+2
    i = I;
    for J = 1:NPJ+2
        fprintf(velu, '%11.5e\t%11.5e\t%11.5e\n',x_u(i), y(J), u(i,J));
    end
    fprintf(velu, '\n');
end
fclose(velu);

% Plot vertical velocity components in velv.txt
velv = fopen('velv.txt','w');
for I = 1:NPI+2
    for J = 2:NPJ+2
        j = J;
        fprintf(velv, '%11.5e\t%11.5e\t%11.5e\n',x(I), y_v(j), v(I,j));
    end
    fprintf(velv,'\n');
end
fclose(velv);
% end output()

%% plot vector map
[X,Y]=meshgrid(y_v, x_u);
quiver(X,Y,u,v,0.05);
%axis equal;
figure;
contour(x,y,f');
