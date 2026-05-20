%% Solves: Steady, compressible convection-diffusion problems.
% Description:
% This program solves steady convection-diffusion problems using the simple
% algorithm described in ch. 6.4 in "Computational Fluid Dynamics" by H.K. Versteeg 
% and W. Malalasekera. Symbols and variables follow exactly the notations in this reference, 
% and all equations cited are from this reference unless mentioned otherwise.

% Converted from C to Matlab by YTANG
% References: 1. Computational Fluid Dynamics, H.K. Versteeg and W. Malalasekera, Longman Group Ltd, 1995

clear
close all
clc

%% constants
NPI         = 14;      % number of grid cells in x-direction [-]
NPJ         = 14;      % number of grid cells in y-direction [-]
XMAX        = 1.;       % width of the domain [m]
YMAX        = 1.;       % length of the domain [m]
T_ITER      = 100;      % maximum iteration for solving T equation
OUTER_ITER  = 1000;     % maximum outer iteration
LARGE       = 1E30;
SMALL       = 1E-30;
U_IN        = 1;      % inlet velocity [m/s]

methods = ["CDS","UDS","HDS"];
method = methods(1);

%% INIT
% initialize all variables with zero values
x   = zeros(1,NPI+2);       % x-dir coordinates of scalars such as T, P
x_u = zeros(1,NPI+2);       % x-dir coordinates of vectors (velocity component u)
y   = zeros(1,NPJ+2);       % y-dir coordinates of scalars such as T, P
y_v = zeros(1,NPJ+2);       % y-dir coordinates of vectors (velocity component v)
u   = zeros(NPI+2,NPJ+2);   % x-dir velocity component
v   = zeros(NPI+2,NPJ+2);   % y-dir velocity component

T   = zeros(NPI+2,NPJ+2);   % temperature [K]
rho = ones(NPI+2,NPJ+2);   % density [kg/m3]
mu  = zeros(NPI+2,NPJ+2);   % dynamic viscosity [Pa*s]
Gamma = zeros(NPI+2,NPJ+2); % thermal conductivity [W/(m*K)]
Cp  = zeros(NPI+2,NPJ+2);   % specific heat capacity [J/(K*kg)]

aP  = zeros(NPI+2,NPJ+2);
aE  = zeros(NPI+2,NPJ+2);
aW  = zeros(NPI+2,NPJ+2);
aN  = zeros(NPI+2,NPJ+2);
aS  = zeros(NPI+2,NPJ+2);
b   = zeros(NPI+2,NPJ+2);

SP  = zeros(NPI+2,NPJ+2);
Su  = zeros(NPI+2,NPJ+2);

F_u = zeros(NPI+2,NPJ+2);
F_v = zeros(NPI+2,NPJ+2);

% making grid
% Length of volume element
Dx = XMAX/NPI;
Dy = YMAX/NPJ;

% Length variable for the scalar points in the x direction
x(1) = 0.;
x(2) = 0.5*Dx;
for i = 3:NPI+1
    x(i) = x(i-1) + Dx;
end
x(NPI+2) = x(NPI+1) + 0.5*Dx;

% Length variable for the scalar points in the y direction
y(1) = 0.;
y(2) = 0.5*Dy;
for j = 3:NPJ+1
    y(j) = y(j-1) + Dy;
end
y(NPJ+2) = y(NPJ+1) + 0.5*Dy;

% Length variable for the velocity components u in the x direction
x_u(1) = 0.;
x_u(2) = 0.;
for i = 3:NPI+2
    x_u(i) = x_u(i-1) + Dx;
end

% Length variable for the velocity components v in the y direction 
y_v(1) = 0.;
y_v(2) = 0.;
for j = 3:NPJ+2
    y_v(j) = y_v(j-1) + Dy;
end

% Initialising all variables
omega     = 0.9;   % Over-relaxation factor for solver
u(:,:)     = 1.;   % Velocity in x-direction
v(:,:)     = 4.;   % Velocity in y-direction
T(:,:)     = 0.;   % Temperature
rho(:,:)   = 1.;   % Density
mu(:,:)    = 0.;   % Turbulent viscosity ( = 10*mu_laminar)
Gamma(:,:) = 1.;   % Thermal conductivity
Cp(:,:)    = 1.;   % J/(K*kg) Heat capacity - assumed constant for this problem

% Definition of first internal node for the different fi variables. See fig. 9.1.
Istart = 2;
Iend   = NPI+1;
Jstart = 2;
Jend   = NPJ+1;

%% BOUNDARY
% set boundary conditions
% Fixed temperature at the upper and lower wall
T(:,1)     = 100.; % Temperature in Kelvin
T(:,NPJ+2) = 0.;   % Temperature in Kelvin

% Fixed temperature at the left and right wall
T(1,:)     = 100.; % Temperature in Kelvin
T(NPI+2,:) = 0.;   % Temperature in Kelvin

%% SOLVE
% outer iteration loop
for iter = 1:OUTER_ITER
    
    % CONVECTIVE MASS FLUX
    % Purpose: To calculate the convective mass flux component pr. unit area defined in eq. 5.7
    for I = 2:NPI+2
        i = I;
        for J = 2:NPJ+2
            j = J;
            F_u(i,J) = (rho(I-1,J)*(x(I)-x_u(i)) + rho(I,J)*(x_u(i)- x(I-1)))*u(i,J)/(x(I)-x(I-1)); % = F(i,j)
            F_v(I,j) = (rho(I,J-1)*(y(J)-y_v(j)) + rho(I,J)*(y_v(j)- y(J-1)))*v(I,j)/(y(J)-y(J-1)); 
        end
    end
    
    % T-EQUATION
    % Purpose: To calculate the coefficients for the T equation.
    for I = Istart:Iend
        i = I;
        for J = Jstart:Jend
            j = J;
            % Geometrical parameters: Areas of the cell faces
            AREAw = y_v(j+1) - y_v(j); % = A(i,j) See fig. 6.2 or fig. 6.5
            AREAe = AREAw;
            AREAs = x_u(i+1) - x_u(i); % = A(i,j)
            AREAn = AREAs;
            
            % The convective mass flux defined in eq. 5.8a
            % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.
            Fw = F_u(i,J)  *Cp(I,J)*AREAw;
            Fe = F_u(i+1,J)*Cp(I,J)*AREAe;
            Fs = F_v(I,j)  *Cp(I,J)*AREAs;
            Fn = F_v(I,j+1)*Cp(I,J)*AREAn;
            
            % The transport by diffusion defined in eq. 5.8b
            % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition
            Dw = 0.5*(Gamma(I-1,J) + Gamma(I,J))/(x(I) - x(I-1))*AREAw;
            De = 0.5*(Gamma(I,J) + Gamma(I+1,J))/(x(I+1) - x(I))*AREAe;
            Ds = 0.5*(Gamma(I,J-1) + Gamma(I,J))/(y(J) - y(J-1))*AREAs;
            Dn = 0.5*(Gamma(I,J) + Gamma(I,J+1))/(y(J+1) - y(J))*AREAn;
            
            % The source terms
            SP(I,J) = -2.*AREAw*AREAs; % "-b*T" part of the source term
            Su(I,J) = 10.*AREAw*AREAs; % "a" part of the source term
            
            % The coefficients (central differencing sheme)
            if method == "CDS"
                aW(I,J) = Dw + Fw/2;
                aE(I,J) = De - Fe/2;
                aS(I,J) = Ds + Fs/2;
                aN(I,J) = Dn - Fn/2;
            elseif method == "UDS"
                aW(I,J) = Dw + max(Fw,0);
                aE(I,J) = De + max(-Fe,0);
                aS(I,J) = Ds + max(Fs,0);
                aN(I,J) = Dn + max(-Fn,0);
            elseif method == "HDS"
                aW(I,J) = max([Fw,(Dw+Fw/2),0]);
                aE(I,J) = max([-Fe,(De-Fe/2),0]);
                aS(I,J) = max([Fs,(Ds+Fs/2),0]);
                aN(I,J) = max([-Fn,(Dn-Fn/2),0]);
            else
                fprintf("Method not defined!");
                return
            end

            
            % eq. 8.31 without time dependent terms (see also eq. 5.14):
            aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J);
            
            % Setting the source term equal to b
            b(I,J) = Su(I,J);                                 
        end
    end
    
    % now the Thomas algorithm can be called to solve the equation.
    % inner iteration loop
    for iter_T = 1:T_ITER
        % solveSOR
        % Purpose: To solve the algebraic equation 7.7. using Successive Over-relaxation
        for I = Istart:Iend
            for J = Jstart:Jend
                T(I,J) = ( aE(I,J)*T(I+1,J)  + aW(I,J)*T(I-1,J) + aN(I,J)*T(I,J+1) + aS(I,J)*T(I,J-1)...
                    + b(I,J))* omega /aP(I,J) - (omega - 1)*T(I,J);
            end
        end
    end
    
    % print intermediate results
    fprintf('%3d T(%d,%d)= %7.3f T(%d,%d)= %7.3f T(%d,%d)= %7.3f T(%d,%d)= %7.3f\n',iter, ...
        round((NPI+2)/4+1), round((NPI+2)/4+1), T(round((NPI+2)/4+1),round((NPI+2)/4+1)),...   % point A
        round(3*(NPI+2)/4), round((NPI+2)/4+1), T(round(3*(NPI+2)/4),round((NPI+2)/4+1)),... % point B
        round((NPI+2)/4+1), round(3*(NPI+2)/4), T(round((NPI+2)/4+1),round(3*(NPI+2)/4)),... % point C
        round(3*(NPI+2)/4), round(3*(NPI+2)/4), T(round(3*(NPI+2)/4),round(3*(NPI+2)/4))); % point D   
end

%% OUTPUT
% write output to txt file
fp = fopen('output.txt','w');
for I = 1:NPI
    for J = 1:NPJ
        fprintf(fp, '%e\t%e\t%e\n',x(I),y(J),T(I,J));
    end
    fprintf(fp,'\n');
end
fclose(fp);

% visulize the temperature profile
pcolor(x,y,T');


