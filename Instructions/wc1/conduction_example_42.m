%% Description:
% This program solves steady conduction problems with Central differences.
% Translated from C to Matlab by YTANG
% References: 1. Computational Fluid Dynamics, H.K. Versteeg and W. Malalasekera, Longman Group Ltd, 1995

clear
close all
clc

%% constants
NPI         = 200;        % number of grid cells in x-direction [-]
XMAX        = 0.5;       % width of the domain [m]
OUTER_ITER  = 100;      % number of outer iterations

%% INIT
% define geometrical variables (See fig. 6.2-6.4 in ref. 1)
Dx = XMAX/NPI;          % length of volume element

% length variable for the scalar points in the x direction
x    = zeros(1,NPI+2);  % x coordinate on pressure points [m]
x(1) = 0.;
x(2) = 0.5*Dx;
for i = 3:NPI+1
    x(i) = x(i-1) + Dx;
end
x(NPI+2) = x(NPI+1) + 0.5*Dx;

% initialise variables with zero values
aE      = zeros(1,NPI+2);
aW      = zeros(1,NPI+2);
aP      = zeros(1,NPI+2);
b       = zeros(1,NPI+2);
SP      = zeros(1,NPI+2);
Su      = zeros(1,NPI+2);

% initialise variables with initial values
T       = zeros(1,NPI+2);       % Temperature
Gamma   = ones(1,NPI+2) * 1000.0; % thermal conductivity

% Definition of first internal node for the different fi variables. See fig. 9.1.
Istart  = 2;
Iend    = NPI+1;

%% BOUNDARY
% Specify boundary conditions for a calculation
T(Istart-1) = 300;
T(Iend+1)   = 500;

%% T-EQUATION
% Calculate coefficients of T equation
for i = Istart:Iend
    % Geometrical parameters: Areas of the cell faces
    AREAw = pi*0.05^2;
    AREAe = pi*0.05^2;
    
    % The diffusion conductance D=(Gamma/Dx)*AREA defined in eq. 5.8b
    Dw = ((Gamma(i-1) + Gamma(i)  )/(2*(x(i)   - x(i-1))))*AREAw;
    De = ((Gamma(i)   + Gamma(i+1))/(2*(x(i+1) - x(i)  )))*AREAe;
    
    % The source terms
    SP(i) = 0.;
    Su(i) = 0.;
    
    % The coefficients (central differencing scheme)
    aW(i) = Dw;
    aE(i) = De;
    aP(i) = aW(i) + aE(i) - SP(i);
    
    %  Setting the source term b = Su
    b(i) = Su(i);
end

%% SOLVE
for iter = 0:OUTER_ITER
    
    % solve T equation
    % Solving from left to right
    for i = Istart:Iend
        T(i) = (aE(i)*T(i+1) + aW(i)*T(i-1) + b(i))/aP(i);
    end
    % Solving from right to left
    for i = Iend:-1:Istart
        T(i) = (aE(i)*T(i+1) + aW(i)*T(i-1) + b(i))/aP(i);
    end
    
    % print intermediate results
    fprintf('%4d T[%d] = %7.3f T[%d] = %7.3f T[%d] = %7.3f T[%d] = %7.3f T[%d] = %7.3f\n',...
        iter, 1, T(2), 2, T(3), 3, T(4), 4, T(5), NPI, T(NPI+1));
end

%% OUTPUT
% write output to .txt file
if iter == OUTER_ITER
    fp = fopen('output.txt','w');
    for i = 1:NPI+2
        fprintf(fp,'%11.4e\t%11.4e\n',x(i),T(i));
    end
    fclose(fp);
end

%% plot results from output
figure(1);
plot(x, T, 's', 'MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6]);
grid on;

