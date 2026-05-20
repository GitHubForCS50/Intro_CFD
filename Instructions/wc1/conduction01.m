%% Description:
% This program solves steady conduction problems with Central differences.
% Translated from C to Matlab by YTANG
% References: 1. Computational Fluid Dynamics, H.K. Versteeg and W. Malalasekera, Longman Group Ltd, 1995

clear
close all
clc

%% constants
NPI         = 5;        % number of grid cells in x-direction [-]
XMAX        = 0.02;       % width of the domain [m]
OUTER_ITER  = 100;      % number of outer iterations
q = 1000.e3;
A = 1;

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
Gamma   = ones(1,NPI+2) * 0.5; % thermal conductivity

% Definition of first internal node for the different fi variables. See fig. 9.1.
Istart  = 2;
Iend    = NPI+1;

%% BOUNDARY
% Specify boundary conditions for a calculation
T(Istart-1) = 100;
T(Iend+1)   = 200;

%% T-EQUATION
% Calculate coefficients of T equation
aW(Istart) = 0;
aE(Istart) = ((Gamma(Istart)+Gamma(Istart+1))/(2*(x(Istart+1)-x(Istart))))*A;
SP(Istart) = -2*Gamma(Istart)*A/Dx;
Su(Istart) = q*A*Dx + 2*Gamma(Istart)*A/Dx*T(Istart-1);
aP(Istart) = aW(Istart) + aE(Istart) - SP(Istart);
b(Istart) = Su(Istart);

for i = Istart+1:Iend-1
    % Geometrical parameters: Areas of the cell faces
    AREAw = 1;
    AREAe = 1;
    
    % The diffusion conductance D=(Gamma/Dx)*AREA defined in eq. 5.8b
    Dw = ((Gamma(i-1) + Gamma(i)  )/(2*(x(i)   - x(i-1))))*AREAw;
    De = ((Gamma(i)   + Gamma(i+1))/(2*(x(i+1) - x(i)  )))*AREAe;
    
    % The source terms
    SP(i) = 0.;
    Su(i) = q*(AREAe + AREAw)/2*Dx;
    
    % The coefficients (central differencing scheme)
    aW(i) = Dw;
    aE(i) = De;
    aP(i) = aW(i) + aE(i) - SP(i);
    
    %  Setting the source term b = Su
    b(i) = Su(i);
end

aW(Iend) = ((Gamma(Iend-1)+Gamma(Iend))/(2*(x(Iend)-x(Iend-1))))*A;
aE(Iend) = 0;
SP(Iend) = -2*Gamma(Iend)*A/Dx;
Su(Iend) = q*A*Dx + 2*Gamma(Iend)*A/Dx*T(Iend+1);
aP(Iend) = aW(Iend) + aE(Iend) - SP(Iend);
b(Iend) = Su(Iend);

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

