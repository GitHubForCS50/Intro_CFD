function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN YMAX
% variables
global y u v T m_in m_out y_v F_u

% Fixed temperature/velocity at the upper and lower wall
for J = 1:NPJ+2
    %       u(2,1:NPJ+2) = U_IN; % inlet
    u(2,J) = U_IN*1.5*(1.-(2.*(y(J)-YMAX/2.)/YMAX)^2); % inlet
end

% Temperature at the walls in Kelvin
T(1:NPI+2,1) = 273.; % bottom wall
T(1:NPI+2,NPJ+2) = 273.; % top wall

% begin: globcont();=======================================================
% Purpose: Calculate mass in and out of the calculation domain to correct for the continuity at outlet.
convect();

m_in  = 0.;
m_out = 0.;

for J = 2:NPJ+1
    j = J;
    AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
    m_in  = m_in  + F_u(2,J)*AREAw;
    m_out = m_out + F_u(NPI+1,J)*AREAw;
end
% end: globcont()==========================================================

% Velocity and temperature gradient at outlet = zero:
% Correction factor m_in/m_out is used to satisfy global continuity
u(NPI+2,2:NPJ+1) = u(NPI+1,2:NPJ+1)*m_in/m_out;
v(NPI+2,2:NPJ+1) = v(NPI+1,2:NPJ+1);
T(NPI+2,1:NPJ+2) = T(NPI+1,1:NPJ+2);

end
