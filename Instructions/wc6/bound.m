function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN YMAX Cmu Ti SMALL h NhJ
% variables
global u v  m_in m_out y_v F_u k eps

% set inlet boundary condition
for J = NhJ+1:NPJ+2
   u(2,J)   = 0; % inlet velocity
   k(1,J)   = 1.5*(U_IN*Ti)^2; % k at inlet
   eps(1,J) = Cmu^0.75 *k(1,J).^1.5/(0.07*(YMAX-h)*0.5); % epsilon at inlet
end

% set boundary condition for upper and down wall, zero velocity by default


% begin: globcont();=======================================================
% Purpose: Calculate mass in and out of the calculation domain to correct for the continuity at outlet.
convect();
m_in = 0.;
m_out = 0.;

for J = 2:NPJ+1
    j = J;
    AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
    m_in  = m_in  + F_u(2,J)*AREAw;
    m_out = m_out + F_u(NPI+1,J)*AREAw;
end
% end: globcont()==========================================================

% Velocity gradient at outlet = zero:
% Correction factor m_in/m_out is used to satisfy global continuity
u(NPI+2,2:NPJ+1) = u(NPI+1,2:NPJ+1)*m_in/(m_out+SMALL);
v(NPI+2,2:NPJ+1) = v(NPI+1,2:NPJ+1);
k(NPI+2,2:NPJ+1) = k(NPI+1,2:NPJ+1);
eps(NPI+2,2:NPJ+1) = eps(NPI+1,2:NPJ+1);
end
