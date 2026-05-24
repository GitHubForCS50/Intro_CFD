function [] = fcoeff()
% Purpose: To calculate the coefficients for the f equation.

% constants
global NPI NPJ Dt Sc Sct
% variables
global x x_u y y_v SP Su F_u F_v  rho Istart Iend ...
    Jstart Jend b aE aW aN aS aP f f_old relax_f mut mu

Istart = 2;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;

convect();

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;
        % Geometrical parameters: Areas of the cell faces
        AREAw = y_v(j+1) - y_v(j); % = A(i,J) See fig. 6.2 or fig. 6.5
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i); % = A(I,j)
        AREAn = AREAs;
        
        % The convective mass flux defined in eq. 5.8a
        % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.    
        Fw = F_u(i,J)*AREAw;
        Fe = F_u(i+1,J)*AREAe;
        Fs = F_v(I,j)*AREAs;
        Fn = F_v(I,j+1)*AREAn;
        
        % The transport by diffusion defined in eq. 5.8b
        % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition
        % Gamma, at the interface is calculated with the use of a harmonic mean.
        Dw = (mu(I-1,J)/Sc+mut(I-1,J)/Sct)*(mu(I,J)/Sc+mut(I,J)/Sct)/...
            ((mu(I-1,J)/Sc+mut(I-1,J)/Sct)*(x(I)-x_u(i)) + ...
            (mu(I,J)/Sc+mut(I,J)/Sct)*(x_u(i)-x(I-1)))*AREAw;
        De = (mu(I+1,J)/Sc+mut(I+1,J)/Sct)*(mu(I,J)/Sc+mut(I,J)/Sct)/...
            ((mu(I+1,J)/Sc+mut(I+1,J)/Sct)*(x_u(i+1)-x(I)) + ...
            (mu(I,J)/Sc+mut(I,J)/Sct)*(x(I+1)-x_u(i+1)))*AREAe;
        Ds = (mu(I,J-1)/Sc+mut(I,J-1)/Sct)*(mu(I,J)/Sc+mut(I,J)/Sct)/...
            ((mu(I,J-1)/Sc+mut(I,J-1)/Sct)*(y(J)-y_v(j)) + ...
            (mu(I,J)/Sc+mut(I,J)/Sct)*(y_v(j)-y(J-1)))*AREAs;
        Dn = (mu(I,J+1)/Sc+mut(I,J+1)/Sct)*(mu(I,J)/Sc+mut(I,J)/Sct)/...
            ((mu(I,J)/Sc+mut(I,J)/Sct)*(y(J+1)-y_v(j+1)) + ...
            (mu(I,J+1)/Sc+mut(I,J+1)/Sct)*(y_v(j+1)-y(J)))*AREAn;
       
        % The source terms
        SP(I,J) = 0.;
        Su(I,J) = 0.;
        
        % The coefficients (hybrid differencing scheme)
        aW(I,j) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,j) = max([-Fe, De - Fe/2, 0.]);
        aS(I,j) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,j) = max([-Fn, Dn - Fn/2, 0.]);
        aPold   = rho(I,J)*AREAe*AREAn/Dt;
        
        % eq. 8.31 without time dependent terms (see also eq. 5.14):
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J) + aPold;
        
        % Setting the source term equal to b        
        b(I,J) = Su(I,J) + aPold*f_old(I,J);
        
        % Introducing relaxation by eq. 6.36 . and putting also the last
        % term on the right side into the source term b(i,J)       
        aP(I,J) = aP(I,J)/relax_f;
        b(I,J)  = b(I,J) + (1.0 - relax_f)*aP(I,J)*f(I,J);
         
        % now the TDMA algorithm can be called to solve the equation.
        % This is done in the next step of the main program.        
    end
end
end

