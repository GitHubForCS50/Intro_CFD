function [] = vcoeff()
% Purpose: To calculate the coefficients for the v equation.

% constants
global NPI NPJ 
% variables
global x x_u y y_v v p mu SP Su F_u F_v d_v relax_v Istart Iend Jstart Jend ...
    b aE aW aN aS aP
% Baffles
global num_baf_bot pos_bot len_bot num_baf_top pos_top len_top

Istart = 2;
Iend = NPI+1;
Jstart = 3;
Jend = NPJ+1;

convect();

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;        
        % Geometrical parameters: Areas of the cell faces
        AREAw = y(J) - y(J-1); % See fig. 6.4
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i);
        AREAn = AREAs;
        
        % eq. 6.11a-6.11d - the convective mass flux defined in eq. 5.8a
        % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.
        Fw = ((F_u(i,J)   + F_u(i,J-1))/2)*AREAw;
        Fe = ((F_u(i+1,J) + F_u(i+1,J-1))/2)*AREAe;
        Fs = ((F_v(I,j)   + F_v(I,j-1))/2)*AREAs;
        Fn = ((F_v(I,j)   + F_v(I,j+1))/2)*AREAn;
        
        % eq. 6.11e-6.11h - the transport by diffusion defined in eq. 5.8b
        % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition
        Dw = ((mu(I-1,J-1) + mu(I,J-1) + mu(I-1,J) + mu(I,J))/(4*(x(I) - x(I-1))))*AREAw;
        De = ((mu(I,J-1) + mu(I+1,J-1) + mu(I,J) + mu(I+1,J))/(4*(x(I+1) - x(I))))*AREAe;
        Ds =  (mu(I,J-1)/(y_v(j) - y_v(j-1)))*AREAs;
        Dn =  (mu(I,J)/(y_v(j+1) - y_v(j)))*AREAn;
        
        % The source terms
        SP(I,j) = 0.;
        Su(I,j) = 0.;
        
        % The coefficients (hybrid differencing scheme)
        aW(I,j) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,j) = max([-Fe, De - Fe/2, 0.]);
        aS(I,j) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,j) = max([-Fn, Dn - Fn/2, 0.]);
        
        % transport of v through the baffles can be switched off by setting the coefficients to zero
        % aE and aW have to go to 0 because the baffles are infinitesimaly
        % thin. However, if the baffles would be thicker, we would set Sp
        % to -LARGE for v as well. Then, we would also keep aW and aE
        % active because it creates a no-slip boundary - it takes the 0
        % velocity into account when calculating the momentum. If we set aW
        % and aE to 0, this would create a free-slip wall

        % % if (I == ceil((NPI+1)/5)-1 && j < ceil((NPJ+1)/3))     % left of baffle #1
        % %     aE(I,j) = 0;
        % % end
        % % if (I == ceil((NPI+1)/5)   && j < ceil((NPJ+1)/3))     % right of baffle #1
        % %     aW(I,j) = 0;
        % % end
        % % if (I == ceil(2*(NPI+1)/5)-1 && j > ceil(2*(NPJ+1)/3)) % left of baffle #2
        % %     aE(I,j) = 0;
        % % end
        % % if (I == ceil(2*(NPI+1)/5)   && j > ceil(2*(NPJ+1)/3)) % right of baffle #2
        % %     aW(I,j) = 0;
        % % end

        for baf = 1:num_baf_bot
            i_baf = max(3, min(NPI+1, round(pos_bot(baf) * NPI) + 1));
            j_max = round(len_bot(baf) * NPJ) + 1;
            if (I == i_baf - 1 && j <= j_max) % Left of baffle
                aE(I,j) = 0;
            end
            if (I == i_baf && j <= j_max)     % Right of baffle
                aW(I,j) = 0;
            end
        end
        
        % Top Baffles
        for baf = 1:num_baf_top
            i_baf = max(3, min(NPI+1, round(pos_top(baf) * NPI) + 1));
            j_min = NPJ + 2 - round(len_top(baf) * NPJ);
            if (I == i_baf - 1 && j >= j_min)
                aE(I,j) = 0;
            end
            if (I == i_baf && j >= j_min)
                aW(I,j) = 0;
            end
        end
        
        % eq. 8.31 without time dependent terms (see also eq. 5.14):
        aP(I,j) = aW(I,j) + aE(I,j) + aS(I,j) + aN(I,j) + Fe - Fw + Fn - Fs - SP(I,J);
        
        % Calculation of d(I,j) = d_v(I,j) defined in eq. 6.23 for use in the
        % equation for pression correction (eq. 6.32) (see subroutine pccoeff).
        d_v(I,j) = AREAs*relax_v/aP(I,j);
        
        % Putting the integrated pressure gradient into the source term b(I,j)
        % The reason is to get an equation on the generalised form
        % (eq. 7.7 ) to be solved by the TDMA algorithm.
        % note: In reality b = a0p*fiP + Su = 0.
        b(I,j) = (p(I,J-1) - p(I,J))*AREAs + Su(I,j);
        
        % Introducing relaxation by eq. 6.37 . and putting also the last
        % term on the right side into the source term b(i,J)
        aP(I,j) = aP(I,j)/relax_v;
        b(I,j)  = b(I,j) + (1 - relax_v)*aP(I,j)*v(I,j);
        
        % now we have implemented eq. 6.37 in the form of eq. 7.7
        % and the TDMA algorithm can be called to solve it. This is done
        % in the next step of the main program.        
    end
end
end

