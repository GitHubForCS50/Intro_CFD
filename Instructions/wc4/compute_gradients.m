function [] = compute_gradients()
% Purpose: Calculate velocity gradients from the old values of u, v

global NPI NPJ

global x y x_u y_v u_old v_old dudx dudy dvdx dvdy

for I = 2:NPI+1
    i = I;
    for J = 2:NPJ+1
        j = J;

        dudx(I,J) = (u_old(i+1,J)-u_old(i,J))/(x_u(i+1)-x_u(i));
        dvdy(I,J) = (v_old(I,j+1)-v_old(I,j))/(y_v(j+1)-y_v(j));
    end
end

for I = 3:NPI+1
    i = I;
    for J = 3:NPJ+1
        j = J;

        dudy(i,j) = (u_old(i,J) - u_old(i, J-1))/(y(J)-y(J-1));
        dvdx(i,j) = (v_old(I,j) - v_old(I-1, j))/(x(I)-x(I-1));
    end
end
