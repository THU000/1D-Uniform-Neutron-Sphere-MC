function [u_m,v_m,w_m] = omega_calculate(S,k)
    %计算散射后的方向
    A = 235;%铀235核
    u_minus = S(5,k);
    v_minus = S(6,k);
    w_minus = S(7,k);
    kesi_1 = rand;
    phi = 2*pi*rand;
    mu_c = 2*kesi_1 - 1;
    mu_L = (1 + A*mu_c) / sqrt(1 + A^2 + 2*A*mu_c);
    %中间变量
    a = mu_L;
    b = sqrt(1 - a^2);
    c = cos(phi);
    d = sin(phi);
    %计算散射后方向
    if 1 - u_minus^2 < 1E-3
        u_m = b*c;
        v_m = b*d;
        w_m = a;
    else
        u_m = a*u_minus + (b*c*w_minus*u_minus - b*d*v_minus)/sqrt(1 - u_minus^2);
        v_m = a*v_minus + (b*c*w_minus*v_minus + b*d*u_minus)/sqrt(1 - u_minus^2);
        w_m = a*w_minus - b*c*sqrt(1 - u_minus^2);
    end
end