function [u_new,v_new,w_new] = omega_new(~)
    %裂变后产生的中子速度方向各项同性
    kesi_1 = rand;
    kesi_2 = rand;
    a = 2*kesi_1 -1;
    b = sqrt(1 - a^2);
    phi = 2*pi*kesi_2;
    %计算速度
    u_new = cos(phi)*b;
    v_new = sin(phi)*b;
    w_new = a;
end