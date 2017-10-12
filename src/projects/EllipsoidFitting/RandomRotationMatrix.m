function [ R ] = RandomRotationMatrix()

    rv = rand(3);

    R = zeros(3);
    theta = rv(1) * pi * 2;
    phi   = rv(2) * pi * 2;
    z     = rv(3) * 2.0;

    r  = sqrt(z);
    Vx = sin(phi) * r;
    Vy = cos(phi) * r;
    Vz = sqrt(2.0 - z);

    st = sin(theta);
    ct = cos(theta);
    Sx = Vx * ct - Vy * st;
    Sy = Vx * st + Vy * ct;

    R(1,1) = Vx * Sx - ct;
    R(1,2) = Vx * Sy - st;
    R(1,3) = Vx * Vz;

    R(2,1) = Vy * Sx + st;
    R(2,2) = Vy * Sy - ct;
    R(2,3) = Vy * Vz;

    R(3,1) = Vz * Sx;
    R(3,2) = Vz * Sy;
    R(3,3) = 1.0 - z;
end