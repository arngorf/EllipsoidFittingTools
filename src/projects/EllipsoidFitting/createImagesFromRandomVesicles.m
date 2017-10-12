numVesicles = 300;

w = 300;
h = 300;
d = 200;

I = zeros(w,h,d);
I(:,:,:) = 173.0/255.0;

fullSavePath = '/home/dith/EllipsoidFittingTools/src/projects/EllipsoidFitting/simed.tif'

ellip_kxs = zeros(d, 1);
ellip_kys = zeros(d, 1);
nVesInSlice = zeros(d, 1);

true_kx = 0.1;
true_ky = -0.25;

for n = 1:numVesicles

    n

    D = zeros(3);

    radii = rand(3);
    center = rand(3);

    cx = center(1) * (2*abs(true_kx)*d+w) - abs(true_kx)*d;
    cy = center(2) * (2*abs(true_ky)*d+h) - abs(true_ky)*d;
    cz = center(3) * d;

    for j = 1:3
        radii(j) = radii(j)*5.0 + 7.5;
    end

    max_rad = max(radii);

    radii = radii.^2;

    D(1,1) = 1.0 / radii(1);
    D(2,2) = 1.0 / radii(2);
    D(3,3) = 1.0 / radii(3);

    R = RandomRotationMatrix()*RandomRotationMatrix()*RandomRotationMatrix();
    H = transpose(R)*D*R;

    A = H(1,1);
    B = H(2,2);
    %C = H(3,3);

    D = H(2,1);
    E = H(3,1);
    F = H(3,2);

    kx_e = (D*F - B*E) / (A*B - D*D) + true_kx;
    ky_e = (D*E - A*F) / (A*B - D*D) + true_ky;

    for k = 1:d

        if (k < (cz - max_rad)) | (k > (cz + max_rad))
            continue;
        end

        found = false;

        for j = 1:h
            if (j < (cy + true_ky*k - max_rad)) | (j > (cy + true_ky*k + max_rad))
                continue;
            end
            for i = 1:w

                x = [i - cx - true_kx*k; j - cy - true_ky*k; k - cz];

                if abs(transpose(x)*H*x - 1) < 0.25
                    I(i,j,k) = 120.0/255.0;
                    found = true;
                end
            end
        end
        if found
            ellip_kxs(k,1) = ellip_kxs(k,1) + kx_e;
            ellip_kys(k,1) = ellip_kys(k,1) + ky_e;
            nVesInSlice(k,1) = nVesInSlice(k,1) + 1;
        end
    end
end

1

I = imgaussfilt3(I,1.5);
I = imnoise(I,'gaussian',0,0.001);

2

for i = 1:d
    if nVesInSlice(i,1) > 0
        ellip_kxs(i,1) = ellip_kxs(i,1) / nVesInSlice(i,1);
        ellip_kys(i,1) = ellip_kys(i,1) / nVesInSlice(i,1);
    else
        ellip_kxs(i,1) = 0;
        ellip_kys(i,1) = 0;
    end
    if i == 1
        imwrite(I(:,:,i), fullSavePath);
    else
        imwrite(I(:,:,i), fullSavePath, 'WriteMode', 'append');
    end
end

3

matlab_kxs = zeros(d-1,1);
matlab_kys = zeros(d-1,1);

4

for i = 2:d
    tform = register(I(:,:,i), I(:,:,i-1), 'translation');
    [x, y] = transformPointsForward(tform, 0, 0);

    matlab_kxs(i-1,1) = x;
    matlab_kys(i-1,1) = y;
end

5

gx=sprintf('%d, ', matlab_kxs);
fprintf('matlab_kxs: %s\n', gx)
gy=sprintf('%d, ', matlab_kys);
fprintf('matlab_kys: %s\n', gy)
gx=sprintf('%d, ', ellip_kxs);
fprintf('kxs: %s\n', gx)
gy=sprintf('%d, ', ellip_kys);
fprintf('kys: %s\n', gy)