% Este script genera los tubos curvos, calcula sus skeletons y guarda los
% resultados en el archivo results.txt

% Data parameters
x = 0: 0.001: 10;
x0 = 5;
y0 = 0;
small_r = 10;

% Results file.
fid = fopen('results.txt', 'w');
fprintf(fid, 'Los resultados\r\n');

% Chosen radii.
r = [5, 6, 8, 15, 300];
for ri = 1: numel(r)
    fprintf(fid, 'r = %d\r\n', r(ri));

    % Get the curve.
    y = sqrt(r(ri)^2 - (x - x0).^2) + (y0 - r(ri));

    % Put the curve somewhere in the matrix.
    mymatrix = zeros(150, 150, 150);
    for i = 1: numel(x)
        mymatrix(25 + round(10 * x(i)), 25 + round(-10 * y(i)), 50) = 1;
    end
    [xx, yy, zz] = ndgrid(-small_r: small_r);
    
    % Pseudo-euclidean distances.
    dist = sqrt(cat(3, ...
        [3 2 3; 2 1 2; 3 2 3], ...
        [2 1 2; 1 0 1; 2 1 2], ...
        [3 2 3; 2 1 2; 3 2 3]));
    
    % Get the actual length of the curve.
    curve = mymatrix;
    mat2xraw(curve, sprintf('generated_models/curve_r%d.xraw', r(ri)));
    curve_idx = find(curve);
    [i, j, k] = ind2sub(size(curve), curve_idx);
    path_length = 0;
    for l = 1: numel(i)
        d = dist .* curve(i(l) - 1: i(l) + 1, j(l) - 1: j(l) + 1, k(l) - 1: k(l) + 1);
        curve(i(l) - 1: i(l) + 1, j(l) - 1: j(l) + 1, k(l) - 1: k(l) + 1) = 0;
        path_length = path_length + sum(d(:));
    end
    fprintf(fid, 'original = %f\r\n', path_length);

    % Dilate the curve.
    nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= small_r;
    dilated_lulo = imdilate(mymatrix, nhood);
    
    % Calculate each skeleton and get the length of each path.
    
    % Thinning.
    skel = thinning3D(dilated_lulo);
    mat2xraw(skel, sprintf('generated_models/curve_th_r%d.xraw', r(ri)));
    skel_idx = find(skel);
    [i, j, k] = ind2sub(size(skel), skel_idx);
    path_length = 0;
    for l = 1: numel(i)
        d = dist .* skel(i(l) - 1: i(l) + 1, j(l) - 1: j(l) + 1, k(l) - 1: k(l) + 1);
        skel(i(l) - 1: i(l) + 1, j(l) - 1: j(l) + 1, k(l) - 1: k(l) + 1) = 0;
        path_length = path_length + sum(d(:));
    end
    fprintf(fid, 'th = %f\r\n', path_length);
    
    % Hamilton-Jacobi.
    skel = thinning3D(hj_skeleton3D(dilated_lulo, -15));
    mat2xraw(skel, sprintf('generated_models/curve_hj_r%d.xraw', r(ri)));
    skel_idx = find(skel);
    [i, j, k] = ind2sub(size(skel), skel_idx);
    path_length = 0;
    for l = 1: numel(i)
        d = dist .* skel(i(l) - 1: i(l) + 1, j(l) - 1: j(l) + 1, k(l) - 1: k(l) + 1);
        skel(i(l) - 1: i(l) + 1, j(l) - 1: j(l) + 1, k(l) - 1: k(l) + 1) = 0;
        path_length = path_length + sum(d(:));
    end
    fprintf(fid, 'hj = %f\r\n', path_length);
    
    % Distance-driven.
    skel = thinning3D(dd_skeleton_v2(dilated_lulo, 4, 0.25));
    mat2xraw(skel, sprintf('generated_models/curve_dd_r%d.xraw', r(ri)));
    skel_idx = find(skel);
    [i, j, k] = ind2sub(size(skel), skel_idx);
    path_length = 0;
    for l = 1: numel(i)
        d = dist .* skel(i(l) - 1: i(l) + 1, j(l) - 1: j(l) + 1, k(l) - 1: k(l) + 1);
        skel(i(l) - 1: i(l) + 1, j(l) - 1: j(l) + 1, k(l) - 1: k(l) + 1) = 0;
        path_length = path_length + sum(d(:));
    end
    fprintf(fid, 'dd = %f\r\n', path_length);
end