function A = multitiff2mat(fname)
    info = imfinfo(fname);
    num_images = numel(info);
    z_scale = 1;
    A = zeros(info(1).Height, info(1).Width, z_scale * num_images);
    for k = 1: z_scale: z_scale * num_images
        slice = imread(fname, floor(1 + k / (z_scale + 1)));
        for z = k: k + z_scale - 1
            A(:, :, z) = slice;
        end
    end
    
    % Invert maybe.
    if nnz(A) > numel(A) / 2
        A = ~A;
    end
end
