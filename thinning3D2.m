% Our goal is programming a 3D thinning algorithm based upon
% MATLAB's morphological operations, replicating (kind of) what's
% done here: http://homepages.inf.ed.ac.uk/rbf/HIPR2/thin.htm
% input: MATLAB 3D matrix (such as read_binvox's output).
% output: MATLAB 3D matrix representing input's skeleton.
function skel = thinning3D2(I)

    % Zero-padding to avoid boundary conditions.
    pad = 2;
    skel = padarray(I, [pad pad pad]);
    dims = size(skel);
    
    % Masks as boolean functions.
    masks = {
        @(N) N(2,2,2) & ~(N(1,1,1)|N(1,2,1)|N(1,3,1)|N(1,1,2)|N(1,2,2)|N(1,3,2)|N(1,1,3)|N(1,2,3)|N(1,3,3)) & N(3,2,2) & (N(2,1,1)|N(2,2,1)|N(2,3,1)|N(3,1,1)|N(3,2,1)|N(3,3,1)|N(2,1,2)|N(2,3,2)|N(3,1,2)|N(3,3,2)|N(2,1,3)|N(2,2,3)|N(2,3,3)|N(3,1,3)|N(3,2,3)|N(3,3,3))
        @(N) N(2,2,2) & ~(N(1,1,1)|N(1,2,1)|N(1,3,1)|N(2,1,1)|N(2,2,1)|N(2,3,1)|N(3,1,1)|N(3,2,1)|N(3,3,1)) & N(2,2,3) & (N(1,1,2)|N(1,2,2)|N(1,3,2)|N(2,1,2)|N(2,3,2)|N(3,1,2)|N(3,2,2)|N(3,3,2)|N(1,1,3)|N(1,2,3)|N(1,3,3)|N(2,1,3)|N(2,3,3)|N(3,1,3)|N(3,2,3)|N(3,3,3))
        @(N) N(2,2,2) & ~(N(1,1,1)|N(1,2,1)|N(1,3,1)|N(2,1,1)|N(2,2,1)|N(2,3,1)|N(3,1,1)|N(3,2,1)|N(3,3,1)|N(1,1,2)|N(1,2,2)|N(1,3,2)|N(1,1,3)|N(1,2,3)|N(1,3,3)) & N(3,2,3) & (N(2,1,2)|N(2,3,2)|N(3,1,2)|N(3,3,2)|N(2,1,3)|N(2,3,3)|N(3,1,3)|N(3,3,3))
        @(N) N(2,2,2) & ~(N(1,1,1)|N(1,2,1)|N(1,3,1)|N(2,2,1)|N(1,2,2)) & N(3,2,2) & N(2,2,3) & ~(N(2,1,1) & N(1,1,2)) & ~(N(2,3,1) & N(1,3,2))
        @(N) N(2,2,2) & ~(N(1,1,1)|N(1,2,1)|N(2,2,1)|N(1,2,2)) & N(1,3,1) & N(3,2,2) & N(2,2,3) & ~(N(2,1,1) & N(1,1,2)) & xor(N(2,3,1), N(1,3,2))
        @(N) N(2,2,2) & ~(N(1,2,1)|N(1,3,1)|N(2,2,1)|N(1,2,2)) & N(1,1,1) & N(3,2,2) & N(2,2,3) & ~(N(2,3,1) & N(1,3,2)) & xor(N(2,1,1), N(1,1,2))
        @(N) N(2,2,2) & ~(N(1,1,1)|N(1,2,1)|N(2,2,1)|N(1,2,2)) & N(2,3,2) & N(3,2,2) & N(2,2,3) & ~(N(2,1,1) & N(1,1,2))
        @(N) N(2,2,2) & ~(N(1,2,1)|N(1,3,1)|N(2,2,1)|N(1,2,2)) & N(2,1,2) & N(3,2,2) & N(2,2,3) & ~(N(2,3,1) & N(1,3,2))
        @(N) N(2,2,2) & ~(N(1,2,1)|N(2,2,1)|N(1,2,2)) & N(1,1,1) & N(2,3,2) & N(3,2,2) & N(2,2,3) & xor(N(2,1,1), N(1,1,2))
        @(N) N(2,2,2) & ~(N(1,2,1)|N(2,2,1)|N(1,2,2)) & N(1,3,1) & N(2,1,2) & N(3,2,2) & N(2,2,3) & xor(N(2,3,1), N(1,3,2))
        @(N) N(2,2,2) & ~(N(1,1,1)|N(1,2,1)|N(1,3,1)|N(2,1,1)|N(2,2,1)|N(3,1,1)|N(3,2,1)|N(1,1,2)|N(1,2,2)|N(1,3,2)|N(1,1,3)|N(1,2,3)|N(1,3,3)) & N(3,3,2) & N(3,2,3)
        @(N) N(2,2,2) & ~(N(1,1,1)|N(1,2,1)|N(1,3,1)|N(2,2,1)|N(2,3,1)|N(3,2,1)|N(3,3,1)|N(1,1,2)|N(1,2,2)|N(1,3,2)|N(1,1,3)|N(1,2,3)|N(1,3,3)) & N(3,1,2) & N(3,2,3)
        @(N) N(2,2,2) & ~(N(1,1,1)|N(1,2,1)|N(1,3,1)|N(2,1,1)|N(2,2,1)|N(2,3,1)|N(3,1,1)|N(3,2,1)|N(3,3,1)|N(1,1,2)|N(1,2,2)|N(1,1,3)|N(1,2,3)) & N(2,3,3) & N(3,2,3)
        @(N) N(2,2,2) & ~(N(1,1,1)|N(1,2,1)|N(1,3,1)|N(2,1,1)|N(2,2,1)|N(2,3,1)|N(3,1,1)|N(3,2,1)|N(3,3,1)|N(1,2,2)|N(1,3,2)|N(1,2,3)|N(1,3,3)) & N(2,1,3) & N(3,2,3)
    };

    % 12 rotations
    mask_rotations = 12;
    
    % Thin until no voxels are removed.
    while(1)
        voxels_removed = 0;
        for rot = 1: mask_rotations
            fprintf('rot %d\n', rot);

            % Linear indices for image.
            skel_li = find(skel);
            [i, j, k] = ind2sub(dims, skel_li);
            marked_for_removal = zeros(dims);

            for m = 1: numel(masks)
                for l = 1: numel(i)

                    % This voxel's (rotated) sorroundings.
                    ri = i(l) - 1: i(l) + 1;
                    rj = j(l) - 1: j(l) + 1;
                    rk = k(l) - 1: k(l) + 1;
                    N = rotate_mask(skel(ri, rj, rk), rot);

                    % Apply the mask.
                    if masks{m}(N)
                        marked_for_removal(i(l), j(l), k(l)) = 1;
                        voxels_removed = voxels_removed + 1;
                    end
                end
            end
            skel = skel - marked_for_removal;
        end

        fprintf('removed %d\n', voxels_removed);
        if voxels_removed == 0
            break;
        end
    end
    
    % Remove padding.
    skel = skel(pad + 1: end - pad, pad + 1: end - pad, pad + 1: end - pad);
end
