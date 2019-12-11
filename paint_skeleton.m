function skel = paint_skeleton(skel)
% Paint the skeleton nicely.

    % The colors.
    COLOR_JOINT = 3;
    COLOR_END = 2;
    COLOR_BORING = 1;
    
    % Silly padding is needed for boundary conditions.
    skel = padarray(skel, [1 1 1]);

    % Find the skeleton voxels.
    [i, j, k] = ind2sub(size(skel), find(skel));
    
    % Paint every end point.
    for l = 1: numel(i)
        N_count = nnz(skel(i(l) - 1: i(l) + 1, j(l) - 1: j(l) + 1, k(l) - 1: k(l) + 1));
        if N_count == 2
            skel(i(l), j(l), k(l)) = COLOR_END;
        elseif N_count > 3
            skel(i(l), j(l), k(l)) = COLOR_JOINT;
        else
            skel(i(l), j(l), k(l)) = COLOR_BORING;
        end
    end
end

