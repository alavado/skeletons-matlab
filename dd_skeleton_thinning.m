function skel = dd_skeleton_thinning(V, step)
    
    dims = size(V);
    
    % Each row is a mask.
    face_direction_masks = [
        -1  0  0  1  0  0  2  0  0
         0 -1  0  0  1  0  0  2  0
         0  0 -1  0  0  1  0  0  2
         ];
         
    edge_direction_masks = [
        -1 -1  0  1  1  0  2  2  0
         1 -1  0 -1  1  0 -2  2  0
         0 -1 -1  0  1  1  0  2  2
         0 -1  1  0  1 -1  0  2 -2
        -1  0 -1  1  0  1  2  0  2
         1  0 -1 -1  0  1 -2  0  2
        ];
    
    % Add mirrored masks
    masks = cat(1, face_direction_masks, -face_direction_masks, ...
        edge_direction_masks, -edge_direction_masks);
    masks_count = size(masks, 1);
    
    % PSS thinning.
    if step == 1
    
    % Apply all masks until no changes are witnessed.
    while true
        voxels_removed = 0;
        for mask_index = 1: masks_count
            
            % Linear indices for current volume.
            I = find(V);
            [i, j, k] = ind2sub(dims, I);
            
            % Voxel type classification.
            type2 = zeros(dims);
            for l = 1: numel(i)
                type2(i(l), j(l), k(l)) = get_pudney_c_star(V, {i(l), j(l), k(l)}) == 1;
            end
            
            % Apply current mask.
            m = masks(mask_index, :);
            for l = 1: numel(i)
                
                % Check mask compliance first.
                if ~type2(i(l), j(l), k(l)) || ...
                   V(i(l) + m(1), j(l) + m(2), k(l) + m(3)) || ...
                   ~type2(i(l) + m(4), j(l) + m(5), k(l) + m(6)) || ...
                   V(i(l) + m(7), j(l) + m(8), k(l) + m(9))
                    continue;
                end

                % Then, we check whether current type-2 voxel
                % satisfies the 3 conditions for removal.
                p = {i(l), j(l), k(l)};

                % First condition: Exactly one 26-connected
                % component in its neighborhood (thus, it's simple).
                % Can't check only C value, because it may not be simple
                % in current volume conditions.
                if ~is_simple(V, p)
                    continue;
                end
                
                % Ranges for current point.
                ri = i(l) - 1: i(l) + 1;
                rj = j(l) - 1: j(l) + 1;
                rk = k(l) - 1: k(l) + 1;

                % Second condition: Exactly one 26-connected
                % component of type-2 voxels in its neighborhood.
                N26 = type2(ri, rj, rk);
                if get_pudney_c(N26, {2, 2, 2}) ~= 1
                    continue;
                end

                % Third condition: More than one object voxel
                % in its neighborhood, to preserve end points.
                N26 = V(ri, rj, rk);
                if sum(N26(:)) < 3
                    continue;
                end

                % Apply the mask.
                V(i(l), j(l), k(l)) = 0;
                type2(i(l), j(l), k(l)) = 0;
                voxels_removed = voxels_removed + 1;
            end
        end
        
        fprintf('%d voxels removed in thinning\n', voxels_removed);
        
        % No changes were observed after applying all masks.
        if voxels_removed == 0
            break;
        end
    end
    
    % NTCS thinning.
    else
        
    % Apply all masks until no changes are witnessed.
    while true
        voxels_removed = 0;
        for mask_index = 1: masks_count
            
            % Apply current mask.
            m = masks(mask_index, :);
            for i = 2: dims(1) - 1
                for j = 2: dims(2) - 1
                    for k = 2: dims(3) - 1
                        
                        % Check mask compliance first.
                        if V(i + m(1), j + m(2), k + m(3)) || ...
                           ~V(i, j, k) || ...
                           ~V(i + m(4), j + m(5), k + m(6)) || ...
                           V(i + m(7), j + m(8), k + m(9))
                            continue;
                        end
                        
                        % Then, we check whether current object voxel
                        % satisfies the 2 conditions for removal.
                        p = {i, j, k};

                        % First condition: not an end point.
                        N26 = V(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1);
                        if nnz(N26) < 3
                            continue;
                        end

                        % Second condition: p is simple.
                        if ~is_simple(V, p)
                            continue;
                        end

                        % Apply the mask.
                        V(i, j, k) = 0;
                        voxels_removed = voxels_removed + 1;
                    end
                end
            end
        end
        
        fprintf('%d voxels removed in thinning\n', voxels_removed);
        
        % No changes were observed after applying all masks.
        if voxels_removed == 0
            break;
        end
        
    end
    end
    
    skel = V;
end

