% Default theta1 = 4
% Default theta2 = 0.25
function [skel, anchors_ntcs] = dd_skeleton_v2(V, theta1, theta2)
    
    % Local convexity sharpness threshold.
    TAU = 7;
    
    % Zero padding to avoid boundary conditions.
    pad = 2;
    V = padarray(double(V), [pad pad pad]);
    dims = size(V);

    % <3, 4, 5> distance transform.
    wf = 3;
    we = 4;
    wv = 5;
    dt = weighted_dist(V, [wf we wv]);
    
    % Neighborhood weights for <3,4,5> DT.
    weights(:, :, 1) = [wv we wv; we wf we; wv we wv];
    weights(:, :, 2) = [we wf we; wf Inf wf; we wf we];
    weights(:, :, 3) = [wv we wv; we wf we; wv we wv];
    
    % 18-connectedness mask.
    mask_18 = cat(3, [0 1 0; 1 1 1; 0 1 0], ...
                     [1 1 1; 1 1 1; 1 1 1], ...
                     [0 1 0; 1 1 1; 0 1 0]);
    
    % Center of maximal balls detection (CMBs).
    disp('Detecting maximal balls...');
    cmb = zeros(dims);
    hr_cmb = zeros(dims);
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~V(i, j, k)
                    continue;
                end
                
                % This voxel's 26-neighborhood.
                N = dt(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1);
                
                % Bordering voxels should not be considered, because
                % certain configurations might yield stupid results. Talk
                % about this in the chapter or something. The configuration
                % in question is a stair with a plain side.
                %if nnz(N) ~= 27
                %    continue;
                %end
                
                % Difference between neighborhood and this voxel.
                dif = N - dt(i, j, k);
                
                % Highly relevant CMBs.
                N_hr_cmb = dif < weights - 1;
                hr_cmb(i, j, k) = sum(N_hr_cmb(:)) == 27;
                
                % Ordinary CMBs.
                N_cmb = dif < weights;
                cmb(i, j, k) = sum(N_cmb(:)) == 27;
            end
        end
    end
    
    % Layer labeling.
    disp('Labeling layers...');
    layers = ceil(dt / wf);
    
    % Core CMB detection.
    % But one caveat: only HR CMB where used.
    disp('Detecting core CMBs...');
    core_cmb = zeros(dims);
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~hr_cmb(i, j, k)
                    continue;
                end
                    
                % This voxel's sorroundings.
                ri = i - 1: i + 1;
                rj = j - 1: j + 1;
                rk = k - 1: k + 1;
                
                % Voxels in N26 in the same layer.
                N = layers(ri, rj, rk);
                N_same_layer = N == layers(i, j, k);
                neighbors_in_layer = sum(N_same_layer(:));
                
                % First step: CMB's with only CMB's as neighbors
                % in the same layer are marked as core CMB's.           
                
                % CMB's in N26 in the same layer.
                N_cmb = hr_cmb(ri, rj, rk);
                N_cmb_layer = N_cmb & N_same_layer;
                cmb_in_layer = sum(N_cmb_layer(:));
                
                % If every voxel in the neighborhood in the same layer is
                % a CMB, this voxel is a core CMB.
                core_cmb(i, j, k) = neighbors_in_layer == cmb_in_layer;
    
                % Second step: neighbor CMB of core CMB's in
                % the same layer are also marked as core CMB's.
                if(~core_cmb(i, j, k))
                    continue;
                end
                
                % Mark as core every CMB in N26 in the same layer,
                % and don't overwrite previously marked cores.
                core_cmb(ri, rj, rk) = N_cmb_layer | core_cmb(ri, rj, rk);
            end
        end
    end
    
    % Non-core CMB filtering.          
    % First criterion: we only keep
    % CMBs under the convexity threshold.
    disp('Filtering non-core CMB...');
    first_criterion = zeros(dims);
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~cmb(i, j, k)
                    continue;
                end
                N = dt(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1);
                first_criterion(i, j, k) = nnz(N > dt(i, j, k)) < TAU;
            end
        end
    end
    
    % Second criterion: we discard CMBs which belong
    % to a connected component consisting solely of
    % non-core CMB.
    second_criterion = zeros(dims);
    discovered = zeros(dims);
    cmp_st = Stack();
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~cmb(i, j, k) || ...
                   discovered(i, j, k)
                    continue;
                end
                
                % Traverse component.
                % DFS-search (26-connectivity is a guess).
                cmp_st.push([i, j, k]);
                only_non_core = true;
                component = zeros(dims);
                while cmp_st.size() > 0
                    
                    % If current voxel belongs to the core,
                    % we can't remove this component.
                    sp = cmp_st.pop();
                    si = sp(1);
                    sj = sp(2);
                    sk = sp(3);
                    component(si, sj, sk) = 1;
                    if core_cmb(si, sj, sk)
                        only_non_core = false;
                    end
                    
                    % We add undiscovered neighbors to the stack.
                    for di = si - 1: si + 1
                        for dj = sj - 1: sj + 1
                            for dk = sk - 1: sk + 1
                                if ~cmb(di, dj, dk) || ...
                                   discovered(di, dj, dk)
                                    continue;
                                end
                                discovered(di, dj, dk) = 1;
                                cmp_st.push([di, dj, dk]);
                            end
                        end
                    end
                end
                
                % Filter out non-core components.
                if ~only_non_core
                    second_criterion = second_criterion | component;
                end
            end
        end
    end
    
    % For curve skeletonization, anchors are the highly
    % relevant CMB's, which OF COURSE include core CMB's.
    % So, no further filtering is required at all. This would
    % be different for surface skeletonization.
    anchors = core_cmb | (first_criterion & second_criterion & hr_cmb);
    
    % Iterative thinning is done at increasing distance values,
    % preserving anchors and linking voxels. This process yields a 
    % pseudo-surface skeleton (PSS).
    disp('Linking voxels...');
    pss = V;
    linking_voxels = zeros(dims);
    
    % We remove simple points for each distance value.
    for d = 1: max(dt(:))
        I = find(V & dt == d);
        [i, j, k] = ind2sub(dims, I);
        for l = 1: numel(i)
            
            % We keep this voxel if it meets the 3 conditions.
            pss(i(l), j(l), k(l)) = ...
                linking_voxels(i(l), j(l), k(l)) || ...
                anchors(i(l), j(l), k(l)) || ...
                ~is_simple(pss, {i(l), j(l), k(l)});
        end
        
        % After current distance value is done, we must
        % connect selected voxels to the inside.
        I = find(pss & dt == d);
        [i, j, k] = ind2sub(dims, I);
        for l = 1: numel(i)
            
            % This voxel's sorroundings.
            ri = i(l) - 1: i(l) + 1;
            rj = j(l) - 1: j(l) + 1;
            rk = k(l) - 1: k(l) + 1;

            % Maximum directional derivative on N26.
            N = dt(ri, rj, rk);
            N(N <= d) = -Inf;
            
            if nnz(N == -Inf) == 27
                continue;
            end
            
            dd = (N - dt(i(l), j(l), k(l))) ./ weights;
            max_dd = dd == max(dd(:));
            
            % Linking voxels are marked as such.
            linking_voxels(ri, rj, rk) = linking_voxels(ri, rj, rk) | max_dd;
        end
    end
    
    % Thinning of the PSS.
    disp('Thinning PSS...');
    pss = dd_skeleton_thinning(pss, 1);

    % Pruning.
    disp('Pruning PSS...');
    pss = dd_skeleton_pruning(pss, anchors, theta1, theta2);
    
    disp('Voxel classification...');
    line_voxels = zeros(dims);
    ending_voxels = zeros(dims);
    branching_voxels = zeros(dims);
    internal_voxels = zeros(dims);
    junction_voxels = zeros(dims);
    extended_junction_voxels = zeros(dims);
    
    % Get PSS' linear indices.
    pss_li = find(pss);
    [i, j, k] = ind2sub(dims, pss_li);
    
    % A voxel with at most two object neighbors that are
    % not neighbors of each other is a line voxel.
    for l = 1: numel(i)
        
        % This voxel's sorroundings.
        ri = i(l) - 1: i(l) + 1;
        rj = j(l) - 1: j(l) + 1;
        rk = k(l) - 1: k(l) + 1;        
        N26 = pss(ri, rj, rk);
        
        % Count neighbors.
        neighbors_count = sum(N26(:)) - 1;
        
        % Ending voxels are initially classified as line.
        if neighbors_count == 1
            line_voxels(i(l), j(l), k(l)) = 1;
        
        % If there are two neighbors, check their connectivity.
        elseif neighbors_count == 2
            neighbor1 = zeros(1, 3);
            neighbor2 = zeros(1, 3);
            for ni = ri
                for nj = rj
                    for nk = rk
                        if ~pss(ni, nj, nk) || ...
                           (ni == i(l) && nj == j(l) && nk == k(l))
                            continue;
                        end
                        if ~neighbor1
                            neighbor1 = [ni nj nk];
                        else
                            neighbor2 = [ni nj nk];
                        end
                    end
                end
            end

            % Check whether the two neighbors are (26-?)adjacent.
            are_adjacent = false;
            for ni = neighbor1(1) - 1: neighbor1(1) + 1
                for nj = neighbor1(2) - 1: neighbor1(2) + 1
                    for nk = neighbor1(3) - 1: neighbor1(3) + 1
                        if ni == neighbor2(1) && ...
                           nj == neighbor2(2) && ...
                           nk == neighbor2(3)
                            are_adjacent = true;
                        end
                    end
                end
            end
            
            % If they aren't, this is a line voxel.
            line_voxels(i(l), j(l), k(l)) = ~are_adjacent;
        end
    end
    
    % Get line voxels' linear indices.
    line_voxels_li = find(line_voxels);
    [i, j, k] = ind2sub(dims, line_voxels_li);
    
    % A line voxel with only one neighbor that is a line
    % voxel is an ending voxel.
    for l = 1: numel(i)
        
        % This voxel's sorroundings.
        ri = i(l) - 1: i(l) + 1;
        rj = j(l) - 1: j(l) + 1;
        rk = k(l) - 1: k(l) + 1;
        N = pss(ri, rj, rk);
        N_line_voxels = line_voxels(ri, rj, rk);
        
        % There must be just 2 line voxels.
        if sum(N(:)) == 2 && sum(N_line_voxels(:)) == 2
            ending_voxels(i(l), j(l), k(l)) = 1;
            line_voxels(i(l), j(l), k(l)) = 0;
        end
    end
    fprintf('Line voxels: %d\n', sum(line_voxels(:)));
    fprintf('Ending voxels: %d\n', sum(ending_voxels(:)));
    
    % Get unclassified voxels.
    unclassified = find(pss & ~(line_voxels | ending_voxels));
    [i, j, k] = ind2sub(dims, unclassified);    
    
    % A not yet classified voxel with a neighboring line
    % voxel is a branching voxel.
    for l = 1: numel(i)
        
        % This voxel's sorroundings.
        ri = i(l) - 1: i(l) + 1;
        rj = j(l) - 1: j(l) + 1;
        rk = k(l) - 1: k(l) + 1;
        
        % We search for a single line voxel.
        N_lv = line_voxels(ri, rj, rk);
        branching_voxels(i(l), j(l), k(l)) = sum(N_lv(:)) > 0;
    end
    fprintf('Branching voxels: %d\n', sum(branching_voxels(:)));
    
    % Get unclassified voxels again.
    unclassified = find(pss & ...
        ~(line_voxels | ending_voxels | branching_voxels));
    [i, j, k] = ind2sub(dims, unclassified);
    fprintf('Still unclassified: %d\n', nnz(unclassified(:)));
    
    % A not yet classified voxel p with c ~= 1 and with a
    % neighboring branching voxel is a branching voxel.
    more_branching_voxels = zeros(dims);
    for l = 1: numel(i)
        if get_pudney_c(pss, {i(l), j(l), k(l)}) == 1
            continue;
        end
        
        % This voxel's sorroundings.
        ri = i(l) - 1: i(l) + 1;
        rj = j(l) - 1: j(l) + 1;
        rk = k(l) - 1: k(l) + 1;
        
        % We search for a single branching voxel.
        N_bv = branching_voxels(ri, rj, rk);
        more_branching_voxels(i(l), j(l), k(l)) = sum(N_bv(:)) > 0;
    end
    branching_voxels = branching_voxels | more_branching_voxels;
    fprintf('More branching voxels: %d\n', sum(branching_voxels(:)));
    
    % Get unclassified voxels once again.
    unclassified = find(pss & ...
        ~(line_voxels | ending_voxels | branching_voxels));
    [i, j, k] = ind2sub(dims, unclassified); 
    
    % A not yet classified voxel whose object neighbors are
    % all branching voxels is a branching voxel.
    more_branching_voxels = zeros(dims);
    for l = 1: numel(i)
        
        % This voxel's sorroundings.
        ri = i(l) - 1: i(l) + 1;
        rj = j(l) - 1: j(l) + 1;
        rk = k(l) - 1: k(l) + 1;
        
        % We search for many branching voxels.
        N = pss(ri, rj, rk);
        N_bv = branching_voxels(ri, rj, rk);
        more_branching_voxels(i(l), j(l), k(l)) = sum(N_bv(:)) == sum(N(:)) - 1;
    end
    branching_voxels = branching_voxels | more_branching_voxels;
    fprintf('Even more branching voxels: %d\n', sum(more_branching_voxels(:)));
    
    % Get unclassified voxels... again.
    unclassified = find(pss & ...
        ~(line_voxels | ending_voxels | branching_voxels));
    [i, j, k] = ind2sub(dims, unclassified);
    fprintf('Still unclassified: %d\n', nnz(unclassified(:)));
    
    % Any not yet classified voxel with c* ~= 1 is an internal voxel.
    for l = 1: numel(i)
        internal_voxels(i(l), j(l), k(l)) = get_pudney_c_star(pss, {i(l), j(l), k(l)}) ~= 1;
    end
    fprintf('Internal voxels: %d\n', sum(internal_voxels(:)));
    
    % Internal voxels' linear indices.
    internal_voxels_li = find(internal_voxels);
    [i, j, k] = ind2sub(dims, internal_voxels_li);    
    
    % An internal voxel p with more than two 6-connected
    % components of background voxels in N(p) face or
    % edge-adjacent to p (condition 1),
    % or being any of the eight internal
    % voxels in a 2x2x2 configuration is reclassified as a
    % junction voxel (condition 2).
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~internal_voxels(i, j, k)
                    continue;
                end
                
                % Condition 1.
                % Background of neighborhood.
                Nb = ~pss(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1);

                % 6-connected components.
                comps = bwconncomp(Nb, 6);

                % 6-adjacent to current voxel.
                % As in c*'s algorithm, PixelIdxList is a cell array
                % where each vector represents a connected component.
                cc_indexed = comps.PixelIdxList;

                % Linear indexes assigned for 18-neighborhood.
                cool_indexes = [2 4 5 6 8 10 11 12 13 14 15 16 17 18 20 22 23 24 26];

                % We must check every found connected component and only
                % regard those 6-adjacent to P.
                cc_count = 0;
                for c = 1: comps.NumObjects()
                    comp_vector = cc_indexed{c};
                    for n = 1: numel(comp_vector)
                        if any(cool_indexes == comp_vector(n))
                            cc_count = cc_count + 1;
                            break;
                        end
                    end
                end
                if cc_count > 2
                    internal_voxels(i, j, k) = 0;
                    junction_voxels(i, j, k) = 1;
                    continue;
                end
                
                % Condition 2
                % We check every 2x2x2 neighbor containing current voxel.
                if nnz(pss(i - 1: i, j - 1: j, k - 1: k)) == 8 || ...
                   nnz(pss(i - 1: i, j - 1: j, k: k + 1)) == 8 || ...
                   nnz(pss(i - 1: i, j: j + 1, k - 1: k)) == 8 || ...
                   nnz(pss(i - 1: i, j: j + 1, k: k + 1)) == 8 || ...
                   nnz(pss(i: i + 1, j - 1: j, k - 1: k)) == 8 || ...
                   nnz(pss(i: i + 1, j - 1: j, k: k + 1)) == 8 || ...
                   nnz(pss(i: i + 1, j: j + 1, k - 1: k)) == 8 || ...
                   nnz(pss(i: i + 1, j: j + 1, k: k + 1)) == 8
                    internal_voxels(i, j, k) = 0;
                    junction_voxels(i, j, k) = 1;
                    continue;
                end
            end
        end
    end
    
    fprintf('internal voxels: %d\n', nnz(internal_voxels));
    fprintf('junction voxels: %d\n', nnz(junction_voxels));
    
    % An internal voxel having a face or an edge neighbor
    % classified as a junction voxel is reclassified as an
    % extended junction voxel.
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~internal_voxels(i, j, k)
                    continue;
                end
                Nj = junction_voxels(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1);
                N18 = mask_18 .* Nj;
                if nnz(N18) > 0
                    internal_voxels(i, j, k) = 0;
                    extended_junction_voxels(i, j, k) = 1;                    
                end
            end
        end
    end
    fprintf('Junction voxels: %d\n', nnz(junction_voxels));
    fprintf('Extended junction voxels: %d\n', nnz(extended_junction_voxels));
    
    % Any not yet classified voxel is a bordering voxel.
    bordering_voxels = pss & ~(line_voxels | ending_voxels | ...
        branching_voxels | internal_voxels | junction_voxels | ...
        extended_junction_voxels);
    fprintf('Bordering voxels: %d\n', nnz(bordering_voxels));
    
    % Distance transform of PSS.
    % Not an ordinary distance transform, but with respect to
    % the set of bordering voxels.
    disp('Calculating cmb in PSS patches...');
    internal_patches = branching_voxels | internal_voxels | ...
        bordering_voxels;
    dt_ip = weighted_dist_referenced( ...
        internal_patches, ...
        bordering_voxels, ...
        junction_voxels, ...
        extended_junction_voxels, ...
        [3 4 5]);
    
    % Replace 3 => 1.
    dt_ip(dt_ip == wf) = 1;
    
    % CMB in PSS.
    cmb_pss = zeros(dims);
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~pss(i, j, k)
                    continue;
                end
                N = dt_ip(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1);
                if nnz(N - dt_ip(i, j, k) < weights - 1) == 27
                    cmb_pss(i, j, k) = 1;
                end
            end
        end
    end
    
    % Remove junction voxels.
    cmb_pss = cmb_pss & ~junction_voxels;
    
    % We take as new anchor points only CMBs that were
    % accepted as anchor points a long time ago.
    % Also, line, branching and junction voxels.
    anchors_pss = (cmb_pss & anchors) | line_voxels | branching_voxels | junction_voxels;
    
    % We remove simple voxels at an increasing distance.
    % This yields the nearly-thin curve skeleton (NTCS).
    ntcs = pss;
    maxd_pss = max(dt_ip(:));
    for d = 1: maxd_pss
        for i = 2: dims(1) - 1
            for j = 2: dims(2) - 1
                for k = 2: dims(3) - 1
                    if ~ntcs(i, j, k) || ...
                       dt_ip(i, j, k) ~= d || ...
                       anchors_pss(i, j, k)
                        continue;
                    end
                    if is_simple(ntcs, {i, j, k})
                        ntcs(i, j, k) = 0;
                    end
                end
            end
        end
    end
    
    % NTCS thinning.
    % Apply mask-based thinning until no changes are observed.
    disp('Thinning NTCS...');
    ntcs = dd_skeleton_thinning(ntcs, 2);
    
    % NTCS pruning.
    % Anchors are PSS's CMBs plus line and junction voxels.
    disp('Pruning NTCS...');
    anchors_ntcs = (anchors & cmb_pss) | line_voxels | junction_voxels;
    skel = ntcs;
    return;
    ntcs = dd_skeleton_pruning(ntcs, anchors_ntcs, theta1, theta2);
    
    skel = ntcs(pad + 1: end - pad, pad + 1: end - pad, pad + 1: end - pad);
    return;
    
    for step = 1: 2
        fprintf('Pruning NTCS: step %d...\n', step);
        deleted = 0;
        
        % We find NTCS' end points before each step.
        ntcs_end_points = zeros(dims);
        for i = 2: dims(1) - 1
            for j = 2: dims(2) - 1
                for k = 2: dims(3) - 1
                    if ~ntcs(i, j, k)
                        continue;
                    end

                    % End points have exactly one neighbor.
                    N26 = ntcs(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1);
                    ntcs_end_points(i, j, k) = sum(N26(:)) == 2;
                end
            end
        end
        
        % We follow each branch starting from found end points.
        branch_number = 1;
        for i = 2: dims(1) - 1
            for j = 2: dims(2) - 1
                for k = 2: dims(3) - 1
                    if ~ntcs_end_points(i, j, k)
                        continue;
                    end
                    
                    fprintf('NTCS pruning step %d: Following branch %d...\n', step, branch_number);
                    branch_number = branch_number + 1;

                    % N: anchor points along the branch.
                    N = 0;

                    % L: branch length.
                    L = 1;

                    % We follow the branch.
                    branch = zeros(1, 3);
                    branch(L, :) = [i, j, k];
                    while true
                        bi = branch(L, 1);
                        bj = branch(L, 2);
                        bk = branch(L, 3);

                        % We count anchor points found along the branch.
                        N = N + anchors_ntcs(bi, bj, bk);

                        % We stop when we find a branching point.
                        N26 = ntcs(bi - 1: bi + 1, bj - 1: bj + 1, bk - 1: bk + 1);
                        neighbors_count = nnz(N26) - 1;
                        if L > 1 && neighbors_count ~= 2
                            break;
                        end

                        % Look for the unvisited neighbor.
                        found = false;
                        for ni = bi - 1: bi + 1
                            for nj = bj - 1: bj + 1
                                for nk = bk - 1: bk + 1
                                    if ~ntcs(ni, nj, nk) || ...
                                       (ni == bi && nj == bj && nk == bk)
                                        continue;
                                    end

                                    % We must also ensure that we aren't
                                    % going back through the branch.
                                    if L > 1
                                        pbi = branch(L - 1, 1);
                                        pbj = branch(L - 1, 2);
                                        pbk = branch(L - 1, 3);
                                        if pbi == ni && pbj == nj && pbk == nk
                                            continue;
                                        end
                                    end

                                    % The next voxel is added to
                                    % current branch voxels list.
                                    L = L + 1;
                                    branch(L, :) = [ni, nj, nk];
                                    found = true;
                                    break;
                                end
                                if found
                                    break;
                                end
                            end
                            if found
                                break;
                            end
                        end
                    end

                    % 1st pruning step: N <= theta1
                    % 2nd pruning step: N/L <= theta2.
                    if (step == 1 && N <= theta1) || ...
                       (step == 2 && N / L <= theta2)
                        for bpi = 1: L - 1
                            bi = branch(bpi, 1);
                            bj = branch(bpi, 2);
                            bk = branch(bpi, 3);

                            % Prune point.
                            ntcs(bi, bj, bk) = 0;
                            deleted = deleted + 1;
                        end
                        
                        % Check branch point.
                        bi = branch(L, 1);
                        bj = branch(L, 2);
                        bk = branch(L, 3);
                        if is_simple(ntcs, {bi, bj, bk})
                            ntcs(bi, bj, bk) = 0;
                            deleted = deleted + 1;
                        end
                    end
                end
            end
        end
        
        % Remove marked internal points.
        %ntcs = ntcs & ~deletable;
        
        % Remove simple marked branching points.
%         for i = 2: dims(1) - 1
%             for j = 2: dims(2) - 1
%                 for k = 2: dims(3) - 1
%                     if ~deletable_branch_points(i, j, k)
%                         continue;
%                     end
%                     if is_simple(ntcs, {i, j, k})
%                         ntcs(i, j, k) = 0;
%                         deleted = deleted + 1;
%                     end
%                 end
%             end
%         end

        fprintf('pruned %d voxels on stage %d\n', deleted, step);
    end
    
    % Remove all remaining simple points.
%     for i = 2: dims(1) - 1
%         for j = 2: dims(2) - 1
%             for k = 2: dims(3) - 1
%                 if ~ntcs(i, j, k)
%                     continue;
%                 elseif is_simple(ntcs, {i, j, k}) && nnz(ntcs(i-1:i+1,j-1:j+1,k-1:k+1)) > 2
%                     ntcs(i, j, k) = 0;
%                 end
%             end
%         end
%     end
    
    % Remove padding.
    skel = ntcs(pad + 1: end - pad, pad + 1: end - pad, pad + 1: end - pad);
    
end   