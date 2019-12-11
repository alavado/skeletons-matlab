% 
function c_star = get_pudney_c_star(V, p)
    
    i = p{1};
    j = p{2};
    k = p{3};
    
    if ~V(i, j, k)
        c_star = 0;
        return;
    end
    
    n_26 = V(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1);
    n_26(2, 2, 2) = 0;

    % P's inverted 18-neighborhood.
    mask_18 = cat(3, [0 1 0; 1 1 1; 0 1 0], ...
                     [1 1 1; 1 1 1; 1 1 1], ...
                     [0 1 0; 1 1 1; 0 1 0]);
    n_18 = ~n_26 .* mask_18;
    n_18(2, 2, 2) = 0;

    % Here, we are interested in 6-connected components.
    comps = bwconncomp(n_18, 6);

    % But also, they must be 6-adjacent to P.
    % PixelIdxList is a cell array where each vector represents
    % a connected component found by function bwconncomp. Indexes
    % are linear, check mask_18's definition.
    cc_C_star = comps.PixelIdxList;

    % Linear indexes assigned to P's 6-neighborhood.
    cool_indexes = [5 11 13 15 17 23];

    % We must check every found connected component and only
    % regard those 6-adjacent to P.
    c_star = 0;
    for i = 1: comps.NumObjects()
        comp_vector = cc_C_star{i};
        CC_counts = 0;
        for j = 1: numel(comp_vector)
            if any(cool_indexes == comp_vector(j))
                CC_counts = 1;
                break;
            end
        end
        c_star = c_star + CC_counts;
    end
    
end

