% Necessary parameter: flux threshold.
% Points with an average outward flux above this are left out.
function skel = hj_skeleton3D(V, THRESH)
    
    % Zero-padding is necessary to avoid boundary value conditions.
    pad = 1;
    V = padarray(V, [pad pad pad]);
    dims = size(V);
 
    % Euclidean distance function.
    edt = bwdist(~V, 'euclidean');
    
    % Gradient vector field.
    gradx = zeros(size(V));
    grady = zeros(size(V));
    gradz = zeros(size(V));
    
    % Gradient calculation via 3D Sobel operator.
    kernel = [1.0 4.0 1.0 4.0 16.0 4.0 1.0 4.0 1.0];    
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~V(i, j, k)
                    continue;
                end
                idx = 1;
                gx = 0.0;
                gy = 0.0;
                gz = 0.0;
                for kk = -1: 1: 1
                    for mm = -1: 1: 1
                        kval = kernel(idx);
                        idx = idx + 1;
                        gx = gx + kval * ...
                            (edt(i + kk, j + mm, k + 1) - ...
                            edt(i + kk, j + mm, k - 1));
                        gy = gy + kval * ...
                            (edt(i + kk, j + 1, k + mm) - ...
                            edt(i + kk, j - 1, k + mm));
                        gz = gz + kval * ...
                            (edt(i + 1, j + kk, k + mm) - ...
                            edt(i - 1, j + kk, k + mm));
                    end
                end
                kval = sqrt(gx * gx + gy * gy + gz * gz);
                if(kval >= 1e-7)
                    gradx(i, j, k) = gx / kval;
                    grady(i, j, k) = gy / kval;
                    gradz(i, j, k) = gz / kval;
                end
            end
        end
    end
    
    % Outward normals at 26-neighborhood.
    bla = sqrt(1/3);
    blo = sqrt(1/2);

    % First layer.
    nx(:, :, 1) = [
        -bla  0.0000  bla
        -blo  0.0000  blo
        -bla  0.0000  bla
    ];
    ny(:, :, 1) = -1 * [
         bla  blo  bla
         0.0000  0.0000  0.0000
        -bla -blo -bla
    ];
    nz(:, :, 1) = -1 * [
         bla  blo  bla
         blo  1.0000  blo
         bla  blo  bla
    ];

    % Second layer.
    nx(:, :, 2) = [
        -blo  0.0000  blo
        -1.0000  0.0000  1.0000
        -blo  0.0000  blo
    ];
    ny(:, :, 2) = -1 * [
         blo  1.0000  blo
         0.0000  0.0000  0.0000
        -blo -1.0000 -blo
    ];
    nz(:, :, 2) = [
         0.0000  0.0000  0.0000
         0.0000  0.0000  0.0000
         0.0000  0.0000  0.0000
    ];

    % Third layer.
    nx(:, :, 3) = [
        -bla  0.0000  bla
        -blo  0.0000  blo
        -bla  0.0000  bla
    ];
    ny(:, :, 3) = -1 * [
         bla  blo  bla
         0.0000  0.0000  0.0000
        -bla -blo -bla
    ];
    nz(:, :, 3) = -1 * [
        -bla -blo -bla
        -blo -1.0000 -blo
        -bla -blo -bla
    ];

    % We calculate flux average at each point of the object.
    fluxes = zeros(size(V));        
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~V(i, j, k)
                    continue;
                end
                
                % Flux is averaged amongst neighbors.
                flux = 0;
                for ni = 1: 3
                    for nj = 1: 3
                        for nk = 1: 3

                            % Neighbor coordinates.
                            gi = i + ni - 2;
                            gj = j + nj - 2;
                            gk = k + nk - 2;

                            % Neighbor flux.
                            dot_product = ...
                                nz(ni, nj, nk) * gradx(gi, gj, gk) + ...
                                nx(ni, nj, nk) * grady(gi, gj, gk) + ...
                                ny(ni, nj, nk) * gradz(gi, gj, gk);
                            flux = flux + dot_product;
                        end
                    end
                end                  
                fluxes(i, j, k) = flux;
            end
        end
    end
    
    disp(min(fluxes(:)));
    
    % Our heap is a min queue.
    pq = PriorityQueue();
    
    % To tell which elements are on the queue.
    in_queue = zeros(size(V));
    
    % We search for simple points in the boundary
    % as our starting point for thinning.
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~V(i, j, k)
                    continue;
                end
                
                % Boundary points have at least one background neighbor.
                N26 = V(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1);

                % If this boundary point is also a simple point,
                % we add it to the heap.
                if nnz(N26) < 27 && is_simple(V, {i, j, k})

                    % Negative flux to use min queue as max queue.
                    pq.push(-fluxes(i, j, k), [i j k]);
                    in_queue(i, j, k) = 1;
                end
            end
        end
    end
        
    % Final skeleton placeholder.
    skel = zeros(size(V));
    
    % Repeat until heap is empty.
    while pq.size() > 0
        
        % Removal of first element.
        [P, flux] = pq.pop();
        
        % Point coordinates.
        i = P(1);
        j = P(2);
        k = P(3);
        in_queue(i, j, k) = 0;
                 
        % We only remove simple points.
        if is_simple(V, {i, j, k})
            
            % If it is not an end point or flux is too high, it just
            % cannot be a skeletal point.
            if -flux > THRESH || ~(skel(i, j, k) || nnz(V(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1)) == 2)
                
                % We remove P.
                V(i, j, k) = 0;
                
                % We insert every simple neighbor of P into the heap.
                for ni = i - 1: i + 1
                    for nj = j - 1: j + 1
                        for nk = k - 1: k + 1
                            
                            if V(ni, nj, nk) && ...
                                ~in_queue(ni, nj, nk) && ...
                                is_simple(V, {ni, nj, nk})
                                in_queue(ni, nj, nk) = 1;
                                pq.push(-fluxes(ni, nj, nk), [ni nj nk]);
                            end
                        end
                    end
                end
            
            % If it's an end point and flux is under threshold,
            % we declare current point as a skeletal point.
            else
                skel(i, j, k) = 1;
            end
        end
        
    end
    
    % Final skeleton is what survives, unpadded.
    skel = V(1 + pad: end - pad, 1 + pad: end - pad, 1 + pad: end - pad);
end