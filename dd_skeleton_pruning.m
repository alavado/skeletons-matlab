function skel = dd_skeleton_pruning(pss, anchors, theta1, theta2)
    % The only difference between the two pruning steps is
    % the final criteria used for pruning each branch.

    dims = size(pss);
    for step = 1: 2
        deletable = zeros(dims);
        discovered = zeros(dims);
        deleted = 0;
        
        % Indices for PSS.
        I = find(pss);
        [i, j, k] = ind2sub(dims, I);

        % We find PSS' end points before each step.
        pss_end_points = zeros(dims);
        for l = 1: numel(i)
            N26 = pss(i(l) - 1: i(l) + 1, j(l) - 1: j(l) + 1, k(l) - 1: k(l) + 1);
            
            % An end point has exactly one neighbor.
            pss_end_points(i(l), j(l), k(l)) = sum(N26(:)) == 2;
        end
        
        % Indices for PSS' end points.
        I = find(pss_end_points);
        [i, j, k] = ind2sub(dims, I);
        
        fprintf('tracking branches...\n');
        
        % We track each end point's branch.
        for l = 1: numel(i)

            % Anchor points and branch length.
            anchor_count = 0;
            branch_length = 1;
        
            fprintf('following branch %d...\n', l);

            % We follow the branch.
            branch = zeros(1, 3);
            branch(branch_length, :) = [i(l), j(l), k(l)];
            while true
                bi = branch(branch_length, 1);
                bj = branch(branch_length, 2);
                bk = branch(branch_length, 3);

                % We count anchor points found along the branch.
                anchor_count = anchor_count + anchors(bi, bj, bk);

                % We stop when we find a branching point.
                N26 = pss(bi - 1: bi + 1, bj - 1: bj + 1, bk - 1: bk + 1);
                neighbors_count = sum(N26(:)) - 1;
                if branch_length > 1 && neighbors_count ~= 2
                    break;
                end
                found = false;

                % Look for the unvisited neighbor.
                for ni = bi - 1: bi + 1
                    for nj = bj - 1: bj + 1
                        for nk = bk - 1: bk + 1
                            if ~pss(ni, nj, nk) || ...
                               (ni == bi && nj == bj && nk == bk)
                                continue;
                            end

                            % We must also ensure that we aren't
                            % going back through the branch.
                            if branch_length > 1
                                pbi = branch(branch_length - 1, 1);
                                pbj = branch(branch_length - 1, 2);
                                pbk = branch(branch_length - 1, 3);
                                if pbi == ni && pbj == nj && pbk == nk
                                    continue;
                                end
                            end

                            % Point is added to branch voxels list.
                            branch_length = branch_length + 1;
                            branch(branch_length, :) = [ni, nj, nk];
                            
                            % Stop this.
                            found = true;
                            break;
                        end
                        if found
                            break
                        end
                    end
                    if found
                        break
                    end
                end
            end
            
            % 1st pruning step: N <= theta1
            % 2nd pruning step: N/L <= theta2.
            if (step == 1 && anchor_count <= theta1) || ...
               (step == 2 && anchor_count / branch_length <= theta2)

                % We prune the whole branch save the last voxel.
                for bpi = 1: branch_length - 1
                    bi = branch(bpi, 1);
                    bj = branch(bpi, 2);
                    bk = branch(bpi, 3);             
                    pss(bi, bj, bk) = 0;
                    deleted = deleted + 1;
                end

                % If the last voxel is still a branch point,
                % we prune all simple connected branch points.
                bi = branch(branch_length, 1);
                bj = branch(branch_length, 2);
                bk = branch(branch_length, 3);
                Nbp = pss(bi - 1: bi + 1, bj - 1: bj + 1, bk - 1: bk + 1);
                if sum(Nbp(:)) > 3
                    
                    % Delete branch point if it's simple.
                    if is_simple(pss, {bi, bj, bk})
                        pss(bi, bj, bk) = 0;
                        deleted = deleted + 1;
                    end
                    
                    % Traverse component.
                    % DFS (26-connectivity is again a guess).
                    cmp_st = Stack();
                    cmp_st.push([bi, bj, bk]);
                    while cmp_st.size() > 0
                        sp = cmp_st.pop();
                        si = sp(1);
                        sj = sp(2);
                        sk = sp(3);                                
                        discovered(si, sj, sk) = 1;

                        % We mark simple branching points for removal.
                        if is_simple(pss, {si, sj, sk})
                            deletable(si, sj, sk) = 1;
                        end

                        % Follow branch voxels component.
                        for di = si - 1: si + 1
                            for dj = sj - 1: sj + 1
                                for dk = sk - 1: sk + 1
                                    if ~pss(di, dj, dk) || ...
                                       discovered(di, dj, dk)
                                        continue;
                                    end

                                    % Adjacent branch points
                                    % go to the stack.
                                    Nbp = pss(di - 1: di + 1, ...
                                        dj - 1: dj + 1, ...
                                        dk - 1: dk + 1);
                                    if nnz(Nbp) > 3
                                        cmp_st.push([di, dj, dk]);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        fprintf('removing branch points...\n');
        
        % Simple branch points removal.
        I = find(deletable);
        [i, j, k] = ind2sub(dims, I);
        for l = 1: numel(i)
            if is_simple(pss, {i(l), j(l), k(l)})
                pss(i(l), j(l), k(l)) = 0;
                deleted = deleted + 1;
            end
        end
        fprintf('deleted %d in step %d\n', deleted, step);
                    
%                     % BFS may be better. NO, DFS is the way to go.
%                     q = Queue();
%                     q.enqueue([bi, bj, bk]);
%                     discovered(bi, bj, bk) = 1;
%                     q2 = Queue();
%                     while q.size() > 0
%                         qp = q.dequeue();
%                         q2.enqueue(qp);
%                         for di = qp(1) - 1: qp(1) + 1
%                             for dj = qp(2) - 1: qp(2) + 1
%                                 for dk = qp(3) - 1: qp(3) + 1
%                                     if ~pss(di, dj, dk) || ...
%                                        discovered(di, dj, dk)
%                                         continue;
%                                     end
%                                     Nbp = pss(di - 1: di + 1, ...
%                                         dj - 1: dj + 1, ...
%                                         dk - 1: dk + 1);
%                                     if sum(Nbp(:)) > 3
%                                         discovered(di, dj, dk) = 1;
%                                         q.enqueue([di, dj, dk]);
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                     
%                     while q2.size() > 0
%                         qp = q2.dequeue();
%                         Nbp = pss(qp(1) - 1: qp(1) + 1, ...
%                                 qp(2) - 1: qp(2) + 1, ...
%                                 qp(3) - 1: qp(3) + 1);
%                         if is_simple(pss, {qp(1), qp(2), qp(3)}) && ...
%                            sum(Nbp(:)) > 3
%                             pss(qp(1), qp(2), qp(3)) = 0;
%                             deleted = deleted + 1;
%                         end
%                     end
%                     
%                 end
%             end
%         end
%         fprintf('deleted %d in step %d\n', deleted, step);
%         skel = pss;
%         return;
%     end
    end
    
    skel = pss;
end

