% Metrics for everyone.
function metrics = get_metrics(path, first_node)

    % Read tiff from path.
    model = tiffstack2mat(path, true);
    
    % Get the biggest component.
    model = get_biggest_connected_component(model);
    
    % Calculate skeletons.
    if nargin > 1
        disp('neurons metrics');
        tic;
        skel_hj = thinning3D2(hj_skeleton3D(model, -20));
        t_hj = toc;
        tic;
        skel_th = thinning3D2(model);
        t_th = toc;
        tic;
        skel_dd = thinning3D2(dd_skeleton_v2(model, 7, 0.25));
        t_dd = toc;
    
        % Read au skeleton from file.
        [~, name, ~] = fileparts(path);
        skel_au = read_skel_obj(strcat('datos/skeletons_au/', name, '_embedding_refinement.obj'), skel_th);
    else
        disp('er metrics');
        
        tic;
        % For 2D models the 2D version is more appropriate.
        skel_hj = thinning3D2(hj_skeleton(model, -40));
        mat2xraw(skel_hj, 'toledo_hj40.xraw');
        t_hj = toc;
        tic;
        skel_th = thinning3D2(model);
        mat2xraw(skel_th, 'toledo_th.xraw');
        t_th = toc;
        tic;
        skel_dd = thinning3D2(dd_skeleton_v2(model, 7, 0.25));
        mat2xraw(skel_dd, 'toledo_dd7.xraw');
        t_dd = toc;
    
        % Read au skeleton from file.
        [~, name, ~] = fileparts(path);
        skel_au = read_skel_obj(strcat('datos/skeletons_au/', name, '_embedding_refinement.obj'), skel_th);
        skel_au = rot90mat(rot90mat(rot90mat(rot90mat(rot90mat(skel_au,3),3),3),2),3);
    end
    
    skels = {skel_th, skel_hj, skel_dd, skel_au};
    skels_codes = {'palagyi', 'siddiqi', 'arcelli', 'au'};
    times = {t_th, t_hj, t_dd, 1};
    metrics_by_algorithm = {'','','',''};
    for s = 1: numel(skels)
        
        mat2xraw(skels{s}, strcat('skels/', name, '_', skels_codes{s}, '.xraw'));
        
        % Other metrics depend on whether first node is specified.
        % Specified => neuron metrics.
        % Format: time, volume, nodes, edges, total_length, paths, angles, degrees, _ary lengths.
        if nargin > 1
            
            % Model name.
            model_name = sprintf('%s', name);

            % Algorithm.
            algorithm = sprintf('%s', skels_codes{s});

            % Execution time.
            model_time = sprintf('%s', num2str(times{s}));

            % Volume (voxels).
            metric_volume = sprintf('%d', nnz(skels{s}));
        
            % Get the matrix and other stuff.
            [A, paths, angles, ~] = get_adjacency_matrix3D(skels{s}, first_node);
            
            % Node count.
            metric_nodes = sprintf('%d', size(A, 1));
            
            % Edges count.
            metric_edges = sprintf('%d', sum(A(:)) / 2);
            
            % Paths and total length.
            total_length = 0;
            paths_str = '';
            for i = 1: numel(paths)
                path = paths{i};
                total_length = total_length + path{4};
                paths_str = strcat(paths_str, ',', num2str(path{4},'%.2f'));
            end
            metric_total_length = sprintf('%s', num2str(total_length,'%.2f'));
            metric_paths = sprintf('[%s]', paths_str(2: end));
            
            % Bifurcation angles.
            angles_str = '';
            for i = 1: numel(angles)
                angles_str = strcat(angles_str, ',', num2str(angles(i),'%.2f'));
            end
            metric_angles = sprintf('[%s]', angles_str(2: end));
            
            % Bifurcation degrees.
            degrees_str = '';
            for i = 1: size(A, 1)
                degrees_str = strcat(degrees_str, ',', num2str(sum(A(i, :))));
            end
            metric_degrees = sprintf('[%s]', degrees_str(2: end));
            
            % Primary, secondary and tertiary lengths.
            
            metrics_by_algorithm{s} = sprintf('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s', ...
                model_name, algorithm, model_time, metric_volume, ...
                metric_nodes, metric_edges, metric_total_length, ...
                metric_paths, metric_angles, metric_degrees);
            
        % Not specified => reticulum metrics.
        % Format: model, time, volume, nodes, edges, total_length, paths, degrees, cycles.
        else
            
            % Model name.
            model_name = sprintf('%s', name);

            % Algorithm.
            algorithm = sprintf('%s', skels_codes{s});

            % Execution time.
            model_time = sprintf('%s', num2str(times{s}));

            % Volume (voxels).
            metric_volume = sprintf('%d', nnz(skels{s}));
        
            % Get the matrix and other stuff.
            [A, paths, ~] = get_adjacency_matrix3D(skels{s});
            
            % Node count.
            metric_nodes = sprintf('%d', size(A, 1));
            
            % Edges count.
            metric_edges = sprintf('%d', sum(A(:)) / 2);
            
            % Paths and total length.
            total_length = 0;
            paths_str = '';
            for i = 1: numel(paths)
                path = paths{i};
                total_length = total_length + path{4};
                paths_str = strcat(paths_str, ',', num2str(path{4},'%.2f'));
            end
            metric_total_length = sprintf('%s', num2str(total_length,'%.2f'));
            metric_paths = sprintf('[%s]', paths_str(2: end));
            
            % Bifurcation degrees.
            degrees_str = '';
            for i = 1: size(A, 1)
                degrees_str = strcat(degrees_str, ',', num2str(sum(A(i, :))));
            end
            metric_degrees = sprintf('[%s]', degrees_str(2: end));
            
            % Chordless cycle count, assuming image is 2D
            comps = bwconncomp(~padarray(skels{s}, [1 1 1]), 4);
            metric_cycles = sprintf('%s', num2str(comps.NumObjects - 1));
            
            % Format: model, time, volume, nodes, edges, total_length, paths, degrees, cycles.
            metrics_by_algorithm{s} = sprintf('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s', ...
                model_name, algorithm, model_time, metric_volume, ...
                metric_nodes, metric_edges, metric_total_length, ...
                metric_paths, metric_degrees, metric_cycles);
        end
    end
    
    % Metrics string.
    metrics = sprintf('%s\n%s\n%s\n%s\n', ...
        metrics_by_algorithm{1}, metrics_by_algorithm{2}, ...
        metrics_by_algorithm{3}, metrics_by_algorithm{4});
end
