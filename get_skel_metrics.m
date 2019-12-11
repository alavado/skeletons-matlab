function [metric_nodes, metric_edges, metric_total_length, metric_paths, metric_angles, metric_degrees, metric_ip, metric_ep] = get_skel_metrics(fname, first_node)
%GET_SKEL_METRICS Skeletonizes and gets the metrics for a multitiff stack

    % Read tiff from path.
    model = multitiff2mat(fname);
    
    % Keep only largest component.
    model = get_biggest_connected_component(model);

    % Rescale. 1st param. indicates dimension (3->Z); 2nd param, factor.
    model = rescale_neuron(model, 3, 5);

    % Skeletonize using the worst algorithm.
    skel = thinning3D2(model);

    % Get the matrix and other stuff.
    [A, paths, angles, ~] = get_adjacency_matrix3D(skel, first_node);

    % Node count.
    metric_nodes = sprintf('%d', size(A, 1));

    % Edges count.
    metric_edges = sprintf('%d', sum(A(:)) / 2);

    % Paths and total length.
    metric_total_length = 0;
    paths_str = '';
    for i = 1: numel(paths)
        path = paths{i};
        metric_total_length = metric_total_length + path{4};
        paths_str = strcat(paths_str, ',', num2str(path{4},'%.2f'));
    end
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
    
    % Internal and end nodes.
    metric_ip = 0;
    metric_ep = 0;
    disp(A)
    for i = 1: size(A, 1)
        if sum(A(i, :)) == 1
            metric_ip = metric_ip + 1;
        else
            metric_ep = metric_ep + 1;
        end
    end
end

