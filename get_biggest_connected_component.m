function I = get_biggest_connected_component(V)
    
    % Label conneted components.
    labels = bwlabeln(V);
    biggest_size = 0;
    for i = 1: max(labels(:))
        label_size = nnz(labels == i);
        if label_size > biggest_size
            biggest_size = label_size;
            biggest_label = i;
        end
    end
    I = labels == biggest_label;
end

