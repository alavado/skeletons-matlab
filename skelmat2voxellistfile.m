function skelmat2voxellistfile(skel, fname)
    f = fopen(fname, 'w');
    [i, j, k] = ind2sub(size(skel), find(skel));
    for l = 1: numel(i)
        fprintf(f, '%d,%d,%d\n', i(l), j(l), k(l));
    end
    fclose(f);
end
