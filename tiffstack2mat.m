% TIFF stack to Matlab matrix parser.
% Stack folder is read in natural alphabetical fashion.
% input: stack directory path.
% output: 3D matrix.
function V = tiffstack2mat(path, crop)

    fnames = dir(path);
    
    % We sort file names using natural sort, instead of
    % ACII sort, so, for instance, tr10.tif precedes tr2.tif.
    tiff_fnames = {};
    tiff_count = 0;
    for i = 1: length(fnames)
        fname = fnames(i).name;
        if ~isempty(strfind(fname, '.tif'))
            tiff_count = tiff_count + 1;
            tiff_fnames{tiff_count} = fname;
        end
    end
    tiff_fnames = sort_nat(tiff_fnames);
    V = [];
    offset = 0;
    for i = 1 + offset: tiff_count
        disp(strcat(path, '\', tiff_fnames{i}));
        
        % We read the image and add it to return variable.
        img = imread(strcat(path, '\', tiff_fnames{i}));
        if ndims(img) > 2
            img = img(:, :, 1);
        end
        img = im2bw(img);
        V(:, :, i - offset) = img;
        if i >= 151 + offset
            break;
        end
    end
    
    % Convert to logical, maybe.
    V = V > 0;
    
    % Crop the stack maybe.
    if nargin > 1 && crop
        
        % Find index ranges.
        [i, j, k] = ind2sub(size(V), find(V));
        irange = min(i): max(i);
        jrange = min(j): max(j);
        krange = min(k): max(k);
        V = V(irange, jrange, krange);
    end
   
end