function mat2tiffstack(V, path)
    
    if nargin > 1
        for i = 1: size(V, 3)
            imwrite(V(:,:,i), strcat(path, '/', num2str(i), '_test.tiff'));
        end        
    else
        for i = 1: size(V, 3)
            imwrite(V(:,:,i), strcat('datos/tmp/model/', num2str(i), '_test.tiff'));
        end
    end

end

