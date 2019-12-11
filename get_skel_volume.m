function vol = get_skel_volume(skel, edt)

    % Find the skel amongst the darkness.
    [i, j, k] = ind2sub(find(skel));
    
    % Inflate the skel with CMBs
    inflated = skel;
    for l = 1: numel(i)
        r = edt(i(l), j(l), k(l));
        for ii = 1: size(skel, 1)
            for jj = 1: size(skel, 2)
                for kk = 1: size(skel, 3)
                    di = (i(l) - ii)^2;
                    dj = (j(l) - jj)^2;
                    dk = (k(l) - kk)^2;
                    if sqrt(di + dj + dk) <= r
                       inflated(ii, jj, kk) = 1;
                    end
                end
            end
        end
    end
    
    % Volume is easy.
    vol = sum(inflated > 0);
end

