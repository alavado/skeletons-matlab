% Number of 26-connected components in N*(p).
% input:
%   V: 3D binary volume.
%   p: point of interest.
% output:
%   c: number of 26-connected components in N*(p).
function c = get_pudney_c(V, p)

    i = p{1};
    j = p{2};
    k = p{3};
    n_26 = V(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1);
    n_26(2, 2, 2) = 0;

    comps = bwconncomp(n_26);
    c = comps.NumObjects();
end

