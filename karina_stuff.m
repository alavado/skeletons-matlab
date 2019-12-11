fname = 'Foxd3-SNAP25-ZT0_34_2_proy';
folder = 'C:\Users\aleja\ownCloud\Skeleton\Karina\TOTAL_Mascaras_proy\ZTmañana_dia1\ZT0_dia1\';
out_folder = 'C:\Users\aleja\ownCloud\Skeleton\Karina\TOTAL_Mascaras_proy\';
M = multitiff2mat(sprintf('%s%s.tif', folder, fname));
%skel = thinning3D2(M);
[skel, anchors_ntcs] = dd_skeleton_v2(M, 4, 0.25);
karina_colored_xraw(skel, sprintf('%s%s_arcelli.xraw', out_folder, fname));
%skel_simple = thinning3D2(M);