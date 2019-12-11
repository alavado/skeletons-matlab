%--------------------------------------------------------------------------
% Author : Rodrigo I. Brevis
% Update : 04-12-2017
% Example: x = 1:nx;
%          y = 1:ny;
%          z = 1:nz;
%          [X,Y,Z] = meshgrid(x,y,z);
%          U = f1(X,Y,Z);
%          V = f2(X,Y,Z);
%          W = f3(X,Y,Z);
%          filename = 'name.vtk';
%          vtk_unstructured_vector_points_3d(X,Y,Z,nx,ny,nz,U,V,W,filename)
%--------------------------------------------------------------------------
function vtk_unstructured_vector_points_3d(X,Y,Z,nx,ny,nz,U,V,W,filename)
num_points = nx*ny*nz;
num_cell = (nx-1)*(ny-1)*(nz-1);
points = [X(:), Y(:), Z(:)]';
%----- Headline
fileID = fopen(filename,'w');
fprintf(fileID, '# vtk DataFile Version 3.0\n');
fprintf(fileID, 'vtk output\n');
fprintf(fileID, 'ASCII\n');
fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(fileID, 'POINTS %12d float\n',num_points);
%--------------------------------------------------------------------------
% DATASET FORMAT
%--------------------------------------------------------------------------
%----- nodes
fprintf(fileID,'%12.8f %12.8f %12.8f\n',points);
%----- cells
fprintf(fileID,'CELLS %12d %12d\n',num_cell,9*num_cell);
vec_cell = 1:num_cell;
vec_points = 1:num_points;
matrix_cell = reshape(vec_cell,[ny-1,nx-1,nz-1]);
matrix_points = reshape(vec_points,[ny,nx,nz]);
matrix_points = matrix_points(1:ny-1,1:nx-1,1:nz-1);
E = 8 + 0*matrix_cell(:);
N1 = matrix_points(:) - 1;
N2 = N1 + 1;
N3 = N1 + ny;
N4 = N3 + 1;
N5 = N1 + ny*nx;
N6 = N5 +1;
N7 = N5 + ny;
N8 = N7 +1;
aux = [E N1 N2 N3 N4 N5 N6 N7 N8]';
fprintf(fileID,'%12d %12d %12d %12d %12d %12d %12d %12d %12d\n',aux);
%----- point type
fprintf(fileID,'CELL_TYPES %12d\n',num_cell);
vec_point_type = 11*ones(num_cell,1);
fprintf(fileID,'%12d\n',vec_point_type);
%--------------------------------------------------------------------------
% DATASET ATTRIBUTE FORMAT
%--------------------------------------------------------------------------
fprintf(fileID,'POINT_DATA %12d\n',num_points);
fprintf(fileID,'VECTORS vectors float\n');
% fprintf(fileID,'LOOKUP_TABLE default\n');
vec_point_data = [U(:) V(:) W(:)]';
fprintf(fileID,'%12.8f %12.8f %12.8f\n',vec_point_data);
%----- colse vkt file
fclose(fileID);