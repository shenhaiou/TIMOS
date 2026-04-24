function tetrahedron_to_vtk(node, tetrahedron, vol_fluence, output_filename, title)

%****************************************************************************
%% Tetrahedron_TO_VTK writes out a 3D FEM fluence result to VTK format.
%  Modified:
%    05 October 2011
%  Author:
%    Haiou Shen
%  Parameters:
%    node(3, num_node).
%    tetrahedron(4, num_element).
%    vol_fluence(num_element, 1).
%    OUTPUT_FILENAME.
%    TITLE: a title for the data.
%
[ num_dim, num_node ] = size ( node );
if(num_dim ~= 3)
    warning('The node should be an array of size [3, n]!');
    return;
end
[ element_dim, num_element,   ] = size ( tetrahedron );
if(element_dim ~= 4)
    warning('The tetrahedron should be an array of size [4, m]!');
    return;
end
%
%  Open the output file.
%
if ( isempty ( output_filename ) )
    output_filename = 'timos_result.vtk';
end

fid = fopen ( output_filename, 'w' );

fprintf ( fid, '# vtk DataFile Version 2.0\n' );
fprintf ( fid, '%s\n', title );
fprintf ( fid, 'ASCII\n' );
fprintf ( fid, '\n' );
fprintf ( fid, 'DATASET UNSTRUCTURED_GRID\n' );

fprintf ( fid, 'POINTS %d double\n', num_node );
for idx = 1 : num_node
    fprintf ( fid, '  %f  %f  %f\n', node(:,idx) );
end
fprintf ( fid, '\n' );

fprintf ( fid, 'CELLS  %d  %d\n', num_element, num_element*5 );
for idx = 1 : num_element
    fprintf ( fid, '  4 %d %d %d %d \n', tetrahedron(:,idx) - 1 );
end

fprintf ( fid, '\n' );
fprintf ( fid, 'CELL_TYPES %d\n', num_element );
for idx = 1:num_element
    fprintf ( fid, '10\n');
end
fprintf (fid, '\n');

fprintf ( fid, 'CELL_DATA %d\n', num_element );
fprintf ( fid, 'SCALARS Fluence double\n' );
fprintf ( fid, 'LOOKUP_TABLE default\n' );
fprintf ( fid, '%f \n', vol_fluence);
fprintf ( fid, '\n');
fclose  ( fid);

fprintf ( 1, '\n' );
fprintf ( 1, 'Tetrahedron_to_VTK ... Done\n' );
return
end
