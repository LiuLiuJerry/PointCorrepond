% formatSpec = 'D:\\GitHub\\data\\Chen\\work\\%s\\mesh_%.4d%s';
% input_file = 'decoded_walk_0.2';
% output_file = 'walk_obj';

formatSpec = 'D:/GitHub/data/%s/mesh_%.4d%s';
input_file = 'human/I_squat/meshes';
output_file = 'human/I_squat/squat';

idx = 0:5:249;
extension1 = '.obj';
extension2 = '.obj';
i = 0;
for ii = idx
    path_in = sprintf(formatSpec, input_file, ii, extension1);
    [V,F]=readOBJ(path_in);
    path_out = sprintf(formatSpec, output_file, i, extension1);
    writeOBJ(path_out,V,F);
    i = i+1;
end
