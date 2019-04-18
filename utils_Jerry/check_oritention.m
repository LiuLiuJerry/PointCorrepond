close all;
clear all;

formatSpec = 'D:/GitHub/data/%s/airplane_%.4d%s';
formatSpec2 = 'D:/GitHub/data/%s/mesh_%.4d%s';
input_file = 'airplane/train';
output_file = 'airplane/airplane';

extension1 = '.2048.ply';
extension2 = '.ply';
delta = 20;
idx = 1:626;

N = 2048;
t = 0;
for i = idx
    path = sprintf(formatSpec, input_file, i, extension1);
    ptcloud = pcread(path);
    
    pts = ptcloud.Location;
    pts2 = pts.*[-1,1,1];
    [ind, d] = knnsearch(pts, pts2, 'K', 1);
    dis = sum(d);
    
    n = N*0.1;
    pts=sortrows(pts, 2);
    tail1 = pts(1:n, :);
    tail2 = pts(end-n:end, :);
    v1 = sum( var(tail1(:, 1)) );
    v2 = sum( var(tail2(:, 1)) );

    if dis < delta && v1<v2
        path_out = sprintf(formatSpec, output_file, t, extension2);
        figure();
        %pcshow(ptcloud);
        pcwrite(ptcloud, path_out);
        t = t+1;
    end
end