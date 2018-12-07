clear all;
close all;
% 和某一个头的刚性对齐
%% sort, oriented
formatSpec = 'D:\\GitHub\\data\\%s\\mesh_%.4d%s';
extension = '.ply';

filename_ort = 'Octopus\octopus2048';
filename_ali = 'Octopus\octopus_HeadAligned';
idx = 0:113;

%匹配源，被匹配的对象
nSource = 2000;
paths = sprintf(formatSpec, filename_ort, nSource, extension)
ptClouds = pcread(paths);
bbox = [ptClouds.XLimits; ptClouds.YLimits; ptClouds.ZLimits];
rs = bbox(:,1)-bbox(:,2);
diameter_S = sqrt(dot(rs,rs));
lamda = 0.8:0.05:1.5

npoints = ptClouds.Count;
Color = jet(npoints);
max_y = max(ptClouds.Location(:,2))

for ii = idx

    patht = sprintf(formatSpec, filename_ort, ii, extension);  
    ptCloudt = pcread(patht);
    xyzt = ptCloudt.Location
    delta_y = max(xyzt(:,2))-max_y
    xyzt(:,2) = xyzt(:,2)-delta_y
    min_mse = -1;
    head_aligned = ptClouds;
    s = 0;
    t_matrix = eye(4)
    for S = lamda
        ptCloudt = pointCloud( (xyzt)*S);
        [tform, new_head, rmse] = pcregistericp(ptClouds, ptCloudt, 'Metric', 'pointToPlane');
        if min_mse == -1 || min_mse > rmse
            min_mse = rmse;
            head_aligned = new_head;
            Cloudt = ptCloudt;
            s = S;
            t_matrix = tform.T;
        end
%         new_head = pctransform(ptClouds, tform);
    end
%     figure(ii);
%     pcshow(head_aligned.Location, 'r', 'MarkerSize', 20); hold on;
%     pcshow(Cloudt.Location, 'b', 'MarkerSize', 20);
    %% write down
    writePath = sprintf(formatSpec, filename_ali, ii, extension);
    pc_aligned = pctransform(Cloudt, affine3d(t_matrix^-1));
    pcwrite(pointCloud(single(pc_aligned.Location)), writePath);
%     figure(ii);
%     pcshow(ptClouds.Location, 'r', 'MarkerSize', 20); hold on;
%     pcshow(pc_aligned.Location, 'b', 'MarkerSize', 20);
    
    
end
