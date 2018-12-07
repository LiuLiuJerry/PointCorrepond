clear all;
close all;
%% sort, oriented
formatSpec = 'D:\\GitHub\\data\\%s\\%s_%.4d.ply';

aligned_name = 'mesh'
leg_name = 'leg'
head_name = 'head'

filename_ort = 'Octopus\octopus2048';
filename_ali = 'Octopus\octopus_HeadAligned';
filename_seg = 'Octopus\octopus_seg'
idx = 0:113;
delta = 100

%匹配源，被匹配的对象
nSource = 2000;
paths = sprintf(formatSpec, filename_ort, aligned_name, nSource)
ptClouds = pcread(paths);

npoints = ptClouds.Count;
Color = jet(npoints);
Source.vertices = double(ptClouds.Location*delta);
Source.normals = ptClouds.Normal;

Options.plot = 1;
Options.stiffFromTarget = 0;
Options.useHungarian = 0;
Options.useAuction = 0;
Options.biDirectional = 0;
r = 3

for ii = idx

    patht = sprintf(formatSpec, filename_ali, aligned_name, ii);  
    ptCloudt = pcread(patht);
    Target.vertices = double(ptCloudt.Location*delta);
    Target.normals = ptCloudt.Normal;
    
    [ vertsTransformed, X, assignment ]  = nricp_head(Source, Target, Options );
    
    [~, dist] = knnsearch(vertsTransformed, Target.vertices, 'K', 1);
    idx_head = find(dist<r);
    head = ptCloudt.Location(idx_head,:);
    legs = ptCloudt.Location(dist>r,:);

    %% write down
    writePath_head = sprintf(formatSpec, filename_seg, head_name, ii);
    writePath_legs = sprintf(formatSpec, filename_seg, leg_name, ii);
    pcwrite(pointCloud(head), writePath_head);
    pcwrite(pointCloud(legs), writePath_legs);
    
%     figure(ii+1);
%     pcshow(head, 'r', 'MarkerSize', 50); hold on;
%     pcshow(legs, 'b', 'MarkerSize', 20);
    
    
end
