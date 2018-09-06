clear all;
close all;
%% sort, oriented
formatSpec = 'D:\\GitHub\\data\\SkeletonsbyJerry\\%s\\airplane_0%d.ply';
filename = 'airplane';
filename_ort = 'airplane_ort';
filename_ali = 'airplane_aligned';
idx = 627:726;
% orientation(formatSpec, filename, filename_ort, idx);

% % ICPSkeleton(formatSpec, filename_ort, filename_ali);
useMunkres = 0;
useAuction = 1;
mex -largeArrayDims auctionAlgorithmSparseMex.cpp -lut
%匹配源，被匹配的对象
nSource = 627;
paths = sprintf(formatSpec, filename_ort, nSource);
ptCloud = pcread(paths);
bbox = [ptCloud.XLimits; ptCloud.YLimits; ptCloud.ZLimits];
rs = bbox(:,1)-bbox(:,2);
diameter_S = sqrt(dot(rs,rs));
scale = 100;

% Source key points
[Source, tagsS, dirSource, partsS, key_idxS, boundary_idxS] = key_points(ptCloud, 1);
xyzs = Source;
Source = double(Source*scale);

npoints = ptCloud.Count;
Color = jet(npoints);

for ii = 628

patht = sprintf(formatSpec, filename_ort, ii);  
ptCloudt = pcread(patht);
% Target key points
[Target, tagsT, dirTarget, partsT, key_idxT, boundary_idxT] = key_points(ptCloudt, 1);
xyzt = Target;
Target = double(Target*scale);

% corresParts(Source, Target, key_idxS, key_idxT,  partsS, partsT);

% Specify that the source deformations should be plotted.
Options.plot = 1;
Options.stiffFromTarget = 0;
Options.useHungarian = 1;
Options.useAuction = 0;
Options.biDirectional = 0;
Options.onlyKeyPoints = 0;

pointsTransformed = Source;
[pointsTransformed, X, assignment] = nricp_ske(Source, Target, Options, dirSource, dirTarget, partsS, partsT, key_idxS, key_idxT, boundary_idxS, boundary_idxT);

writePath = sprintf(formatSpec,filename_ali, ii);
Target = Target(assignment, :);
xyzt = xyzt(assignment, :);
alignedCloud = pointCloud(xyzt);

keyptsS_transformed = pointsTransformed(key_idxS, :);
% keyptsT = Target(key_idxT, :);


pcwrite(alignedCloud, writePath);

%% write down
% writePath = sprintf(formatSpec,folderout, 627);
% trCloud = pointCloud(pointsTransformed/scale);
% pcwrite(trCloud, writePath);

% figure(ii);
% % plotlines(alignedCloud.Location, Color);
% pcshow(alignedCloud.Location, Color, 'MarkerSize', 50);
% figure(nSource);
% pcshow(ptCloud.Location, Color, 'MarkerSize', 50);
% % plotlines(ptCloud.Location, Color);
figure(ii);
    for i = unique(tagsS)'
        seg_xyz = Target(tagsS==i, :);
        plot3(seg_xyz(:, 1), seg_xyz(:, 2), seg_xyz(:, 3), '.', 'MarkerSize', 6); hold on;
    end
    for i = partsS'
        plot3(Target(i,1), Target(i,2), Target(i,3)); hold on;
    end
% key points
scatter3(Target(key_idxS, 1), Target(key_idxS, 2), Target(key_idxS, 3), 30, jet(size(key_idxS,1)), 'filled'); hold on;
% boundary points
scatter3(Target(boundary_idxS, 1), Target(boundary_idxS, 2), Target(boundary_idxS, 3), 30, jet(size(boundary_idxS,1)), 'filled'); hold on;
    axis equal;
    hold off;
    
 figure(nSource);
    for i = unique(tagsS)'
        seg_xyz = Source(tagsS==i, :);
        plot3(seg_xyz(:, 1), seg_xyz(:, 2), seg_xyz(:, 3), '.', 'MarkerSize', 6); hold on;
    end
    for i = partsS'
        plot3(Source(i,1), Source(i,2), Source(i,3)); hold on;
    end
% key points
scatter3(Source(key_idxS, 1), Source(key_idxS, 2), Source(key_idxS, 3), 30, jet(size(key_idxS,1)), 'filled'); hold on;
% boundary points
scatter3(Source(boundary_idxS, 1), Source(boundary_idxS, 2), Source(boundary_idxS, 3), 30, jet(size(boundary_idxS,1)), 'filled'); hold on;

    axis equal;
    hold off;
    
end
