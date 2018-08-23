function [] = ICPSkeleton(formatSpec, folder, folderout)
% clear all;
% close all;
% formatSpec = 'D:\\GitHub\\data\\%s\\airplane_0%d.ply';
% folder = 'Skeletons';
% ind = max(strfind(folder, '_'));
% folderout = ['skeleton_aligned_', folder(ind+1 : length(folder))]
% folder = 'gt_skeleton';
% folderout = '.';

useMunkres = 0;
useAuction = 1;
mex -largeArrayDims auctionAlgorithmSparseMex.cpp -lut

paths = sprintf(formatSpec, folder, 627);
ptClouds = pcread(paths);
scale = 1000;
Source = double(ptClouds.Location*scale);

%tree
dist = pdist2(Source, Source);
gra = graph(dist);
[minitree,pred] = minspantree(gra);% 最小生成树
searchIdx = dfsearch(minitree, 1);
Source = Source(searchIdx, :);
ptClouds = pointCloud(Source);

npoints = ptClouds.Count;
Color = jet(npoints);

for ii = 650 : 705

patht = sprintf(formatSpec, folder, ii);  
ptCloudt = pcread(patht);
Target = double(ptCloudt.Location*scale);
xyzt = ptCloudt.Location;

% Specify that the source deformations should be plotted.
Options.plot = 1;
Options.stiffFromTarget = 0;
Options.useHungarian = 0;
Options.useAuction = 1;
Options.biDirectional = 1;

pointsTransformed = Source;
[pointsTransformed, X, assignment] = nricp_ske(Source, Target, Options);

if useAuction == 1
      dis = pdist2(pointsTransformed, Target);
      dis = dis.*dis;
      %把求最大值的问题变成最小值
      disMax = -dis + max(max(dis));
      min(min(disMax));
      costmat = double(disMax*100);    
      assignment = sparseAssignmentProblemAuctionAlgorithm(costmat);    
elseif useMunkres == 1
      dis = pdist2(pointsTransformed, Target);
      assignment = munkres(dis);
end

writePath = sprintf(formatSpec,folderout, ii);
alignedCloud = pointCloud(xyzt(assignment, :));
pcwrite(alignedCloud, writePath);

% writePath = sprintf(formatSpec,folderout, 627);
% trCloud = pointCloud(pointsTransformed/scale);
% pcwrite(trCloud, writePath);

figure(ii);
% plotlines(alignedCloud.Location, Color);
pcshow(alignedCloud.Location, Color, 'MarkerSize', 50);
figure();
pcshow(ptClouds.Location, Color, 'MarkerSize', 50);
% plotlines(ptClouds.Location, Color);

end

end
