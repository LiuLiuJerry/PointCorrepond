clear all;
close all;
formatSpec = 'D:\\GitHub\\pointCloud\\data\\%s\\airplane_0%d.ply';
folder = 'ort_skeleton_256';
ind = max(strfind(folder, '_'));
folderout = ['skeleton_aligned_', folder(ind+1 : length(folder))]
% folder = 'gt_skeleton';
% folderout = '.';
useMunkres = 0;
useAuction = 1;

paths = sprintf(formatSpec, folder, 627);
ptClouds = pcread(paths);
scale = 1000;
Source = double(ptClouds.Location*scale);

npoints = ptClouds.Count;
Color = jet(npoints);

for ii = 628 : 726

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
      dis = dis.*dis;
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
pcshow(alignedCloud.Location, Color, 'MarkerSize', 50);
figure();
pcshow(ptClouds.Location, Color, 'MarkerSize', 50);

end