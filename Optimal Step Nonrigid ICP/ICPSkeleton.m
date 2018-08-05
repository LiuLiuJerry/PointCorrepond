clear all;
close all;
formatSpec = 'D:\\GitHub\\pointCloud\\data\\%s\\airplane_0%d.ply';
path1 = sprintf(formatSpec, 'ort_skeleton_2048', 628);
path2 = sprintf(formatSpec, 'ort_skeleton_2048', 627);
%path2 = sprintf(formatSpec, 'ort_skeleton_2048', 683);
ptCloud1 = pcread(path1);
ptCloud2 = pcread(path2);

xyz1 = ptCloud1.Location;
%xyz1=sortrows(xyz1, 3); xyz1=sortrows(xyz1, 1); xyz1=sortrows(xyz1, 2);
xyz2 = ptCloud2.Location;
%xyz2=sortrows(xyz2, 3); xyz2=sortrows(xyz2, 1); xyz2=sortrows(xyz2, 2);
Source.vertices = double(xyz1*100);
Target.vertices = double(xyz2*100);

% Specify that the source deformations should be plotted.
Options.plot = 1;
Options.stiffFromTarget = 0;
Options.elasticity = 1;

[pointsTransformed, X] = nricp_ske(Source, Target, Options);