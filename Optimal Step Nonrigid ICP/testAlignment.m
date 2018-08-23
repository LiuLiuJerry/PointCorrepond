%%单纯进行一对一匹配，不进行变形
clear all;
useMunkre = 0;

formatSpec = 'D:\\GitHub\\data\\%s\\airplane_0%d.2048.ply';
folder = 'SkeletonsbyJerry\\airplane_ort';
% folder = 'gt_skeleton';
path1 = sprintf(formatSpec, folder, 627);
path2 = sprintf(formatSpec, folder, 638);
%path2 = sprintf(formatSpec, 'ort_skeleton_2048', 683);
ptCloud1 = pcread(path1);
ptCloud2 = pcread(path2);

vertsSource = ptCloud1.Location;
vertsTarget = ptCloud2.Location;
vertsSource=sortrows(vertsSource, 3); vertsSource=sortrows(vertsSource, 1); vertsSource=sortrows(vertsSource, 2);

dis = pdist2(vertsSource*100000, vertsTarget*100000);

% dis =[
% 
%      1     2     3     4     5
%      5     4     3     2     4
%      4     5     6     3     2
%      2     1     3     0     2
%      3     4     5     1     2];

if useMunkre == 0
    mex -largeArrayDims auctionAlgorithmSparseMex.cpp -lut
    disMax = -dis + max(max(dis));
    min(min(disMax))
    costmat = double(disMax);
    res = sparseAssignmentProblemAuctionAlgorithm(costmat);
else
    res = munkres(dis);
end
targetId = res;
U = vertsTarget(targetId,:);
size(U)
ptCloudsorted = pointCloud(U);

len = size(dis, 1);
Color = jet(len);
figure(1)
    xlabel('x')
    ylabel('y')
    zlabel('z')
pcshow(ptCloudsorted.Location, Color, 'MarkerSize', 50);
figure(2)
pcshow(vertsSource, Color, 'MarkerSize', 50);

