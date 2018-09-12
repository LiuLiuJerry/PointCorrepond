function [xyz, seg, dir, parts, keypts_idx, boundary_idx] = key_points(ptCloud, useAuction)
% %读文件
% clear all;
% close all;
% formatSpec = 'D:\\GitHub\\data\\SkeletonsbyJerry\\%s\\airplane_0%d.ply';
% filename = 'airplane';
% filename_ort = 'airplane_ort';
% filename_ali = 'airplane_aligned';
% useAuction = 1;
% mex -largeArrayDims auctionAlgorithmSparseMex.cpp -lut
% [636, 639, 642, 647];
% for ii= [636, 639, 642, 647]
%     
%     path = sprintf(formatSpec, filename_ort, ii);
%     name = sprintf('airplane_0%d,ply', ii);
%     ptCloud = pcread(path);
    xyz = ptCloud.Location;
    bbox = [ptCloud.XLimits; ptCloud.YLimits; ptCloud.ZLimits];
    rs = bbox(:,1)-bbox(:,2);
    diameter = sqrt(dot(rs,rs));
    xyz_mid = symmetry(xyz);

    %区域生长算法
    [seg, xyz, dir, parts, importantpts, keypts_idx, boundary_idx] = findKeyPoints(xyz_mid, 0.8, diameter*0.035);


% close all; clear all;
end