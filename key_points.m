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
%     figure();
%     pcshow(xyz, jet(size(xyz, 1)), 'MarkerSize', 50);
    
    %对称处理
    xyz_sym = xyz.*[-1, 1, 1];
    
    if useAuction == 1
      dis = pdist2(xyz_sym, xyz);
      dis = dis.*dis;
      %把求最大值的问题变成最小值
      disMax = -dis + max(max(dis));
      min(min(disMax));
      costmat = double(disMax*100);    
      assignment = sparseAssignmentProblemAuctionAlgorithm(costmat);    
    else 
      dis = pdist2(xyz_sym, xyz);
      assignment = munkres(dis);
    end
    xyz_sym = xyz_sym(assignment, :);
    %太远的匹配点删掉
    dist = xyz_sym-xyz;
    dist = sqrt(dot(dist', dist')); %1*N
    idx = find(dist > diameter*0.08);
    
    % plot
    figure();
    tmp = find(dist < diameter*0.04);
    for ii=tmp
        pts = [xyz(ii,:); xyz_sym(ii, :)];
        plot3(pts(:,1), pts(:,2), pts(:,3)); hold on;
    end
    axis equal;
    
    xyz_mid = (xyz+xyz_sym)./2;
    xyz_mid(idx) = xyz(idx);
    
    pcshow(xyz_mid, jet(size(xyz, 1)), 'MarkerSize', 50); hold on;
    title('Symmetry');
    hold off;
    
%     figure();pcshow(xyz_mid, jet(size(xyz_mid, 1)), 'MarkerSize', 50);

    %区域生长算法
    [seg, xyz, dir, parts, importantpts, keypts_idx, boundary_idx] = segment(xyz_mid, 5, 0.8, diameter*0.035);


% close all; clear all;
end