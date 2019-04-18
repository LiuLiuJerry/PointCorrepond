function xyz = symmetry(xyz)
%% 把提取到的骨架变成完全左右对称的骨架
%  即左右侧信息互相补全
%     xyz = ptCloud.Location;
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

end