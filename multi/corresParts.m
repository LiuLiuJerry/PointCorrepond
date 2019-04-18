function []=corresParts(Source, Target,  key_idxS, key_idxT, partsS, partsT)
%同时沿着最长的地方缩放

% 刚性变换
pidxS = Source(key_idxS, :);
pidxT = Target(key_idxT, :);
ptS = pointCloud(pidxS);
ptT = pointCloud(pidxT);
tform = pcregistericp(ptS, ptT);
ptCloudAligned = pctransform(pointCloud(Source),tform);
xyzs = ptCloudAligned.Location;
xyzt = Target;
figure();
    for ii = partsS'
        plot3(xyzs(ii,1), xyzs(ii,2), xyzs(ii,3)); hold on;
    end
    for ii = partsT'
        plot3(xyzt(ii,1), xyzt(ii,2), xyzt(ii,3)); hold on;
    end
% key points
plot3(xyzs(key_idxS, 1), xyzs(key_idxS, 2), xyzs(key_idxS, 3), '.', 'MarkerSize', 30); hold on;
% key points
plot3(xyzt(key_idxT, 1), xyzt(key_idxT, 2), xyzt(key_idxT, 3), '.', 'MarkerSize', 30 ); hold on;
axis equal;
% knn 找到备选集合
% 整个图中的位置
% 包围盒（归一化）的距离
% 计算边的向量夹角，若节点连了多条边，则找最近的匹配
end