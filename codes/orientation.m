function []=orientation(formatSpec, filename, outfile, N)
% 遍历找到对称的方向然后旋转到x轴
% 若y<0 的方向是机尾，则再旋转180°
% close all;
% clear all;
% formatSpec = 'D:\\pointCloud\\%s\\airplane_0%d.ply';
for ii=679:679
    
    path = sprintf(formatSpec, filename, ii);
    name = sprintf('airplane_0%d,ply', ii);
    ptCloud2 = pcread(path);
    
    Color = jet(N);
    %init symmerty
    xyz2 = ptCloud2.Location;
    xyz2(:, 1) = xyz2(:, 1)*-1;
    [D, I] = pdist2(ptCloud2.Location, xyz2,'euclidean', 'Smallest', 1);
    %D = pdist2(ptCloud2.Location, xyz2,'euclidean');
    mindis = sum(D);
    xyz = ptCloud2.Location;
    xyz_ = xyz2;
    %找对称轴
    rot = 32;
    pii = pi/rot;
    for jj = 1:1:rot
        [ptCloudTformed, xyz2] = rotate(pii*jj, ptCloud2);
        %rad = 180/rot*jj
        len = sum(ptCloudTformed.XLimits);
        xyz2(:, 1) = xyz2(:, 1)-len/2;
        xyz_2 = xyz2;
        xyz_2(:, 1) = xyz2(:, 1)*-1;

        [D, I] = pdist2(xyz2, xyz_2, 'euclidean', 'Smallest', 1);
        %D = pdist2(xyz_2, xyz2,'euclidean');
        dis = sum(D);
        if(dis < mindis)
            mindis = dis;
            xyz = xyz2;
            xyz_ = xyz_2;
        end
    end
    
    %[tform, ptCloudTformed] = pcregistericp(ptCloudTformed, ptCloud1, "InlierRatio", 1);
    ptCloud2 = pointCloud(xyz_);
    ptCloud1 = pointCloud(xyz);
     
    xyz1 = ptCloud1.Location;
    xyz2 = ptCloud2.Location;
    xyz1=sortrows(xyz1, 3); xyz1=sortrows(xyz1, 1); xyz1=sortrows(xyz1, 2);
    xyz2=sortrows(xyz2, 3); xyz2=sortrows(xyz2, 1); xyz2=sortrows(xyz2, 2);
    ptCloud1 = pointCloud(xyz1);
    ptCloud2 = pointCloud(xyz2);
    
    %检查是不是尾巴在前面
    n = N*0.1;
    tail1 = xyz2(1:n, :);
    tail2 = xyz2(end-n:end, :);
    v1 = sum( var(tail1(:, 1)) );
    v2 = sum( var(tail2(:, 1)) );
    if(v1 > v2)
        [ptCloud2, xyz] = rotate(pi, ptCloud2);
    end
    figure(ii-626)
    %pcshow(ptCloud1.Location, Color);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    pcshow(ptCloud2.Location, Color);
    title(name);
    hold off
    writePath = sprintf(formatSpec, outfile, ii);
    pcwrite(ptCloud2, writePath);
end


%% function for rotate
function [ptCloud, xyz2]=rotate(rad, ptCloud2)
        A = [cos(rad) sin(rad) 0 0; ...
        -sin(rad) cos(rad) 0 0; ...
                0         0  1 0; ...
                0         0  0 1];
        tform1 = affine3d(A);
        ptCloudTformed = pctransform(ptCloud2,tform1);
        xyz2 = ptCloudTformed.Location;
        ptCloud = pointCloud(xyz2);
end

end
