function []=corresParts(Source, Target,  key_idxS, key_idxT, partsS, partsT)
%ͬʱ������ĵط�����

% ���Ա任
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
% knn �ҵ���ѡ����
% ����ͼ�е�λ��
% ��Χ�У���һ�����ľ���
% ����ߵ������нǣ����ڵ����˶����ߣ����������ƥ��
end