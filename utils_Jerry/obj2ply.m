close all;

formatSpec = 'D:\\GitHub\\data\\%s\\mesh_%.4d%s';
formatSpec2 = 'D:\\GitHub\\data\\%s\\mesh_%.4d%s';
input_file = 'octopus\octopus_aligned';
output_file = 'octopus\octopus_aligned_ply';
% input_file = 'human\I_squat\squat_sample'
% output_file = 'human\I_squat\squat_2048'
% input_file = 'human\human_aligned';
% output_file = 'human\human_aligned_ply';

extension1 = '.obj';
extension2 = '.ply';
idx = [0:49];
start = 0;
%0-49 50-62 63-113
for ii = idx
    path = sprintf(formatSpec, input_file, ii, extension1);
    [pts, ~] = read_mesh(path);
    xyz0 = pts;
    if size(pts,2) > size(pts,1)
        xyz0 = pts';
    end
    n = size(xyz0,1);
    k = 2048;
  
%     for mesh_0000.ply octopus
    alphaX = -pi/4;
    alphaY = -pi;
    transformX = [1, 0, 0; 0, cos(alphaX), sin(alphaX); 0, -sin(alphaX), cos(alphaX)];
    transformY = [cos(alphaY), 0, sin(alphaY); 0, 1, 0; -sin(alphaY),  0, cos(alphaY)];
    tmp = transformX * xyz0';
    xyz0 = tmp';

% % %     for Cartoon style octopus
%     alphaX = 0;
%     alphaY = -pi-1/10*pi;
%     alphaZ = -pi;
%     transformX = [1, 0, 0; 0, cos(alphaX), sin(alphaX); 0, -sin(alphaX), cos(alphaX)];
%     transformY = [cos(alphaY), 0, sin(alphaY); 0, 1, 0; -sin(alphaY),  0, cos(alphaY)];
%     tmp = transformY * xyz0';
%     xyz0 = tmp';

    % 归一化
    ptCloud = pointCloud(xyz0);
    figure()
    pcshow(ptCloud); %hold on;
    x1=xlabel('X轴');        
    x2=ylabel('Y轴');      
    x3=zlabel('Z轴');  
    
    xyz = xyz0;
    Limits = [ptCloud.XLimits; ptCloud.YLimits; ptCloud.ZLimits];
    max_len = norm(Limits(:,2)-Limits(:,1));
    center = (Limits(:,2)+Limits(:,1))./2;
    center(3,:) = 0;
    xyz = xyz-center';
    xyz = xyz./max_len; %斜对角长度是1
    ptCloud = pointCloud(xyz);
    
%     figure()
%     pcshow(ptCloud); %hold on;
    pathout = sprintf(formatSpec2, output_file, ii+start, extension2);
    pcwrite(ptCloud, pathout);
end