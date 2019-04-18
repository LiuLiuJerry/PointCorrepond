close all;
clear all;

formatSpec = 'D:\\GitHub\\data\\%s\\mesh_%d%s';
formatSpec2 = 'D:\\GitHub\\data\\%s\\mesh_%.4d%s';
input_file = 'multi\lioness\lioness_sample';
output_file = 'multi\lioness\lioness_Normalized';
% input_file = 'multi\Cats\cat_sample';
% output_file = 'multi\Cats\cat_Normalized';
% input_file = 'human\I_squat\squat_sample'
% output_file = 'human\I_squat\squat_2048'
% input_file = 'human\human_aligned';
% output_file = 'human\human_aligned_ply';

extension1 = '.ply';
extension2 = '.ply';

% read data 1 
path = sprintf(formatSpec, input_file, 0, extension1);
pcResp = pcread(path);
Limits = [pcResp.XLimits; pcResp.YLimits; pcResp.ZLimits];
max_len = norm(Limits(:,2)-Limits(:,1));
center = (Limits(:,2)+Limits(:,1))./2;

idx = [0:9];
start = 0;
for ii = idx
    path = sprintf(formatSpec, input_file, ii, extension1);
%     [pts, ~] = read_mesh(path);
    ptCloud = pcread(path);
    xyz0 = ptCloud.Location;
    if size(xyz0,2) > size(xyz0,1)
        xyz0 = xyz0';
    end
    n = randperm(size(xyz0,1)); %随机数
    k = 2048;
    idx_pts = n(1:k);
    xyz0 = xyz0(idx_pts, :);
  
% %     for mesh_0000.ply octopus
%     alphaX = -pi/4;
%     alphaY = -pi;
%     transformX = [1, 0, 0; 0, cos(alphaX), sin(alphaX); 0, -sin(alphaX), cos(alphaX)];
%     transformY = [cos(alphaY), 0, sin(alphaY); 0, 1, 0; -sin(alphaY),  0, cos(alphaY)];
%     tmp = transformX * xyz0';
%     xyz0 = tmp';

% % %     for Cartoon style octopus
%     alphaX = 0;
%     alphaY = -pi-1/10*pi;
%     alphaZ = -pi;
%     transformX = [1, 0, 0; 0, cos(alphaX), sin(alphaX); 0, -sin(alphaX), cos(alphaX)];
%     transformY = [cos(alphaY), 0, sin(alphaY); 0, 1, 0; -sin(alphaY),  0, cos(alphaY)];
%     tmp = transformY * xyz0';
%     xyz0 = tmp';

    % 归一化 
    xyz = xyz0;

% % read data 1 每个模型单独做归一化 
%     Limits = [ptCloud.XLimits; ptCloud.YLimits; ptCloud.ZLimits];
%     max_len = norm(Limits(:,2)-Limits(:,1));
%     center = (Limits(:,2)+Limits(:,1))./2;
%     xyz = xyz-center';
%     xyz = xyz./max_len; %斜对角长度是1
%     xyz = [xyz(:,1)  xyz(:,3)  -xyz(:,2)];
%     ptCloud = pointCloud(xyz);

% read data 1 同第一个模型对齐做归一化
    xyz = xyz-center';
    xyz = xyz./max_len; %斜对角长度是1
    xyz = [xyz(:,1)  xyz(:,3)  -xyz(:,2)];
    ptCloud = pointCloud(xyz);
    
    figure()
    pcshow(ptCloud); %hold on;
    x1=xlabel('X轴');        
    x2=ylabel('Y轴');      
    x3=zlabel('Z轴');  
    
    pathout = sprintf(formatSpec2, output_file, ii+start, extension2);
    pcwrite(ptCloud, pathout);
end