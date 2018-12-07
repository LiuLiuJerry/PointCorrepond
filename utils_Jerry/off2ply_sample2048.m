close all;
clear all;

formatSpec = 'D:\\GitHub\\data\\%s\\mesh_%d%s';
formatSpec2 = 'D:\\GitHub\\data\\%s\\mesh_%.4d%s';
% input_file = 'Octopus\octopus_sample';
% output_file = 'Octopus\octopus2048';
input_file = 'Chen\horse\horse_sample4096';
output_file = 'Chen\horse\horse4096';

extension1 = '.ply';
extension2 = '.ply';
idx = [0:47];
start = 0
%0-49 50-62 63-113
for ii = idx
    path = sprintf(formatSpec, input_file, ii, extension1)
    [pts, ~] = read_mesh(path);
    xyz0 = pts;
    if size(pts,2) > size(pts,1)
        xyz0 = pts';
    end
    n = size(xyz0,1);
    k = 4096;

%     if ii > 82
%         if ii == 83 
%             alphaZ = pi/2-20./180*pi;
%             transformZ = [cos(alphaZ), sin(alphaZ), 0; -sin(alphaZ), cos(alphaZ), 0; 0, 0, 1];
%             tmp = transformZ * xyz0';
%             xyz0 = tmp';
%         end
%         if ii == 86|| ii==85
%             alphaZ = -10./180*pi;
%             transformZ = [cos(alphaZ), sin(alphaZ), 0; -sin(alphaZ), cos(alphaZ), 0; 0, 0, 1];
%             tmp = transformZ * xyz0';
%             xyz0 = tmp';
%         end
%         alphaX = pi/2;
%         transformX = [1, 0, 0; 0, cos(alphaX), sin(alphaX); 0, -sin(alphaX), cos(alphaX)];
%         tmp = transformX * xyz0';
%         xyz0 = tmp';
%     end
%     if ii >= 47 && ii <= 82
%         alphaX = -pi/4;
%         transformX = [1, 0, 0; 0, cos(alphaX), sin(alphaX); 0, -sin(alphaX), cos(alphaX)];
%         tmp = transformX * xyz0';
%         xyz0 = tmp';        
%     end
%     if ii< 47
%         alphaY = -pi-20./180*pi;
%         transformY = [cos(alphaY), 0, sin(alphaY); 0, 1, 0; -sin(alphaY),  0, cos(alphaY)];
%         tmp = transformY * xyz0';
%         xyz0 = tmp'; 
% 
%     end

    
    % ¹éÒ»»¯
    ptCloud = pointCloud(xyz0);
    
    xyz = xyz0(randperm(n, k), :);
    Limits = [ptCloud.XLimits; ptCloud.YLimits; ptCloud.ZLimits];
    max_len = norm(Limits(:,2)-Limits(:,1));
    center = (Limits(:,2)+Limits(:,1))./2;
%     center(3,:) = 0;
    xyz = xyz-center';
    xyz = xyz./max_len;
    ptCloud = pointCloud(xyz);
    
    figure();
    pcshow(ptCloud.Location); %hold on;
    x1=xlabel('XÖá');        
    x2=ylabel('YÖá');      
    x3=zlabel('ZÖá');  
    pathout = sprintf(formatSpec2, output_file, ii+start, extension2);
    pcwrite(ptCloud, pathout);
end