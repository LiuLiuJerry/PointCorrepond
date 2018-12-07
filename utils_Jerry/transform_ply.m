formatSpec = 'D:\\GitHub\\data\\%s\\mesh_%.4d%s';
input_file = 'Octopus\octopus';
output_file = 'Octopus\octopus2048';
%input_file = 'Chen\samba'
%output_file = 'Chen\samba2048_ply'
extension1 = '.ply';
extension2 = '.ply';
idx = [53];

for ii = idx
    path = sprintf(formatSpec, output_file, ii, extension1);
    ptCloud = pcread(path)
    xyz0 = ptCloud.Location
    n = size(xyz0,1);
    k = 2048;
    
    alphaX = 0;
    alphaY = pi/2;
    transformX = [1, 0, 0; 0, cos(alphaX), sin(alphaX); 0, -sin(alphaX), cos(alphaX)];
    transformY = [cos(alphaY), 0, sin(alphaY); 0, 1, 0; -sin(alphaY),  0, cos(alphaY)];
    tmp = transformY * xyz0';
    xyz0 = tmp';
    
    % πÈ“ªªØ
    ptCloud = pointCloud(float(xyz0));
    figure()
    pcshow(ptCloud); %hold on;
    x1=xlabel('X÷·');        
    x2=ylabel('Y÷·');      
    x3=zlabel('Z÷·');  
    
    pathout = sprintf(formatSpec, output_file, ii, extension2);
     pcwrite(ptCloud, pathout);
end