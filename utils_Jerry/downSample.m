formatSpec = 'D:\\pointCloud\\%s\\airplane_0%d.ply';
ii=686
Color = jet(2048*0.125);
for ii = 627:726
    path = sprintf(formatSpec,'ort_skeleton_2048', ii);
    name = sprintf('airplane_0%d,ply', ii);
    ptCloudIn = pcread(path);

    ptCloud = pcdownsample(ptCloudIn,'random',0.125) 

    xyz1 = ptCloud.Location;
    xyz1=sortrows(xyz1, 3); xyz1=sortrows(xyz1, 1); xyz1=sortrows(xyz1, 2);
    ptCloudOut = pointCloud(xyz1);
    
    figure(ii-626)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    pcshow(ptCloudOut.Location, Color);
    title(name);
    writePath = sprintf(formatSpec,'ort_skeleton_256', ii);
    pcwrite(ptCloudOut, writePath);

end