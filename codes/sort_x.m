formatSpec = 'D:\\pointCloud\\%s\\airplane_0%d.ply';
ii=686
    
path = sprintf(formatSpec,'gt_skeleton', ii);
name = sprintf('airplane_0%d,ply', ii);
ptCloud2 = pcread(path);
xyz2 = ptCloud2.Location;
xyz2 = sortrows(xyz2, 2);
ptCloud = pointCloud(xyz2);
    
writePath = sprintf(formatSpec,'gt_skeleton', 10000);
pcwrite(ptCloud, writePath);
