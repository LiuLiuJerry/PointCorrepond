formatSpec = 'D:\\GitHub\\data\\SkeletonsbyJerry\\%s\\airplane_0%d.2048.ply';
filename = 'airplane';
ii=686
    
path = sprintf(formatSpec,filename, ii);
ptCloud2 = pcread(path);
xyz2 = ptCloud2.Location;
xyz2 = sortrows(xyz2, 1);
ptCloud = pointCloud(xyz2);
    
writePath = path;
pcwrite(ptCloud, writePath);
