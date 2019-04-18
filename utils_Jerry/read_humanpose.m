clear all;
close all;
data_path = 'D:/I_squat/I_squat/poses/template_poses.txt';
path_out = 'D:/I_squat/I_squat/poses/ske1.ply';
n = 250;
k = 42;
data = importdata(data_path);
data = reshape(data, n,k);
data1 = data(1,:);

pose = reshape(data1, 3, 14);
pcshow(pose', 'MarkerSize', 50);
pcwrite(pointCloud(pose'), path_out);
