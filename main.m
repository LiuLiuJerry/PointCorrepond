clear all;
close all;
%% sort, oriented
formatSpec = 'D:\\GitHub\\data\\SkeletonsbyJerry\\%s\\airplane_0%d.2048.ply';
filename = 'airplane';
filename_ort = 'airplane_ort';
filename_ali = 'airplane_aligned';
N = 512;
% orientation(formatSpec, filename, filename_ort, N);

ICPSkeleton(formatSpec, filename_ort, filename_ali);

clear all;
close all;