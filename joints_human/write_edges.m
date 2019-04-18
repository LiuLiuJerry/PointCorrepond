%% write edges down
% 定义骨架的结构
n_imptpts = 2;
nbones_h = 1;
nbones_leg = 3;
nbones_arm = 2;
n_pts = n_imptpts + nbones_h + nbones_leg*2 + nbones_arm*2 ;
cnt_parts = 6;

% 写的顺序： xyz = [kpts; xyz_heads; xyz_legsf_l; xyz_legsf_r; xyz_legsb_l; xyz_legsb_r; ];
nsamples = [n_imptpts, nbones_h, nbones_arm, nbones_arm, nbones_leg, nbones_leg];

M = zeros(n_pts, n_pts);
n2 = 0;
n1 = 0;
for j = 1:cnt_parts
        n1 = n2+1;
        n2 = n2+nsamples(j);
        for i = n1:n2-1
            M(i,i+1) = 1;
        end
end
M(2,3) = 1;
M(2,4) = 1;
M(2,6) = 1;
M(1,8) = 1;
M(1,11) = 1;

path_txt = "D:/GitHub/data/human/I_squat/squat_kpts/edges.txt";
[m,n] = find(M>0);

% %test
% ptCloud = pcread('D:/GitHub/data/human/I_squat/squat_kpts/mesh_0002.ply');
% pts = ptCloud.Location;
% E = [m,n];
% figure();
% for i = 1:size(E,1)
%     e = E(i,:);
%     plot3(pts(e,1), pts(e,2), pts(e,3));hold on;
% end

fp=fopen(path_txt,'a');%'A.txt'为文件名；'a'为打开方式：在打开的文件末端添加数据，若文件不存在则创建。
for i = 1:size(E,1)
    e = E(i,:);
    fprintf(fp,'%d %d\n',e(1), e(2));%fp为文件句柄，指定要写入数据的文件。注意：%d后有空格。
end
fclose(fp);%关闭文件。