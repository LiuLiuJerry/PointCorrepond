function f=fun_impdist(X, xyz_pts, xyz_kpts, nbones)
% X : n *(nbones * 3) 2维数组
% xyz_pts:分割好的从0开始的腿的坐标，元胞数组
n = size(X, 1);
dist_kpts = zeros(n, 1);
dist_models = zeros(n,1);
for i = 1:n %n 个模型
    xyz = X(i,:);
    xyz = reshape(xyz, [nbones, 3]);

    kpts = xyz_kpts{i};
    nkpts = size(kpts, 1);    
    [~, dist] = knnsearch(xyz, kpts);
    dist = dist(:,1);
    dist_kpts(i) = sum(dist);   
    
    pts = xyz_pts{i};
    npts = size(pts, 1);
 
    lines = [xyz(1:nbones-1, :), xyz(2:nbones, :)];
    disToBones = zeros(npts, nbones-1);
    for j = 1:nbones-1
        AB = lines(j,4:6)-lines(j,1:3);
        lAB = norm(AB);
        for p = 1:npts
            % 和骨头的距离
            pt = pts(p, :);
            AP = pt-lines(j,1:3);
            BP = pt-lines(j,4:6);            
            lAP = norm(AP);            
            lBP = norm(BP);
            
            r = dot(AP,AB)./lAB./lAB;
            lAC = norm(r*lAB);
            lCP = sqrt(lAP.*lAP - lAC.*lAC);
            if r <= 0
                disToBones(p,j) = lAP;
            elseif r >= 1
                disToBones(p,j) = lBP;
            else 
                disToBones(p,j) = lCP;
            end
        end
    end
    disToBones = min(disToBones, [], 2);
    dist_models(i) = sum(disToBones) +  dist_kpts(i); 
end

f = sum(dist_models);

end