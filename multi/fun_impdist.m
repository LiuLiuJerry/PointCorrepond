function f=fun_impdist(X, xyz_pts, xyz_kpts, nbones)
% X : n *(nbones * 3) 2ά����
% xyz_pts:�ָ�õĴ�0��ʼ���ȵ����꣬Ԫ������
n = size(X, 1);
dist_kpts = zeros(n, 1);
dist_models = zeros(n,1);
for i = 1:n %n ��ģ��
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
            % �͹�ͷ�ľ���
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