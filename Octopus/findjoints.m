function [f] = findjoints(X, legs_idx, xyzs)
%%
% legs: 各个腿上的采样点，按顺序排列, 列向量
% A: 三段骨节的长度
nlegs = size(legs_idx, 1);
p = zeros(1, 3);
cost = zeros(nlegs, 1);
nbones = 4;
pos = zeros(nbones+1, 3);
for i = 1:nlegs
    xyz = xyzs{floor((i+1)./2)}; %两条腿一起
    leg = legs_idx{i};
    leg = xyz(leg',:);
    pos(1,:) = leg(1,:);
    dir = zeros(4, 3);
    for j = 1:nbones %三节骨头    
        %各个模型的骨架长度一样
        leglen = leg - pos(j,:);
        leglen = sqrt(dot(leglen, leglen, 2));
        tmp = find(leglen > X(j));
        if isempty(tmp)
            pts = [pos(j, :); leg];
        else
            pts = [pos(j, :); leg(1:tmp(1)-1, :)];
        end
        
        if size(pts, 1) == 1 %使用和上一个点连线的方向
            n = dir(j-1, :);
        elseif size(pts, 1) == 2
            n = pts(2,:) - pts(1,:);
        else
            %用得到的点拟合直线       
            x = pts(:,1); % x = az+b  A = [a,b]' F = [Z,1];
            y = pts(:,2);  % y = bz+d  B = [c,d]'
            z = pts(:,3);
            F = [z, ones(size(z))];
            M = F'*F; % F'FA = F'X; M*A = N
            N = F'*x; % F'FB = F'Y; M*B = O;
            O = F'*y;
            A=(M\N); %A = [a,b],
            B=(M\O); %B = [c,d]'; 
    %         x1=A(1)*z+A(2);
    %         y1=B(1)*z+B(2);
    %         z1=z;        
            n = [A(1)\1, B(1)\1, 1];
        
        end
        n = n./norm(n);
        dir(j, :) = n;
        p = X(j).*n;
        pos(j+1, :) = pos(j, :) + p; % x.*n = X(j)
        leg = leg(tmp(1):end, :);    
    end
    %计算三条直线之间的两个夹角
    ang = dot(dir(2:nbones, :), dir(1:nbones-1, :), 2);
    cost(i) = sum(ang);
    %拟合误差一定小的情况下， 夹角最大，或者cos最大
end

f = sum(cost);
end 
