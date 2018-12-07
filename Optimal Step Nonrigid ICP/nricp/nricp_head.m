function [ vertsTransformed, X, assignment ] = nricp_head( Source, Target, Options)
% nricp performs an adaptive stiffness variant of non rigid ICP.
%
% 按照关键点分割后的，带有法向的点之间的配准，多了一些约束条件，效果好像不怎么样
% 约束：关键点的要匹配，权重比较大，鼓励方向相同的配准， 分割好的部分有刚性的约束，即倾向于不变
%
% Inputs:
%   Source: structured object with fields - 
%       Source : V x 3 vertices of template model
%       vertices: vertexes 注意必须是double
%       Normals
% 
%   Target : stuctured object as above for target model.
% 
%   Options : structured object with fields:
%       gamm : real valued, weights differences in the rotational and skew
%           part of the deformation against the translational part.
%       epsilon : real values, tolerence for change in transformation.
%       lambda : If using the bi-directional distance metric this weights
%           the contribution of the target -> source term.
%       alphaSet : decreasing vector of real-valued stiffness parameters. 
%           High stiffness parameters force global transformations whereas 
%           low values allow for local deformations.
%       beta : force the deformed shape to have same stiffness as target
%       biDirectional : logical, specifies that a bi-directional distance 
%           is used.
%       plot : logical, specifies that the transformations should be
%           plotted.
%       rigidInit : logical, specifies that rigid ICP should be performed
%           first before allowing non-rigid and non-global deformations.
%       useHungarian : one-to-one alignment. when you use this as search
%           meethod, you don't need to set biDirectional or stiffFromTarget as
%           true. use knn when it is set 0
%       useAuction : another one-to-one alignment method, which is faster
%
% Outputs:
%   vertsTransformed : N X 3 vertices of transformed source mesh,
%   X : (4N) X 3 stacked matrix of transformations.
%   assignment : 1 * N Matrix of corresponded target index

% Set default parameters
if ~isfield(Options, 'gamm')
    Options.gamm = 1;
end
if ~isfield(Options, 'epsilon')
    Options.epsilon = 1e-3;
end
if ~isfield(Options, 'lambda')
    Options.lambda = 10;
end
if ~isfield(Options, 'alphaSet')
    Options.alphaSet = linspace(1.5, 1, 1);%stiffness parameters.平均分两份 
end
if ~isfield(Options, 'beta')
    Options.betaSet = linspace(1, 2, 2);
end
if ~isfield(Options, 'stiffFromTarget')
    Options.stiffFromTarget = 0;
end
if ~isfield(Options, 'useHungarian ')
    Options.useHungarian  = 0;
end
if ~isfield(Options, 'useAuction ')
    Options.useAuction  = 0;
end
if ~isfield(Options, 'biDirectional')
    Options.biDirectional = 1;
end
if ~isfield(Options, 'plot')
    Options.plot = 0;
end
if ~isfield(Options, 'rigidInit')
    Options.rigidInit = 1;
end
if ~isfield(Options, 'SearchRadius')
    Options.SearchRadius = 0.035;
end
if ~isfield(Options, 'onlyKeyPoints')
    Options.onlyKeyPoints = 1;
end
if ~isfield(Options, 'dirWeighting')
    Options.dirWeighting = 0;
end
% Get source vertices 
vertsSource = Source.vertices;
nVertsSource = size(vertsSource, 1);

% Get target vertices
vertsTarget = Target.vertices;
nVertsTarget = size(vertsTarget, 1);

dirSource = Source.normals;
dirTarget = Target.normals;

% Optionally plot source and target surfaces
if Options.plot == 1
    figure();
    clf;
%     p = scatter3(vertsTarget(:,1), vertsTarget(:,2), vertsTarget(:,3), 36, jet(nVertsTarget), '+');
    p = scatter3(vertsTarget(:,1), vertsTarget(:,2), vertsTarget(:,3), 36, 'r', '+');
    hold on;

    h = scatter3(vertsSource(:,1),vertsSource(:,2), vertsSource(:,3), 36, 'b', 'filled');

    material dull; light; grid on; xlabel('x'); ylabel('y'); zlabel('z');
    view([60,30]); axis equal; axis manual;
    legend('Target', 'Source', 'Location', 'best')
    drawnow;
end

% Get subset of target vertices if Options.biDirectional == 1
if Options.biDirectional == 1
    %samplesTarget = sampleVerts(Target, 15);
    samplesTarget = vertsTarget;
    nSamplesTarget = size(samplesTarget, 1);
end
%% stiffness 
% Set matrix G (equation (3) in Amberg et al.) 
% 在旋转和倾斜部分的变形对抗平移部分的变形
G = diag([1 1 1 Options.gamm]);

% Set incidence matrix M 
%source上邻居点距离相近的约束
ptCloud = pointCloud(vertsSource);
bbox = [ptCloud.XLimits; ptCloud.YLimits; ptCloud.ZLimits];
rs = bbox(:,1)-bbox(:,2);
diameter_S = sqrt(dot(rs,rs));
DistS = pdist2(vertsSource, vertsSource);
r = diameter_S*0.04;
A = DistS < r;

% 图关联矩阵 如果 s 和 t 是 G 中第 j 条边的源和目标节点的节点 ID，则 I(s,j) = -1 且 I(t,j) = 1。即 I 的每一列指示 G 中单条边的源和目标节点。
M = adjacency2incidence(A)';

%target的边对应的source距离近的约束
if Options.stiffFromTarget == 1
    DistT = pdist2(vertsTarget, vertsTarget);
    GraT = graph(DistT);
    [Tt,pred] = minspantree(GraT);
    At = adjacency(Tt);
end

% Precompute kronecker product of M and G
kron_M_G = kron(M, G);%M中每个值都变成4*4的小矩阵
%%
% 稀疏矩阵原矩阵的大小为 （nVertsSource, 4*nVertsSource)
% Set matrix D (equation (8) in Amberg et al.)
I = (1:nVertsSource)';
J = 4*I;
D = sparse([I;I;I;I],[J-3;J-2;J-1;J],[vertsSource(:);ones(nVertsSource,1)],nVertsSource, 4*nVertsSource);

% Set weights vector
wVec = ones(nVertsSource,1);

% Set target points matrix tarU and target weights matrix tarU
% if Options.biDirectional == 1.
if Options.biDirectional == 1
    tarU = samplesTarget;
    tarW = eye(nSamplesTarget);
end

%% Do rigid iterative closest point if Options.rigidInit == 1
if Options.rigidInit == 1
    disp('* Performing rigid ICP...');
    % logical: true or false
     [R, t] = icp(vertsTarget', vertsSource', 50, 'Verbose', true, 'Matching', 'kDtree');

    X = repmat([R'; t'], nVertsSource, 1);
    vertsTransformed = D*X;
    
    % Update plot
    if Options.plot == 1
        set(h, 'XData', vertsTransformed(:,1),'YData', vertsTransformed(:,2),'ZData', vertsTransformed(:,3));
        drawnow;
    end
else
    % Otherwise initialize transformation matrix X with identity matrices
    X = repmat([eye(3); [0 0 0]], nVertsSource, 1);
end
%%
% get number of element in the set of stiffness parameters Options.alphaSet
nAlpha = numel(Options.alphaSet);

% Enter outer loop of the non-rigid iterative closest point algorithm. The
% outer loop iterates over stiffness parameters alpha.
disp('* Performing non-rigid ICP...');
for i = 1:nAlpha
    
    % Update stiffness
    alpha = Options.alphaSet(i);
    beta = Options.betaSet(i);
    
 
    % set oldX to be very different to X so that norm(X - oldX) is large on 
	% first iteration
	oldX = 10*X;
    
    % Enter inner loop. For each stiffness setting alternate between 
    % updating correspondences and getting optimal transformations X. 
    % Break the loop when consecutive transformations are similar.
    while norm(X - oldX) >= Options.epsilon 
        
        % Transform source points by current transformation matrix X
        vertsTransformed = D*X;
        
        % Update plot 
        if Options.plot == 1
            set(h, 'XData', vertsTransformed(:,1),'YData', vertsTransformed(:,2),'ZData', vertsTransformed(:,3));
            drawnow;
        end
        
        % Determine closest points on target U to transformed source points
       %% To adapt the cost function to fixed correspondences, only the first term has to be changed. 

        if i == nAlpha && Options.useAuction == 1
            if nVertsSource ~= nVertsTarget
                disp('Sorry! the number of S and T must be same in Auction!');
                return;
            end
            dis = pdist2(vertsTransformed, vertsTarget);
            disMax = -dis + max(max(dis));
            min(min(disMax));
            costmat = double(disMax*100);    
            targetId = sparseAssignmentProblemAuctionAlgorithm(costmat);    
            disp('arrived Auction here');
%         elseif i == nAlpha && Options.useHungarian == 1 %最后一次迭代的时候才使用
        elseif Options.useHungarian == 1 
            dis = pdist2(vertsTransformed, vertsTarget);
            targetId = munkres(dis);

        else
            targetId = knnsearch(vertsTarget, vertsTransformed);
           
        end
%         uncorresponded1 = nUncorrPnts(targetId)
        U = vertsTarget(targetId,:);
        kron_M_G = alpha .* kron_M_G;
        % Update weight matrix  通过获取 B 的列并沿 d 指定的对角线放置它们，来创建一个 m×n 稀疏矩阵
        % Optionally transform surface normals to compare with target and
        % give zero weight if surface and transformed normals do not have
        % similar angles.
        if Options.dirWeighting == 1
            corNormalsTarget = dirTarget(targetId,:);
            angle = abs(dot(dirSource, corNormalsTarget, 2));
            wVec = wVec .* (angle>0.5);
        end
        W = spdiags(wVec, 0, nVertsSource, nVertsSource);

        % Get closest points on source tarD to target samples samplesTarget
        if Options.biDirectional == 1
            transformedId = knnsearch(vertsTransformed, samplesTarget);
            uncorresponded2 = nUncorrPnts(transformedId)
            tarD = sparse(nSamplesTarget, 4 * nVertsSource);
            for j = 1:nSamplesTarget
                cor = transformedId(j); % cor是source对应的索引，应该和X的顺序对应起来
                tarD(j,(4 * cor-3):(4 * cor)) = [vertsSource(cor,:) 1];
            end
           %% target stiffness
           if Options.stiffFromTarget == 1
                [row,col] = find(At); %找到非零数值的坐标和值组成三元组
                nedge = size(row, 1);
                Mt2s = sparse(size(row,1), size(At, 1));
                for j = 1:nedge
                    cor1 = transformedId(row(j, 1));
                    cor2 = transformedId(col(j, 1));
                    if(cor1 ~= cor2)
                        Mt2s(j, cor1) = 1;
                        Mt2s(j, cor2) = -1;
                    end
                end
                kron_Mt2s_G = kron(Mt2s, G);
                kron_M_G = cat(1, kron_M_G, beta .*kron_Mt2s_G); %cat(1,A,B) is the same as [A;B]
           end
        end
       %% When correspondences are fixed, the cost function becomes a sparse quadratic system which can be minimised exactly
        % Specify B and C (See equation (12) from paper)
        %点和点之间的距离约束，平移后和target的距离约束
        A = [...
            kron_M_G; 
            W * D;
            ];
        B = [...
            zeros(size(kron_M_G,1), 3);
            W * U;
            ];
        
        % Concatentate additional terms if Options.biDirectional == 1.
        if Options.biDirectional == 1
            A = [...
                A;
                Options.lambda .* tarW * tarD
                ];
            B = [...
                B;
                Options.lambda .* tarW * tarU
                ];
        end

        % Get optimal transformation X and remember old transformation oldX
        oldX = X;
        X = (A' * A) \ (A' * B);
    end
end
%sort(targetId)
% Compute transformed points 
vertsTransformed = D*X;
assignment = targetId;

%% plot key points
% Update plot and remove target mesh
if Options.plot == 1
    set(h, 'XData', vertsTransformed(:,1),'YData', vertsTransformed(:,2),'ZData', vertsTransformed(:,3));
    drawnow;
    pause(2);
%     delete(p);
end

function n = nUncorrPnts(Idx)
    n1 = size(Idx, 1);
    n2 = size(unique(Idx), 1);
    n = n1-n2;
end
end


