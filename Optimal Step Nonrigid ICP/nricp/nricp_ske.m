function [ vertsTransformed, X ] = nricp( Source, Target, Options )
% nricp performs an adaptive stiffness variant of non rigid ICP.
%
% This function deforms takes a dense set of landmarks points from a template
% template model and finds a deformation which matches a target shape. 
% The deformations are encouraged to be natural and smooth by means of a 
% stiffness constraint, which is relaxed in increments.
% 
% For details on the stiffness constraint and optimization procedure see: 
% 'Optimal Step Nonrigid ICP Algorithms for Surface Registration', 
% Amberg, Romandhani and Vetter, CVPR, 2007.
%
% Inputs:
%   Source: structured object with fields - 
%       Source.vertices: V x 3 vertices of template model
%       Source.faces: F x 3 list of connected vertices.
%       Source.normals: (Optional) FV x 3 list of surface normals. Make
%           sure to set Options.normals = 1 if using normals.
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
%
% Outputs:
%   vertsTransformed : N X 3 vertices of transformed source mesh,
%   X : (4N) X 3 stacked matrix of transformations.

% Set default parameters
if ~isfield(Options, 'gamm')
    Options.gamm = 1;
end
if ~isfield(Options, 'epsilon')
    Options.epsilon = 1e-4;
end
if ~isfield(Options, 'lambda')
    Options.lambda = 2;
end
if ~isfield(Options, 'alphaSet')
    Options.alphaSet = linspace(10, 1, 1);
end
if ~isfield(Options, 'beta')
    Options.betaSet = linspace(1, 2, 1);
end
if ~isfield(Options, 'stiffFromTarget')
    Options.stiffFromTarget = 0;
end
if ~isfield(Options, 'elasticity')
    Options.elasticity = 0;
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

% Get source vertices 
vertsSource = Source.vertices;
nVertsSource = size(vertsSource, 1);

% Get target vertices
vertsTarget = Target.vertices;


% Optionally plot source and target surfaces
if Options.plot == 1
    clf;
    p = scatter3(vertsTarget(:,1), vertsTarget(:,2), vertsTarget(:,3),'+', 'y');
    hold on;
    
    h = scatter3(vertsSource(:,1),vertsSource(:,2), vertsSource(:,3),'filled', 'b');
    material dull; light; grid on; xlabel('x'); ylabel('y'); zlabel('z');
    %view([60,30]); axis equal; axis manual;
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
DistS = pdist2(vertsSource, vertsSource);
GraS = graph(DistS);
[Ts,pred] = minspantree(GraS);% 最小生成树
% 图关联矩阵 如果 s 和 t 是 G 中第 j 条边的源和目标节点的节点 ID，则 I(s,j) = -1 且 I(t,j) = 1。即 I 的每一列指示 G 中单条边的源和目标节点。
M = incidence(Ts)';

%target的边对应的source距离近的约束
if Options.stiffFromTarget == 1
    DistT = pdist2(vertsTarget, vertsTarget);
    GraT = graph(DistT);
    [Tt,pred] = minspantree(GraT);
    At = adjacency(Tt);
end

% Precompute kronecker product of M and G
kron_M_G = kron(M, G);
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
        % 将target按照最近点的距离重新排序 最近点可能重复
        targetId = knnsearch(vertsTarget, vertsTransformed);
        nUncorrPnts(targetId)
        U = vertsTarget(targetId,:);
        kron_M_G = alpha .* kron_M_G;
        % Update weight matrix  通过获取 B 的列并沿 d 指定的对角线放置它们，来创建一个 m×n 稀疏矩阵
        % 可以用来强调或删除某些点的影响
        W = spdiags(wVec, 0, nVertsSource, nVertsSource);

        % Get closest points on source tarD to target samples samplesTarget
        if Options.biDirectional == 1
            transformedId = knnsearch(vertsTransformed, samplesTarget);
            nUncorrPnts(transformedId)
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
                kron_M_G = cat(1, kron_M_G, beta .*kron_Mt2s_G);
           end
        end
       %% When correspondences are fixed, the cost function becomes a sparse quadratic system which can be minimised exactly
        % Specify B and C (See equation (12) from paper)
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

% Update plot and remove target mesh
if Options.plot == 1
    set(h, 'XData', vertsTransformed(:,1),'YData', vertsTransformed(:,2),'ZData', vertsTransformed(:,3));
    drawnow;
    pause(2);
%     delete(p);
end

function [ samples ] = sampleVerts( Mesh, radius )
% sampleVerts sub samples the vertices of a mesh. Vertices are selected 
% so that no other nodes lie within a pre-determined radius.
% 
% Inputs:
%   Mesh : structured object with fields:
%                   Mesh.vertices: N x 3 vertices of Mesh.
%                   Mesh.faces: M x 3 list of connected vertices.
%   radius : controls the spacing of the vertices.
    samples = [];
    vertsLeft = Mesh.vertices;
    itt = 1;
while size(vertsLeft, 1) > 0
        nVertsLeft = size(vertsLeft, 1);

        % pick a sample from remaining points
        vertN = randsample(nVertsLeft, 1);
        vert = vertsLeft(vertN, :);

        % Add it to sample set
        samples(itt,:) = vert;

        % Remove nearby vertices
        idx = rangesearch(vertsLeft, vert, radius);
        idRemove = idx{1};
        vertsLeft(idRemove, :) = [];

        itt = itt + 1;
end


function n = nUncorrPnts(Idx)
    n1 = size(Idx, 1);
    n2 = size(unique(Idx), 1);
    n = n1-n2;


