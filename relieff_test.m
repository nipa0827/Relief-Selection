function [weight] = relieff_test(X,Y,xa,K,varargin)
%RELIEFF Importance of attributes (predictors) using ReliefF algorithm.
%   [RANKED,WEIGHT] = RELIEFF(X,Y,K) computes ranks and weights of
%   attributes (predictors) for input data matrix X and response vector Y
%   using ReliefF algorithm for classification or RReliefF for regression
%   with K nearest neighbors. For classification, RELIEFF uses K nearest
%   neighbors per class. RANKED are indices of columns in X ordered by
%   attribute importance, meaning RANKED(1) is the index of the most
%   important predictor. WEIGHT are attribute weights ranging from -1 to 1
%   with large positive weights assigned to important attributes.
%
%   If Y is numeric, RELIEFF by default performs RReliefF analysis for
%   regression. If Y is categorical, logical, a character array, or a cell
%   array of strings, RELIEFF by default performs ReliefF analysis for
%   classification.
%
%   Attribute ranks and weights computed by RELIEFF usually depend on K. If
%   you set K to 1, the estimates computed by RELIEFF can be unreliable for
%   noisy data. If you set K to a value comparable with the number of
%   observations (rows) in X, RELIEFF can fail to find important
%   attributes. You can start with K=10 and investigate the stability and
%   reliability of RELIEFF ranks and weights for various values of K.
%
%   [RANKED,WEIGHT] = RELIEFF(X,Y,K,'PARAM1',val1,'PARAM2',val2,...)
%   specifies optional parameter name/value pairs:
%
%       'method'         - Either 'regression' (default if Y is numeric) or
%                          'classification' (default if Y is not numeric).
%       'prior'          - Prior probabilities for each class, specified as
%                          a string ('empirical' or 'uniform') or as a
%                          vector (one value for each distinct group name)
%                          or as a structure S with two fields:  S.group
%                          containing the group names as a categorical
%                          variable, character array, or cell array of
%                          strings; and S.prob containing a vector of
%                          corresponding probabilities. If the input value
%                          is 'empirical' (default), class probabilities
%                          are determined from class frequencies in Y. If
%                          the input value is 'uniform', all class
%                          probabilities are set equal.
%       'updates'        - Number of observations to select at random for
%                          computing the weight of every attribute. By
%                          default all observations are used.
%       'categoricalx'   - 'on' or 'off', 'off' by default. If 'on', treat
%                          all predictors in X as categorical. If 'off',
%                          treat all predictors in X as numerical. You
%                          cannot mix numerical and categorical predictors.
%       'sigma'          - Distance scaling factor. For observation I,
%                          influence on the attribute weight from its
%                          nearest neighbor J is multiplied by
%                          exp((-rank(I,J)/SIGMA)^2), where rank(I,J) is
%                          the position of J in the list of nearest
%                          neighbors of I sorted by distance in the
%                          ascending order. Default is Inf (all nearest
%                          neighbors have the same influence) for
%                          classification and 50 for regression.
%
%   Example:
%       % Identify important predictors in the ionosphere dataset:
%       load ionosphere;
%       [ranked,weights] = relieff(X,Y,10);
%       bar(weights(ranked));
%       xlabel('Predictor rank');
%       ylabel('Predictor importance weight');
%
% See also KNNSEARCH, PDIST2.

%   Copyright 2010-2013 The MathWorks, Inc.
%

% Check number of required arguments
if nargin<3
    error(message('stats:relieff:TooFewInputs'));
end

% Check if the predictors in X are of the right type
if ~isnumeric(X)
    error(message('stats:relieff:BadX'));
end

% Parse input arguments
validArgs = {'method' 'prior' 'updates'   'categoricalx' 'sigma'};
defaults  = {      ''      []     'all'            'off'      []};

% Get optional args
[method,prior,numUpdates,categoricalX,sigma] ...
    = internal.stats.parseArgs(validArgs,defaults,varargin{:});

% Classification or regression?
isRegression = [];
if ~isempty(method)
    method = internal.stats.getParamVal(method,{'regression' 'classification'},'''Method''');
    isRegression = strcmp(method,'regression');
end

% Check the type of Y
if isnumeric(Y)
    if isempty(isRegression)
        isRegression = true;
    end
elseif iscellstr(Y) || ischar(Y) || isa(Y,'categorical') || islogical(Y)
    if     isempty(isRegression)
        isRegression = false;
    elseif isRegression
        error(message('stats:relieff:BadYTypeForClass'));
    end
else
    error(message('stats:relieff:BadYType'));
end

% Reject prior for regression
if isRegression && ~isempty(prior)
    error(message('stats:relieff:NoPriorForRegression'));
end

% Check if the input sizes are consistent
if (~ischar(Y) && length(Y)~=size(X,1)) ...
        || (ischar(Y) && size(Y,1)~=size(X,1))
    error(message('stats:relieff:XYSizeMismatch'));
end

% Prepare data for classification or regression
if isRegression
    [X,Y] = removeNaNs(X,Y);
    
else % Group Y for classification. Get class counts and probabilities.
    % Get groups and matrix of class counts
    if isa(Y,'categorical')
        Y = removecats(Y);
    end
    yy = Y;
    [Y,grp] = grp2idx(Y);
    [X,Y] = removeNaNs(X,Y);
    Ngrp = numel(grp);
    N = size(X,1);
    C = false(N,Ngrp);
    C(sub2ind([N Ngrp],(1:N)',Y)) = true;
    
    % Get class probs
    if isempty(prior) || strcmpi(prior,'empirical')
        classProb = sum(C,1);
    elseif strcmpi(prior,'uniform')
        classProb = ones(1,Ngrp);
    elseif isstruct(prior)
        if ~isfield(prior,'group') || ~isfield(prior,'prob')
            error(message('stats:relieff:PriorWithMissingField'));
        end
        if iscell(prior.group)
            usrgrp = prior.group;
        else
            usrgrp = cellstr(prior.group);
        end
        [tf,pos] = ismember(grp,usrgrp);
        if any(~tf)
            error(message('stats:relieff:PriorWithClassNotFound', grp{ find( ~tf, 1 ) }));
        end
        classProb = prior.prob(pos);
    elseif isnumeric(prior)
        if ~isfloat(prior) || length(prior)~=Ngrp || any(prior<0) || all(prior==0)
            error(message('stats:relieff:BadNumericPrior', Ngrp));
        end
        classProb = prior;
    else
        error(message('stats:relieff:BadPrior'));
    end
    
    % Normalize class probs
    classProb = classProb/sum(classProb);
    
    % If there are classes with zero probs, remove them
    zeroprob = classProb==0;
    if any(zeroprob)
        t = zeroprob(Y);
        if sum(t)==length(Y)
            error(message('stats:relieff:ZeroWeightPrior'));
        end
        Y(t) = [];
        X(t,:) = [];
        C(t,:) = [];
        C(:,zeroprob) = [];
        classProb(zeroprob) = [];
    end
end

% Do we have enough observations?
if length(Y)<2
    error(message('stats:relieff:NotEnoughObs'));
end

% Check the number of nearest neighbors
if ~isnumeric(K) || ~isscalar(K) || K<=0
    error(message('stats:relieff:BadK'));
end
K = ceil(K);

% Check number of updates
if (~ischar(numUpdates) || ~strcmpi(numUpdates,'all')) && ...
        (~isnumeric(numUpdates) || ~isscalar(numUpdates) || numUpdates<=0)
    error(message('stats:relieff:BadNumUpdates'));
end
if ischar(numUpdates)
    numUpdates = size(X,1);
else
    numUpdates = ceil(numUpdates);
end

% Check the type of X
if ~ischar(categoricalX) || ...
        (~strcmpi(categoricalX,'on') && ~strcmpi(categoricalX,'off'))
    error(message('stats:relieff:BadCategoricalX'));
end
categoricalX = strcmpi(categoricalX,'on');

% Check sigma
if ~isempty(sigma) && ...
        (~isnumeric(sigma) || ~isscalar(sigma) || sigma<=0)
    error(message('stats:relieff:BadSigma'));
end
if isempty(sigma)
    if isRegression
        sigma = 50;
    else
        sigma = Inf;
    end
end

% The # updates cannot be more than the # observations
numUpdates = min(numUpdates, size(X,1));

% Choose the distance function depending upon the categoricalX
if ~categoricalX
    distFcn = 'cityblock';
else
    distFcn = 'hamming';
end

% Find max and min for every predictor
p = size(X,2);
Xmax = max(X);
Xmin = min(X);
Xdiff = Xmax-Xmin;

Fmax = max(xa);
Fmin = min(xa);
Fdiff = Fmax - Fmin;

% Exclude single-valued attributes
isOneValue = Xdiff < eps(Xmax);
if all(isOneValue)
    ranked = 1:p;
    weight = NaN(1,p);
    return;
end
X(:,isOneValue) = [];
Xdiff(isOneValue) = [];
rejected = find(isOneValue);
accepted = find(~isOneValue);

% Scale and center the attributes
if ~categoricalX
    X = bsxfun(@rdivide,bsxfun(@minus,X,mean(X)),Xdiff);
    xa = bsxfun(@rdivide,bsxfun(@minus,xa,mean(xa)),Fdiff);
end

% Get appropriate distance function in one dimension.
% thisx must be a row-vector for one observation.
% x can have more than one row.
if ~categoricalX
    dist1D = @(thisx,x) cityblock(thisx,x);
else
    dist1D = @(thisx,x) hamming(thisx,x);
end

% Call ReliefF. By default all weights are set to NaN.
weight = NaN(1,p);
if ~isRegression
    weight = RelieffClass(X,xa,yy,C,classProb,numUpdates,K,distFcn,dist1D,sigma);
    % weight = weight/size(X,2);
else
    weight =   RelieffReg(X,Y,          numUpdates,K,distFcn,dist1D,sigma);
end



% -------------------------------------------------------------------------
function attrWeights = RelieffClass(scaledX,xa,yy,C,classProb,numUpdates,K,...
    distFcn,dist1D,sigma)
% ReliefF for classification

[numObs,numAttr] = size(scaledX);
attrWeights = 0;
Nlev = size(C,2);

% Choose the random instances
rndIdx = randsample(numObs,numUpdates);
idxVec = (1:numObs)';

% Make searcher objects, one object per class.
searchers = cell(Nlev,1);
for c=1:Nlev
    searchers{c} = createns(xa(C(:,c),:),'Distance',distFcn);
end

classH = zeros(size(unique(yy),1),1);
classM = zeros(size(unique(yy),1),1);
% Outer loop, for updating attribute weights iteratively
for i = 1:numUpdates
    thisObs = rndIdx(i);
    
    % Choose the correct random observation
    selectedX = xa(thisObs,:);
    % Find the class for this observation
    thisC = C(thisObs,:);
    
    % Find the k-nearest hits
    sameClassIdx = idxVec(C(:,thisC));
    
    % we may not always find numNeighbor Hits
    lenHits = min(length(sameClassIdx)-1,K);
    
    % find nearest hits
    % It is not guaranteed that the first hit is the same as thisObs. Since
    % they have the same class, it does not matter. If we add observation
    % weights in the future, we will need here something similar to what we
    % do in ReliefReg.
    Hits = [];
    if lenHits>0
        idxH = knnsearch(searchers{thisC},selectedX,'K',lenHits+1);
        idxH(1) = [];
        Hits = sameClassIdx(idxH);
    end
    
    [neighbour,D_neigh]= knnsearch(xa,selectedX,'K',lenHits+1,'Distance','cityblock');
    neighbour(1) = [];
    D_neigh(1) = [];
    
    zz = ismember(neighbour, sameClassIdx);
    
    rval = sum(zz)/7;
    cls = find(thisC);
    
    % Process misses
    missClass = find(~thisC);
    Misses = [];
    
    if ~isempty(missClass) % Make sure there are misses!
        % Find the k-nearest misses Misses(C,:) for each class C ~= class(selectedX)
        % Misses will be of size (no. of classes -1)x(K)
        Misses = zeros(Nlev-1,min(numObs,K+1)); % last column has class index
        
        for mi = 1:length(missClass)
            
            % find all observations of this miss class
            missClassIdx = idxVec(C(:,missClass(mi)));
            
            % we may not always find K misses
            lenMiss = min(length(missClassIdx),K);
            
            % find nearest misses
            idxM = knnsearch(searchers{missClass(mi)},selectedX,'K',lenMiss);
            Misses(mi,1:lenMiss) = missClassIdx(idxM);
            
        end
        
        % Misses contains obs indices for miss classes, sorted by dist.
        Misses(:,end) = missClass;
    end
    
    %***************** ATTRIBUTE UPDATE *****************************
    % Inner loop to update weights for each attribute
    
    dH = diffH(scaledX,thisObs,Hits,dist1D,sigma)/numUpdates;
    dM = diffM(scaledX,thisObs,Misses,dist1D,sigma,classProb)/numUpdates;
    attrWeights = attrWeights - dH + dM;
    
    %     tt = find(thisC==1);
    %
    %
    %         classH(tt) = classH(tt) + diffH(scaledX,thisObs,Hits,dist1D,sigma,classProb(cls));
    %         classM(tt) = classM(tt) + diffM(scaledX,thisObs,Misses,dist1D,sigma,classProb);
    
end

% c = arrayfun(@(x)length(find(yy == x)), unique(yy), 'Uniform', false);
% gc = cell2mat(c);
%
% attrWeights = -mean(classH./gc) + mean(classM./gc);

%Helper functions for RelieffReg and RelieffClass

%--------------------------------------------------------------------------
% DIFFH (for RelieffClass): Function to calculate difference measure
% for an attribute between the selected instance and its hits

function distMeas = diffH(X,thisObs,Hits,dist1D,sigma,classProb)

% If no hits, return zero by default
if isempty(Hits)
    distMeas = 0;
    return;
end

% Get distance weights
distWts = exp(-((1:length(Hits))/sigma).^2)';
distWts = distWts/sum(distWts);

% Calculate weighted sum of distances
distMeas = sum(dist1D(X(thisObs,:),X(Hits,:)).*distWts);


%--------------------------------------------------------------------------
% DIFFM (for RelieffClass) : Function to calculate difference measure
% for an attribute between the selected instance and its misses
function distMeas = diffM(X,thisObs,Misses,dist1D,sigma,classProb)

distMeas = 0;

% If no misses, return zero
if isempty(Misses)
    return;
end

% Loop over misses
for mi = 1:size(Misses,1)
    
    ismiss = Misses(mi,1:end-1)~=0;
    
    if any(ismiss)
        cls = Misses(mi,end);
        nmiss = sum(ismiss);
        
        distWts = exp(-((1:nmiss)/sigma).^2)';
        distWts = distWts/sum(distWts);
        
        distMeas = distMeas + ...
            sum(dist1D(X(thisObs,:),X(Misses(mi,ismiss),:)).*distWts(1:nmiss)) ...
            *classProb(cls);
    end
end
% Normalize class probabilities.
% This is equivalent to P(C)/(1-P(class(R))) in ReliefF paper.
totProb = sum(classProb(Misses(:,end)));
distMeas = distMeas/totProb;


function [distMeas] = diffT(X,thisObs,Hits,Misses,dist1D,sigma,classProb)

distMeas = 0;

% If no misses, return zero
if isempty(Misses)
    return;
end

D_hit_new = dist1D(X(thisObs,:),X(Hits,:));

% Loop over misses
for mi = 1:size(Misses,1)
    
    ismiss = Misses(mi,1:end-1)~=0;
    
    if any(ismiss)
        cls = Misses(mi,end);
        nmiss = sum(ismiss);
        distWts = exp(-((1:nmiss)/sigma).^2)';
        distWts = distWts/sum(distWts);
        
        D_miss = dist1D(X(thisObs,:),X(Misses(mi,ismiss),:));
        
        %[D_miss_sorted,tt_miss] = sort(D_miss,'ascend');
        
        penalty = find(D_miss<=max(D_hit_new));
        
        
        distMeas = distMeas + ...
            (sum(sqrt(sum((X(Hits,:)-X(Misses(mi,ismiss),:)).^2,2)).*distWts(1:nmiss)) ...
            *classProb(cls));
    end
end

% Normalize class probabilities.
% This is equivalent to P(C)/(1-P(class(R))) in ReliefF paper.
totProb = sum(classProb(Misses(:,end)));
distMeas = distMeas/totProb;


function count = countRval(X,thisObs,Hits,Misses, dist1D)

% If no misses, return zero
if isempty(Misses)
    return;
end

D_hit_new = dist1D(X(thisObs,:),X(Hits,:));
count = 0;

% Loop over misses
for mi = 1:size(Misses,1)
    
    ismiss = Misses(mi,1:end-1)~=0;
    
    if any(ismiss)
        D_miss = dist1D(X(thisObs,:),X(Misses(mi,ismiss),:));
        %[D_miss_sorted,tt_miss] = sort(D_miss,'ascend');
        penalty = find(D_miss<=max(D_hit_new));
        count = count + size(penalty,1)/7;
    end
end

% Normalize class probabilities.
% This is equivalent to P(C)/(1-P(class(R))) in ReliefF paper.
count = count / size(Misses,1);



function [X,Y] = removeNaNs(X,Y)
% Remove observations with missing data
NaNidx = bsxfun(@or,isnan(Y),any(isnan(X),2));
X(NaNidx,:) = [];
Y(NaNidx,:) = [];


function d = cityblock(thisX,X)
B = [];
for i =1:size(X,1)
    B(i,:)= abs(thisX-X(i,:)).^2;
end

B = (sum(B,2));
d = (B).^(1/2);
%d = sqrt(sum((thisX-X).^2,2));

function d = hamming(thisX,X)
d = thisX~=X;
