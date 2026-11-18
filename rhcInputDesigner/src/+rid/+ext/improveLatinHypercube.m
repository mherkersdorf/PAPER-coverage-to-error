function [X,Progress,vectorMinDis,vectorD2,vectorPhiP,vectorPhiPNN] = improveLatinHypercube(varargin)
% improveLatinHypercube calculates an optimal Latin Hypercube Design.
%
% improveLatinHypercube generates an optimal Latin Hypercube design with respect
% to the intended distance criterion. Therefore it swaps design points to
% maximize the distance between the design points.
%
% [X,vectorMinDis,vectorD2,vectorPhiP,Progress] = improveLatinHypercube(varargin)
%
% [~] = improveLatinHypercube - Starts the demo: 30 Points, 2D;
% [~] = improveLatinHypercube(X, options) - Give any latin hypercube plus some options.
% [~] = improveLatinHypercube(#Samples, dim, options) - Give number of samples and number of dimenions plus some options.
%
% INPUT:
%   X:          (N x p) - Initial Latin Hypercube Matrix.
%   #Samples:   (1 x 1) - Number of samples in the Latin Hypercube.
%   dim:        (1 x 1) - Number of dimension.
%   options:  Structure containing the options.
%             'plot' - positive rational number, defines how long the is not cahnged during optimization; 0 means no plot.
%             'distance' - 'euclidean' the only implemented distance measure so far.
%             'normalize'- true or false (default: true) switch on/off nomalization of the Latin Hypercube points between 0 and 1.
%             The next two options define how many worst points should be tested for swapping.
%             WARNING: Only one can be used at a time!!
%             'worstDistance2improve' - natural number, defines how many distance levels should be tested to swap.
%             'worstPoints2improve' - natural number, defines how many worst points should be tested to swap.
%             'simultaneousDimensions' - natural number, defines how many dimensions can be changed simultaneously.
%             'fixedSamples' - Index vector that defines the locked points.
%             'maxCalculationTime' - Scalar defines the maximum calculation time in seconds.
%             'timeSpacing' - Scalar defines the time after which the progress of the current optimization is stored.
%             <0> storage after every imporvement of the minimal distance in the Latin Hypercube.
%             between <0> and <inf> storage after defined time step.
%             <inf> storage of the initial values and and the output values iteration.
%
% OUTPUT:
%   X:            N x p - Finished Latin Hypercube Matrix.
%   Progress:     Structure containing the important values for each m improvement defined by the option 'timeSpacing'.
%                 'X'       (N x p x m) - Finished Latin Hypercube Matrix.
%                 'MinDis'  (1 x m)     - Minimum distance of the Latin Hypercube Matrix.
%                 'D2'      (1 x m)     - Squared unnormalized distances.
%                 'PhiP'    (1 x m)     - Phi_p value (t=2,p=50).
%                 'time'    (1 x m)     - calculation time ins sec.
%   vectorMinDis: 1 x k - Minimum distances for all k point swaps.
%   vectorD2:     1 x k - Squared unnormalized distances for all k point swaps.
%   vectorPhiP:   1 x k - Phi_p value (t=2,p=50) for all k point swaps.
%   vectorPhiPNN: 1 x k - Phi_p value (t=2,p=50) for all k point swaps, only the nearest neighbours.
%
%   Torsten Fischer, 04-March-2015
%   Institute of Mechanics & Automatic Control, University of Siegen, Germany
%   Copyright (c) 2015 by Torsten Fischer


%% check input arguments and initialze the latin hypercube matrix
if isempty(varargin)
    % test mode!
    X = repmat(1:5,2,1)';
    numberOfSamples = 5;
    numberOfDimensions = 2;
    % activate plot
    varargin{1} = 'plot'; varargin{2} = 0.1;
elseif numel(varargin{1}) > 1
    % X matrix is given
    X = varargin{1};
    [numberOfSamples, numberOfDimensions] = size(X);
elseif (numel(varargin{1})) == 1 && (numel(varargin{2}))
    % n and dim are given
    numberOfSamples = varargin{1};
    numberOfDimensions = varargin{2};
    X = repmat(1:numberOfSamples,numberOfDimensions,1)';
    if numberOfDimensions>numberOfSamples
        warning('lhs:inputArguments','Number of dimensions is greater than number of samples!')
    end
else
    error('improveLatinHypercube:varargin','Give a matrix <X> with samples (numberOfSamples x numberOfDimension) as first input argument, or <numberOfSamples> and <numberOfDimension> as two seperate input argument!')
end

%% set defaults and options
options = [];
options.plot = 0;
options.distance = 'euclidean';
options.normalize = true;
options.worstDistance2improve = [];
options.worstPoints2improve = [];
options.simultaneousDimensions = 1;
options.fixedSamples = [];
options.maxCalculationTime = [];
options.timeSpacing = 0;
options.display = false;

% Compare defaults and given options
while ~isempty(varargin)
    % Loop over all fieldnames
    if ischar(varargin{1})
        % If it is a string
        if any(strcmp(fieldnames(options), varargin{1}))
            % If field is specified in the inputs use it instead of defaults
            options.(varargin{1}) = varargin{2};
            % Erase the used entries.
            varargin([1 2]) = [];
        else
            error('improveLatinHypercube:options','%s is an unknown option.',varargin{1})
        end
    else
        % Erase the non-option entries.
        varargin(1) = [];
    end
end

% Check which dimension are allowed to be swapped at a time.
if options.simultaneousDimensions >= numberOfDimensions
    options.simultaneousDimensions = numberOfDimensions-1;
    warning('improveLatinHypercube:options','<simultaneousDimensions> must be smaller than <numberOfDimensions>.\n Set <simultaneousDimensions> to <numberOfDimensions>-1');
elseif options.simultaneousDimensions < 1
    options.simultaneousDimensions = 1;
    warning('improveLatinHypercube:options','<simultaneousDimensions> must be at least 1.\n Set <simultaneousDimensions> to 1.');
end

% specify dimensions which should be testet
dimensionToChange = num2cell(1:numberOfDimensions);
for couterSimultaneousDimensions = 2:options.simultaneousDimensions
    % if more than one dimension should be switched
    dimensionToChange = [dimensionToChange, mat2cell(nchoosek(1:numberOfDimensions,couterSimultaneousDimensions),ones(1,factorial(numberOfDimensions)/(factorial(numberOfDimensions-couterSimultaneousDimensions) * factorial(couterSimultaneousDimensions))),couterSimultaneousDimensions)'];
end

% normalize if necessary
if options.normalize
    minX = min(X);
    maxX = max(X);
    X = bsxfun(@rdivide, bsxfun(@minus, X, min(X)), max(X)-min(X));
end

%% Initialize the optimization
% Initialize the <Progress> structure to store minimum distance improvement.
Progress.X = [];
Progress.MinDis = zeros(1);
Progress.D2 = [];
Progress.PhiP = [];
Progress.time = [];
iterUpdate = 0;

if options.timeSpacing < 0
    error('improveLatinHypercube:options','<timeSpacing> must be between zero and inf!')
end

% Initialize develop flag true to start optimization
develop = true;

% Calculate the nearest neighbours and distances of the initial latin hypercube design.
[NNPointDistance,NNPoints] = pdist2(X,X,options.distance,'Smallest',2);
NNPoints = NNPoints(2,:)';
NNPointDistance = NNPointDistance(2,:)';

% Additional outputs to review the optimisation.
vectorMinDis = [];
vectorD2 = [];
vectorPhiP = [];
vectorPhiPNN = [];


%% Precalculate the number of worst points to consider during the optimization.
% Check if both are empty.
if isempty(options.worstDistance2improve) && isempty(options.worstPoints2improve)
    options.worstDistance2improve = 1;
    % Check if both are given.
elseif ~isempty(options.worstDistance2improve) && ~isempty(options.worstPoints2improve)
    error('improveLatinHypercube:options','Only one can be given: Dicide <worstDistance2improve> or <worstPoints2improve>.')
else
    % Check if the number of distance levels are given and if all points must be considered.
    if isempty(options.worstPoints2improve) && options.worstDistance2improve >= numel(NNPointDistance)
        idxWorstPoints2Calculate = 1:numel(NNPointDistance);
    elseif isempty(options.worstDistance2improve) % Check if the number of worst points are given.
        % Number of must be greater one.
        if options.worstPoints2improve < 1
            error('improveLatinHypercube:options','<worstPoints2improve> must be a natural number.')
        end
        % Check if all points must be considered.
        if options.worstPoints2improve >= numel(NNPointDistance)
            idxWorstPoints2Calculate = 1:numel(NNPointDistance);
        else % Calculate only the demanded number of worst points
            idxWorstPoints2Calculate = 1:options.worstPoints2improve;
        end
    end
end

% Get the starting time, if needed.
startTime = tic;
idxChanged = [];
timeIsOver = false;
plotCounter = 1;
countSwaps = 0;

%% Main optimization loop
while develop
    
    % Update history of the optimization
    if nargout > 5
        [minDis, D2, PhiP] = distributionQuality(X,min(NNPointDistance),options.normalize);
        vectorMinDis = [vectorMinDis minDis];
        vectorD2 = [vectorD2  D2];
        vectorPhiP = [vectorPhiP PhiP];
        PhiPNN = norm(1./NNPointDistance,50);
        vectorPhiPNN = [vectorPhiPNN PhiPNN];
        
    elseif nargout > 4
        [minDis, D2, PhiP] = distributionQuality(X,min(NNPointDistance),options.normalize);
        vectorMinDis = [vectorMinDis minDis];
        vectorD2 = [vectorD2  D2];
        vectorPhiP = [vectorPhiP PhiP];
    elseif nargout > 2
        [minDis, D2] = distributionQuality(X,min(NNPointDistance),options.normalize);
        vectorMinDis = [vectorMinDis minDis];
        vectorD2 = [vectorD2  D2];
    elseif nargout > 1
        [minDis] = distributionQuality(X,min(NNPointDistance),options.normalize);
        vectorMinDis = [vectorMinDis minDis];
    end
    
    % Setting develop to 'false' aborts optimization.
    develop = false;
    
    % Rank the nearest neighbours
    [MinPointDistance,IdxRankPointDistance] = sort(NNPointDistance);
    
    if options.display && ~isempty(idxChanged)
        fprintf('\n\n Swap of %i and %i leads to a minimum distance of %4.4f .',idxChanged(1),idxChanged(2),MinPointDistance(1))
    end
    
    % Clear change index
    idxChanged = [];
    
    % Store history of the distance improvement
    if nargout >1
        if iterUpdate == 0
            iterUpdate = iterUpdate + 1;
            Progress.X(:,:,iterUpdate) = X;
            [Progress.MinDis(iterUpdate), Progress.D2(iterUpdate), Progress.PhiP(iterUpdate)] = distributionQuality(X,MinPointDistance(1),options.normalize);
            Progress.time(iterUpdate) = toc(startTime);
            if options.timeSpacing ~= 0 && options.timeSpacing ~= inf
                timeProgress = tic;
            end
        else
            if options.timeSpacing == 0
                % Store only after each improvement of the critical distance.
                if nargout > 1 && (MinPointDistance(1) - Progress.MinDis(end) > eps)
                    iterUpdate = iterUpdate + 1;
                    Progress.X(:,:,iterUpdate) = X;
                    [Progress.MinDis(iterUpdate), Progress.D2(iterUpdate), Progress.PhiP(iterUpdate)] = distributionQuality(X,MinPointDistance(1),options.normalize);
                    Progress.time(iterUpdate) = toc(startTime);
                end
            elseif options.timeSpacing ~= inf && toc(timeProgress) > options.timeSpacing
                % Store at each intended time step.
                iterUpdate = iterUpdate + 1;
                Progress.X(:,:,iterUpdate) = X;
                [Progress.MinDis(iterUpdate), Progress.D2(iterUpdate), Progress.PhiP(iterUpdate)] = distributionQuality(X,MinPointDistance(1),options.normalize);
                Progress.time(iterUpdate) = toc(startTime);
                timeProgress = tic;
            end
        end
    end
    
    % Check the time, if intended.
    if ~isempty(options.maxCalculationTime)
        % Get the current time
        currentTime = toc(startTime);
        if options.maxCalculationTime < currentTime
            if options.display
                fprintf('\n\n Maximum calcualtion time of %4.4f is reached: %4.4f.',options.maxCalculationTime,currentTime)
            end
            % Time is over abort optimization.
            timeIsOver = true;
            break
        end
    end
    
    if ~isempty(options.fixedSamples)
        idxFixedSamples = any(cell2mat(arrayfun(@(x) IdxRankPointDistance == x , options.fixedSamples,'uniformOutput',false)),2);
        MinPointDistance(idxFixedSamples) = [];
        IdxRankPointDistance(idxFixedSamples) = [];
        idxWorstPoints2Calculate = 1:numel(MinPointDistance);
%         idxWorstPoints2Calculate(idxFixedSamples) = [];
    end
    
    % Check if the number of distance levels is given.
    if isempty(options.worstPoints2improve) && options.worstDistance2improve < numel(NNPointDistance)
        % Check for natural number.
        if options.worstDistance2improve < 1
            error('improveLatinHypercube:options','<worstDistance2improve> must be a natural number.')
        end
        % Get the points with equal nearest neighbour distance
        idxWorstPoints2Calculate= [];
        MinPointDistanceCalc = MinPointDistance;
        IdxRestCalc = 1:length(MinPointDistance);
        % Get all the points of the considered distance levels.
        for counterDistanceLevel = 1:options.worstDistance2improve
            idxCurrentDistanceLevel = find(abs(MinPointDistanceCalc - MinPointDistanceCalc(1)) < eps);
            idxWorstPoints2Calculate = [idxWorstPoints2Calculate IdxRestCalc(idxCurrentDistanceLevel)];
            MinPointDistanceCalc(idxCurrentDistanceLevel) =[];
            IdxRestCalc(idxCurrentDistanceLevel) = [];
            if isempty(MinPointDistanceCalc)
                break
            end
        end
    end
    
    % Loop over all worst points
    for idxWorstPoint = idxWorstPoints2Calculate
        
        % Get the current worst point.
        MinNNPoint = IdxRankPointDistance(idxWorstPoint);
        MinNNPartner = NNPoints(MinNNPoint);
        
        if options.plot
            figure(99)
            clf
            if size(X,2) == 2
                plot(X(:,1),X(:,2),'.','Markersize',25)
                hold on
                plot(X(MinNNPoint,1),X(MinNNPoint,2),'r.','Markersize',25)
                plot(X(MinNNPartner,1),X(MinNNPartner,2),'ro','Markersize',25)
                hold off
                for kk = 1:size(X,1)
                    text(X(kk,1)+0.025, X(kk,2)+0.025,num2str(kk))
                end
            else
                plot3(X(:,1),X(:,2),X(:,3),'.','Markersize',25)
                hold on
                plot3(X(MinNNPoint,1),X(MinNNPoint,2),X(MinNNPoint,3),'r.','Markersize',25)
                plot(X(MinNNPartner,1),X(MinNNPartner,2),X(MinNNPartner,3),'ro','Markersize',25)
                hold off
            end
            if options.normalize
                set(gca,'XTick',linspace(0,1,numberOfSamples),'YTick',linspace(0,1,numberOfSamples))
            else
                set(gca,'XTick',linspace(1,numberOfSamples,numberOfSamples),'YTick',linspace(1,numberOfSamples,numberOfSamples))
            end
            grid on
        end
        
        % Those points who have a coordinate between the worst case pair are bad trading partners.
        maxWorstCase = max(X([MinNNPoint MinNNPartner],:));
        minWorstCase = min(X([MinNNPoint MinNNPartner],:));
        
        for idxCurrentDimension = 1:numberOfDimensions % loop over all inputs
            
            % Specify which coordinate of which dimension should be switched.
            currentDimensionToChange = dimensionToChange{idxCurrentDimension};
            if numel(currentDimensionToChange) > 1
                keyboard
            end
            
            % Specify which cycling partner should be used.
            partner2try = 1:numberOfSamples;
            
            % Neither the NN nor the cases that have already been tested are possible trading partners.
            partner2try([NNPoints(MinNNPoint), IdxRankPointDistance(1:idxWorstPoint)' options.fixedSamples]) = [];
            
            % If no trading partner can be found, abort calculation.
            if isempty(partner2try);
                develop = false;
                break
            end
            
            for idxCurrentSampleCounter = 1:length(partner2try)% Loop over all possible trading partners.
                
                % specify which cycling partner should be used
                idxCurrentSample = partner2try(idxCurrentSampleCounter);
                
                % Check if the partner is a possible good trading partner.
                if ~all((X(idxCurrentSample,currentDimensionToChange)<=maxWorstCase(currentDimensionToChange)) & (X(idxCurrentSample,currentDimensionToChange)>=minWorstCase(currentDimensionToChange)))
                    
                    % switch LH levels
                    XTemp = X;
                    XTemp(idxCurrentSample,currentDimensionToChange) = X(MinNNPoint,currentDimensionToChange);
                    XTemp(MinNNPoint,currentDimensionToChange) = X(idxCurrentSample,currentDimensionToChange);
                    
                    % Calculate the distance of the points in the range of the swapped worst point.
                    idxTryMinNNPoint = [1:MinNNPoint-1  MinNNPoint+1:size(XTemp,1)];
                    distance2AllPointsMinNNPoint = sqrt(sum(bsxfun(@minus,XTemp(idxTryMinNNPoint,:),XTemp(MinNNPoint,:)).^2,2));
                    [NNChangedMinNNPoint,Index] = min(distance2AllPointsMinNNPoint);
                    NNChangeIdxMinNNPoint = idxTryMinNNPoint(Index);
                    
                    % Check if the distance of the swapped worst point is greater then the current worst distance.
                    if NNChangedMinNNPoint - MinPointDistance(idxWorstPoint) > eps
                        
                        % Calculate the distance of the points in the range of the swapped current samples point.
                        idxTryCurrentSample = [1:idxCurrentSample-1  idxCurrentSample+1:size(XTemp,1)];
                        distance2AllPointsCurrentSample = sqrt(sum(bsxfun(@minus,XTemp(idxTryCurrentSample,:),XTemp(idxCurrentSample,:)).^2,2));
                        [NNChangedCurrentSample, Index] = min(distance2AllPointsCurrentSample);
                        NNChangeIdxCurrentSample = idxTryCurrentSample(Index);
                        
                        % Check if the distance of the swapped current sample point is greater then the current worst distance.
                        if NNChangedCurrentSample - MinPointDistance(idxWorstPoint) > eps
                            % Update the lhs design.
                            develop = true;
                            idxChanged = [MinNNPoint idxCurrentSample];
                            
                            if options.plot
                                figure(99)
                                if size(X,2) == 2
                                    hold on
                                    plot(X(idxChanged,1),X(idxChanged,2),'go','Markersize',25)
                                    hold off
                                else
                                    hold on
                                    plot3(X(idxChanged,1),X(idxChanged,2),X(idxChanged,3),'go','Markersize',25)
                                    hold off
                                end
                                if options.normalize
                                    set(gca,'XTick',linspace(0,1,numberOfSamples),'YTick',linspace(0,1,numberOfSamples))
                                else
                                    set(gca,'XTick',linspace(1,numberOfSamples,numberOfSamples),'YTick',linspace(1,numberOfSamples,numberOfSamples))
                                end
                                grid on
                            end
                            
                            % Calculate the nearest neigbour of the neighbour points of the swapped samples.
                            [NNPointsCalcNew, NNPointDistanceCalcNew, IdxCalcNew] = calculateNNPoints(XTemp,NNPoints,options.distance,MinNNPoint,idxCurrentSample);
                            
                            % Update the sample matrix.
                            X = XTemp;
                            countSwaps = countSwaps +1;
                            
                            % Update the nearest neighbour relations.
                            NNPoints([MinNNPoint idxCurrentSample IdxCalcNew]) = [NNChangeIdxMinNNPoint NNChangeIdxCurrentSample NNPointsCalcNew]';
                            NNPointDistance([MinNNPoint idxCurrentSample IdxCalcNew]) = [NNChangedMinNNPoint NNChangedCurrentSample NNPointDistanceCalcNew]';
                            
                            % Check if the distances of the changed points to all other points are smaller than their current nearest neighbours, than do the necessary update.
                            idxUpdateMinNNPoint = find(NNPointDistance(idxTryMinNNPoint) - distance2AllPointsMinNNPoint >= eps,1);
                            if ~isempty(idxUpdateMinNNPoint)
                                NNPoints(idxTryMinNNPoint(idxUpdateMinNNPoint)) = MinNNPoint;
                                NNPointDistance(idxTryMinNNPoint(idxUpdateMinNNPoint)) = distance2AllPointsMinNNPoint(idxUpdateMinNNPoint);
                            end
                            
                            idxUpdateCurrentSample = find(NNPointDistance(idxTryCurrentSample) - distance2AllPointsCurrentSample >= eps,1);
                            if ~isempty(idxUpdateCurrentSample)
                                NNPoints(idxTryCurrentSample(idxUpdateCurrentSample)) = idxCurrentSample;
                                NNPointDistance(idxTryCurrentSample(idxUpdateCurrentSample)) = distance2AllPointsCurrentSample(idxUpdateCurrentSample);
                            end
                            
                            % permute the dimension vector.
                            dimensionToChange = [dimensionToChange(2:end) dimensionToChange(1)];
                            
                            break
                        end
                    end
                end
            end
            if develop
                break
            end
        end
        if develop
            break
        end
    end
    
    if ~isempty(idxChanged)
        if options.plot
            figure(99)
            if size(X,2) == 2
                hold on
                plot(X(idxChanged,1),X(idxChanged,2),'g.','Markersize',25)
                hold off
            else
                hold on
                plot3(X(idxChanged,1),X(idxChanged,2),X(idxChanged,3),'g.','Markersize',25)
                hold off
            end
            if options.normalize
                set(gca,'XTick',linspace(0,1,numberOfSamples),'YTick',linspace(0,1,numberOfSamples))
            else
                set(gca,'XTick',linspace(1,numberOfSamples,numberOfSamples),'YTick',linspace(1,numberOfSamples,numberOfSamples))
            end
            grid on
            [minDis, D2] = distributionQuality(X,min(NNPointDistance),options.normalize);
            title(['minDis ' num2str(minDis) '; D^2: ' num2str(D2)],'Fontsize',12)
            %             pause(options.plot)
            %             saveas(gcf,[num2str(plotCounter) '_DetOptMinDis.fig'])
            %             plotCounter = plotCounter + 1;
            drawnow
        end
    end
    
end

if options.display
    if timeIsOver
        fprintf('\n\n Optimization is finished. No calculation time left.  Distance reached: %4.4f.\n\n',min(NNPointDistance))
    else
        fprintf('\n\n Optimization is finished. No further swaps possible. Distance reached: %4.4f.\n\n',min(NNPointDistance))
    end
end

if options.plot
    figure(99)
    if options.normalize
        set(gca,'XTick',linspace(0,1,numberOfSamples),'YTick',linspace(0,1,numberOfSamples))
    else
        set(gca,'XTick',linspace(1,numberOfSamples,numberOfSamples),'YTick',linspace(1,numberOfSamples,numberOfSamples))
    end
    grid on
    [minDis, D2] = distributionQuality(X,min(NNPointDistance),options.normalize);
    title(['minDis ' num2str(minDis) '; D^2: ' num2str(D2)],'Fontsize',12)
    drawnow
end

if nargout > 1
    % Store at each intended time step.
    iterUpdate = iterUpdate + 1;
    Progress.X(:,:,iterUpdate) = X;
    [Progress.MinDis(iterUpdate), Progress.D2(iterUpdate), Progress.PhiP(iterUpdate)] = distributionQuality(X,MinPointDistance(1),options.normalize);
    Progress.time(iterUpdate) = toc(startTime);
end

% normalize if necessary
if options.normalize
    X = bsxfun(@plus, bsxfun(@times, X, (maxX-minX)), minX);
end

end

function [NNPointsNew,NNPointDistanceNew,IdxCalcNew] = calculateNNPoints(X,NNPoints,typeDistance,Point1,Point2)
% Updates the distance and index vector only for the necessary points.
if nargin < 4
    IdxCalcNew = 1:numel(NNPoints);
else
    % Just update the changed points and thier nearest neighbors
    IdxCalcNew = [find(NNPoints == Point1 | NNPoints == Point2)]';
end

% Get real number idicies
IdxCalcNew = unique(IdxCalcNew);

% Calculate all modified nearest neighbour relationships
NNPointsNew = NaN(size(IdxCalcNew));
NNPointDistanceNew = NaN(size(IdxCalcNew));

% Loop over all points to calculate.
for k = 1:size(IdxCalcNew,2)
    switch typeDistance
        case 'euclidean'
            Distance = sum(bsxfun(@minus,X([1:IdxCalcNew(k)-1 IdxCalcNew(k)+1:size(X,1)],:),X(IdxCalcNew(k),:)).^2,2);
        otherwise
            error('Only <euclidean> is implemented so far.')
    end
    
    [Value,Index] = min(Distance);
    
    if Index>=IdxCalcNew(k)
        % Add one for the correct index.
        Index = Index+1;
    end
    % Update nearest neighbour relation.
    NNPointsNew(k) = Index;
    NNPointDistanceNew(k) = sqrt(Value);
end
end