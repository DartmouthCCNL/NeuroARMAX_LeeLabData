function [vars, stat, ret] = Fitting_Function(minFunc, inpData, varsNum, iterNum, crossPer, crossIter, posInit, varsRule, noZSfields, display)
% A general porpouse fitting function.
% It accepts any general function and dataset plus number of variables, restrictions for
% fitting and initial values and returns solution and stats.
% 
%   [USAGE]
%   [vars, stat, ret] = Fitting_Function(minFunc, inpData, varsNum, iterNum, crossPer, crossIter, posInit, varsRule, noZSfields)
%
%   [INPUT]
%   minFunc:    Handler to the cost function for finding the solution.
%               Fitting_Function calls minFun(vars, inpData) to calculate
%               the final value. Considering vars as solution and inpData
%               as the data needed for cost function. 
%
%   inpData:    Input dataset needed to cal minFunc. Struct contain
%               different fields each contain Total_Data_Sample rows.
%               Should contain an 'output' field
%
%   varsNum:    Number of variables in solution.
%
%   iterNum:    Number of iteration (initial points) for running a
%               cross-validation instance.
%
%   crossPer:   Percentage of test data considred for each cross-validation
%               instance. [0:1]
% 
%   crossIter:  Number of cross-validation instances.
%
%   posInit:    Logical vector with the size of varsNum, indicating the
%               variables which need to be initiated with a positive value.
%               e.g. [false true false false, true, false]
%
%   varsRule:   A cell structure with the size of varsNum, contaning
%               minimum and maximum values of each variable. Empty cells
%               indicate no limit. e.g. {[], [], [-2 4], [], [0 1], []}
%
%   noZSfields: A cell contains the label of cells which does not need
%               zscoring. Empty cell means nothing. e.g. {'timeReward'}
%
%   display:    Display iterations for cross-validation process 
%               (true, false)
%               
%
%
%   [OUTPUT]
%   vars:       Final solution.
%   
%   stat:       Containing information about the solution. 
%               stat.all:     Stats of all data
%               stat.train:     Stats of train data
%               stat.test:      Stats of test data
%                       .N:     Number of data (all, test, train)
%                       .K:     Number of variables (all, test, train)
%                       .TSS:   Total sum of squares (all, test, train)
%                       .ESS:   Estimated sum of squares (all, test, train) (R2)
%                       .RSS:   Residual sum of squares (all, test, train)
%                       .AIC:   AIC of the data (all, test, train)
%                       .BIC:   BIC of the data (all, test, train)
%
%                       .P:     P-values of the data (all)
%                       .T:     T-values of the data (all)
%                       .P_test:P-values of running signtest over set of
%                               parameters (all)
%
%                       .fOutputMean:   Mean of all outputs (all, test, train)
%                       .fOutput:       All outputs (all, test, train)
%
%   ret:        Extra return values
%                       .fvalFinal      Best f-value
%                       .fvarsFinal     Best solution
%
%                       .fvalTest       F-values for test data (all instances and iterations)
%                       .fvalTrain      F-values for train data (all instances and iterations)
%                       .fvarsAll       Vars (all instances and iterations)
% 
%                       .idxSortFinal   Best models in each instance
% 
%                       .fvalFinalset   F-values for test data (all instances)
%                       .fvarsFinalset  Vars (all instances)
%
% Copyright 2019 Mehran M. Spitmaan (mehran.m.spitman@gmail.com).
% 

%% Setup optimization function
options             = optimoptions('lsqnonlin', 'Algorithm','trust-region-reflective','MaxIter',5000, 'MaxFunEvals',10000);
% options             = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','MaxIter',5000, 'MaxFunEvals',10000);

options.Display     = 'off';

%% setup upper lower limits
upperLim            = [];
lowerLim            = [];

for cntVar = 1:varsNum
    
    if isempty(varsRule{cntVar})
        lowerLim(cntVar)    = -inf;
        upperLim(cntVar)    = inf;
    else
        lowerLim(cntVar)    = varsRule{cntVar}(1);
        upperLim(cntVar)    = varsRule{cntVar}(2);
    end
    
end

%% Setup cross-validation instances
fldNames            = fieldnames(inpData);
totNum              = size(inpData.(fldNames{1}),1); % Total number of samples
foldNum             = floor(totNum * crossPer); % Number of samples in each fold

idxTemp = [];
for cntIter             = 1:crossIter
    rndtemp             = randperm(totNum);
    idxTemp{cntIter}    = rndtemp(1:foldNum);
end

idxTrain            = []; % Indecies of train dataset
idxTest             = []; % Indecies of test dataset
for cntIter = 1:crossIter    
    idxTrain{cntIter}   = setdiff([1:totNum],idxTemp{cntIter});
    idxTest{cntIter}    = [idxTemp{cntIter}];
end

%% Iteration over cross-validation instances

for cntInstance = 1:crossIter

    for cntF = 1:length(fldNames)
        dataSetBinTrain.(fldNames{cntF})         = logical(zeros(size(inpData.(fldNames{cntF}))));
        dataSetBinTrain.(fldNames{cntF})(idxTrain{cntInstance},:)         = true;
    end      
    
    for cntF = 1:length(fldNames)
        dataSetBinTest.(fldNames{cntF})         = logical(zeros(size(inpData.(fldNames{cntF}))));
        dataSetBinTest.(fldNames{cntF})(idxTest{cntInstance},:)         = true;
    end 
   
        
    % Trianing Part
    for cntIter = 1:iterNum        
        %% Initiating variables
        
        varsTemp                            = rand(1,varsNum)*2-1;               
        
        varsTemp(posInit)                   = ...
            abs(varsTemp(posInit));        
        
        minFuncTrain                        = @(v) minFunc(v,dataSetBinTrain);
        
        [solution,RESNORML,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] =...
            lsqnonlin(minFuncTrain,varsTemp,lowerLim,upperLim,options);        
        
        varsTrain(:,cntIter,cntInstance)    = solution';
        fvalTrain(:,cntIter,cntInstance)    = RESNORML';
        
        
        minFuncTest                         = @(v) minFunc(v,dataSetBinTest);
        
        fvalTest(:,cntIter,cntInstance)     = sse(minFuncTest(solution));
        
        if display
            [cntIter,cntInstance]
        end
    end
    
end

%% calculating best model based on test SSE
[~, idxTemp]        = min(fvalTest,[],2);
idxTemp             = idxTemp(:);

fvalTestBest        = [];
varsTrainBest       = [];
for cntInstance = 1:crossIter
    fvalTestBest(cntInstance)       =  fvalTest(:,idxTemp(cntInstance),cntInstance);
    varsTrainBest(:,cntInstance)    = varsTrain(:,idxTemp(cntInstance),cntInstance);
end

% Calculating median of best values over cross-validation instances
fvalFinal           = median(fvalTestBest);
fvarsFinal          = median(varsTrainBest,2);

% setting return values
ret.fvalFinal       = fvalFinal;        % best f-value
ret.fvarsFinal      = fvarsFinal;       % best solution

ret.fvalTest        = fvalTest;         % f-values for test data (all instances and iterations)
ret.fvalTrain       = fvalTrain;        % f-values for train data (all instances and iterations)
ret.fvarsAll        = varsTrain;        % vars (all instances and iterations)

ret.idxSortFinal    = idxTemp;          % best models in each instance

ret.fvalFinalset    = fvalTestBest;     % f-values for test data (all instances)
ret.fvarsFinalset   = varsTrainBest;    % vars (all instances)

% setting final return value
vars                = fvarsFinal';

%% Calculating STATS for TRAIN dataset
statTemp = [];
for cntInstance = 1:crossIter

    for cntF = 1:length(fldNames)
        dataSetBinTrain.(fldNames{cntF})         = logical(zeros(size(inpData.(fldNames{cntF}))));
        dataSetBinTrain.(fldNames{cntF})(idxTrain{cntInstance},:)         = true;
    end      
    
    minFuncTrain                        = @(v) minFunc(v,dataSetBinTrain);
        
    varsTemp                            = varsTrainBest(:,cntInstance);
    statTemp{cntInstance}.RSS           = minFuncTrain(varsTemp');   
    statTemp{cntInstance}.N             = length(statTemp{cntInstance}.RSS);
    statTemp{cntInstance}.K             = size(varsTrainBest,1);
    
%     statTemp{cntInstance}.fOutputMean   = mean(dataSetBinTrain.(outputField));
%     statTemp{cntInstance}.fOutput       = dataSetBinTrain.(outputField);
    
%     statTemp{cntInstance}.TSS           = sse(statTemp{cntInstance}.fOutput - statTemp{cntInstance}.fOutputMean);
    statTemp{cntInstance}.RSS           = sse(statTemp{cntInstance}.RSS);
%     statTemp{cntInstance}.ESS           = 1-(statTemp{cntInstance}.RSS/statTemp{cntInstance}.TSS);
        
    statTemp{cntInstance}.AIC           = 2*statTemp{cntInstance}.K + statTemp{cntInstance}.N * log(statTemp{cntInstance}.RSS);               

end

stat.train = statTemp;


%% Calculating STATS for TEST dataset
statTemp = [];
for cntInstance = 1:crossIter
    
    for cntF = 1:length(fldNames)
        dataSetBinTest.(fldNames{cntF})         = logical(zeros(size(inpData.(fldNames{cntF}))));
        dataSetBinTest.(fldNames{cntF})(idxTest{cntInstance},:)         = true;
    end  
    
    minFuncTest                        = @(v) minFunc(v,dataSetBinTest);
        
    varsTemp                            = varsTrainBest(:,cntInstance);
    statTemp{cntInstance}.RSS           = minFuncTest(varsTemp');   
    statTemp{cntInstance}.N             = length(statTemp{cntInstance}.RSS);
    statTemp{cntInstance}.K             = size(varsTrainBest,1);
    
%     statTemp{cntInstance}.fOutputMean   = mean(dataSetBinTest.(outputField));
%     statTemp{cntInstance}.fOutput       = dataSetBinTest.(outputField);
    
%     statTemp{cntInstance}.TSS           = sse(statTemp{cntInstance}.fOutput - statTemp{cntInstance}.fOutputMean);
    statTemp{cntInstance}.RSS           = sse(statTemp{cntInstance}.RSS);
%     statTemp{cntInstance}.ESS           = 1-(statTemp{cntInstance}.RSS/statTemp{cntInstance}.TSS);
        
    statTemp{cntInstance}.AIC           = 2*statTemp{cntInstance}.K + statTemp{cntInstance}.N * log(statTemp{cntInstance}.RSS);               

end

stat.test = statTemp;



%% Calculating STATS for ALL dataset
statTemp = [];

for cntF = 1:length(fldNames)
    dataSetBinAll.(fldNames{cntF})      = inpData.(fldNames{cntF});
end


minFuncAll                        = @(v) minFunc(v,dataSetBinAll);

varsTemp                                = fvarsFinal;
statTemp.RSS                            = minFuncAll(varsTemp');
statTemp.N                              = length(statTemp.RSS);
statTemp.K                              = size(fvarsFinal,1);

% statTemp.fOutputMean                    = mean(dataSetBinAll.(outputField));
% statTemp.fOutput                        = dataSetBinAll.(outputField);

% statTemp.TSS                            = sse(statTemp.fOutput - statTemp.fOutputMean);
statTemp.RSS                            = sse(statTemp.RSS);
% statTemp.ESS                            = 1-(statTemp.RSS/statTemp.TSS);

statTemp.AIC                            = 2*statTemp.K + statTemp.N * log(statTemp.RSS);


% Claculate P-vals
for cntVars = 1:size(varsTemp,1)
    varsTTemp = varsTemp;
    
    filterTemp = zeros(size(varsTTemp));
    filterTemp(cntVars) = 1;
    varsTTemp = varsTTemp.*filterTemp;

    rssTempt = minFunc(varsTTemp',dataSetBinAll); 
    nTemp = length(rssTempt);        
    mseTemp(cntVars) = sse(rssTempt)/(nTemp-2);
            
%     tvalT(cntVars) = varsTemp(cntVars) / sqrt(mseTemp(cntVars)/(statTemp.TSS));
%     pvalT(cntVars) = 2*tcdf(-abs(tvalT(cntVars)),nTemp-2);
    
    [pvalTsignTest(cntVars)] = signtest(ret.fvarsFinalset(cntVars,:));    
end

% statTemp.P             = pvalT;
% statTemp.tVal          = tvalT;
statTemp.P_test        = pvalTsignTest;


stat.all = statTemp;


end

