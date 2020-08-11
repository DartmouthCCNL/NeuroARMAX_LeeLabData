classdef ARMAXNeuro < handle
    %   ARMAXNEURO is a class for estimating short, seasonal, and exogenous
    %   (reward and choice) memories plus exogenous selctivity of neurons.
    %   
    %   For costruct an instance use "ARMAXNeuro(neuronData)"
    %
    %   Copyright 2019 Mehran M. Spitmaan (mehran.m.spitman@gmail.com).
    
    properties (SetAccess = private)        
        
        ITI = [];                       % Storage for ITI
        
        timeIntervarlsAbsolute          % Absolute time values for each bin and trial
       
        startPoint = 2;                 % starting point for trials
        
        firingRateMat                   % Total Firing rate for trial and bins
        
        maxOrder = 0;                   % maximum memory order of component
        
        boundries                       % boudries for trials and bins
        
        totalDataPoint                  % Total number of datapoint
        
        firingRateMatMean               % Total Firing rate for bins (neural profile)
        
        firingRateMatMeanArranged       % Total Firing rate for bins (neural profile) rearranged for all data points
                
        fittingModelData                % Final Data needed for fitting the model
        
        display                         % Display messages in output
        
        paramSet                        % Data structure for parameter set and their boundries
        
        fittingSet                      % Data structure for fitting procedure
        
        fittingResults                  % Data structure for fitting results
        
        saveAddress                     % Saving folder address
        
        saveNamePrefix                  % Prefix of saving file name
        
        zscoring = 1;                   % Indicates needs for zscoring the data
        
        %% Data Structure for the Input Data
        neuronData
        % Contain all the input information we want to create the model
        %
        %   neuronData (Struct)
        %       .spikeTime:         All the spike time for entire experiment
        %
        %       .signalTime:        Structure contain absolute time of all
        %                           signal for each trial.
        %                           Name of each field represents name of each
        %                           signal
        %
        %       .signalValue:       Structure contain values of all the possible
        %                           signal for each trial. There might be signals
        %                           with signalTime but without signalValue
        %                           Name of each field represents name of each
        %                           signal
        %
        %       .endOfTrialSignal:  Name of the signal that indicates the end of the
        %                           trial. There might be a time distance
        %                           before/after "endOfTrialSignal" to the real end
        %                           of the trial which can be set by
        %                           "endOfTrialDist". (String)
        %
        %       .endOfTrialDist:    Time distanse from "endOfTrialSignal" that
        %                           indicated the real end of the trial. (msec)
        %
        %       .maxTrialLen:       Maximum time for each trial. Useful for
        %                           arranging data into bins. (msec)
        %                           The real value for trial length would be the
        %                           caculated by
        %                               "(signalTime.(endOfTrialSignal)(currTrial) +
        %                               endOfTrialDist) - ...
        %                               (signalTime.(endOfTrialSignal)(currTrial-1) +
        %                               endOfTrialDist);"
        %
        %       .beginOfTrialSignal:    Name of the signal that indicates the begin of the
        %                               trial. (String)
        %
        %       .binSize:           Size of each bin (msec)
        %
        %       .trialRange:        Range of trials; Total number of trials.
        
        
        %% Data Structure for Model       
        
        component
        % Components contains all 3 possible memory types and exogenous
        % selectivity, each have their on data structure. 
        %
        %   .memShort (Struct for short memory)
        %       .defined:               Whether memory is defined or not
        %                               (Boolean)
        %
        %       .order:                 Indicates the memory order
        %
        %       .effectiveTimeWin:      Time window that this memory is
        %                               effective (refer to timeMask function)
        %
        %       .PACF_transfer:         Whether transfer AR Coeffs into
        %                               PACF (Partial AutoCorrelation Function) domain (Boolean)
        %
        %       .varIdx:                Index of variables in total vriable
        %                               vector
        %
        %       .Data:                  Data stored for this memory
        %               .allARCoeff:    All real values of AR Coeffs
        %               .allPAACFCoeff: All real values of PACF Coeffs
        %               .allTau:        All real values of converted taus
        %               .meanTau:       Mean of converted taus
        %               .maxTau:        Max of converted taus
        %
        %                   .p:         Each field has p-val
        %                   .val:       Each field has real value
        %                   .stat:      Each field has stat (t-val, etc.)
        %
        %        
        %
        %   .memSeason (Struct for seasonal memory)
        %       .defined:               Whether memory is defined or not
        %                               (Boolean)
        %
        %       .order:                 Indicates the memory order
        %
        %       .effectiveTimeWin:      Time window that this memory is
        %                               effective (refer to timeMask function)
        %
        %       .PACF_transfer:         Whether transfer AR Coeffs into
        %                               PACF (Partial AutoCorrelation Function) domain (Boolean)
        %
        %       .varIdx:                Index of variables in total vriable
        %                               vector
        %
        %       .Data:                  Data stored for this memory
        %               .allARCoeff:    All real values of AR Coeffs
        %               .allPAACFCoeff: All real values of PACF Coeffs
        %               .allTau:        All real values of converted taus
        %               .meanTau:       Mean of converted taus
        %               .maxTau:        Max of converted taus
        %
        %                   .p:         Each field has p-val
        %                   .val:       Each field has real value
        %                   .stat:      Each field has stat (t-val, etc.)
        %
        %
        %
        %   .memExo (Struct(s) for Exogenous memory)
        %       .name:                  Name of the exogenous memory
        %
        %           .signalName:        A set contains name(s) of the
        %                               signals involve in this memory       
        %
        %           .interaction:       Whether memory is based on
        %                               interaction of signals
        %
        %           .interactionType:   Type of interactio
        %                               1) simple product
        %                               2) seperate values
        %
        %           .signalTime:        Name of the 'signalTime'
        %                               corresponding for this exo memory                                      
        %
        %           .memFunc:           Type of memory function
        %               .expo:          Whether it is a exponential
        %                               function
        %               .expoOrder:     Order of exponential function
        %
        %           .memOrder:          Indicates the memory order
        %
        %           .effectiveTimeWin:  Time window that this memory is
        %                               effective (refer to timeMask function)
        %
        %           .exoSelEff:         Whether consider the effect of exo
        %                               signal (Boolean)
        %                               We should have similar name exo
        %                               signal in "exoSignal" component
        %
        %
        %
        %
        %   .exoSignal (Struct(s) for Exogenous signals)
        %       .name:                  Name of the exogenous memory
        %
        %           .signalName:        A set contains name(s) of the
        %                               signals involve in this memory       
        %
        %           .interaction:       Whether memory is based on
        %                               interaction of signals        
        %
        %           .interactionType:   Type of interactio
        %                               1) simple product
        %                               2) seperate values
        %
        %           .effectiveTimeWin:  Time window that this memory is
        %                               effective (refer to timeMask function)
        %

        
    end
    
    methods
        % Constructor
        function obj = ARMAXNeuro(neuronData, display, saveAddress, saveNamePrefix)
            %   ARMAXNEURO Construct an instance of this class
            %   You can intriduce "neuronData" based on structure below:
            %            
            
            if nargin <2
                display = 1;
                saveAddress = '~/MyNeuralDataStorage/';
                saveNamePrefix = 'MyNeuralModel';
            elseif nargin <3
                saveAddress = '~/MyNeuralDataStorage/';
                saveNamePrefix = 'MyNeuralModel';
            elseif nargin <4
                saveNamePrefix = 'MyNeuralModel';
            end                
            
            obj.neuronData = neuronData;
            obj.display = display;
            obj.saveAddress = saveAddress;
            obj.saveNamePrefix = saveNamePrefix;

            % Initializing components
            
            % Short Memory
            component.memShort.defined = false;
            component.memShort.order = 0;
            component.memShort.effectiveTimeWin = nan;
            component.memShort.PACF_transfer = false;
            component.memShort.varIdx = [];
            
            lblTemp = {'allARCoeff','allPAACFCoeff','allTau','meanTau','maxTau'};
            for cntLBL = 1:length(lblTemp)
               lblCurr =  lblTemp{cntLBL};
               component.memShort.Data.(lblCurr).p = [];
               component.memShort.Data.(lblCurr).val = [];
               component.memShort.Data.(lblCurr).stat = [];
            end
            
            % Seasonal Memory
            component.memSeason.defined = false;
            component.memSeason.order = 0;
            component.memSeason.effectiveTimeWin = nan;
            component.memSeason.PACF_transfer = false;
            component.memSeason.varIdx = [];
            
            lblTemp = {'allARCoeff','allPAACFCoeff','allTau','meanTau','maxTau'};
            for cntLBL = 1:length(lblTemp)
               lblCurr =  lblTemp{cntLBL};
               component.memSeason.Data.(lblCurr).p = [];
               component.memSeason.Data.(lblCurr).val = [];
               component.memSeason.Data.(lblCurr).stat = [];
            end
                        
            % Exo Memory 
            % #################################################
            % ############ Need Update ########################
            % #################################################
            component.memExo = struct;   
            
            % Exo Signal
            % #################################################
            % ############ Need Update ########################
            % ################################################# 
            component.exoSignal = struct; 
            
            obj.component = component;
        end
        
        % Set whether model should display messages or not
        function setDisplay(obj, display)
            obj.display = display;
        end
        
        % Initializing Data
        function initData(obj)
            %   Initiates the data and make it ready for further
            %   computation: 
            %   - Computing trial_alignment based on signals
            %   - Computing the firing rates
            %   - Computing ITI
            
            obj.displayMsg('Start initiation process...');
            
            % Computing ITI
            obj.compITI();
            
            % Compute absolute time intevals
            startPointEst = obj.neuronData.maxTrialLen-obj.neuronData.endOfTrialDist;
            alignSignTime = obj.neuronData.signalTime.(obj.neuronData.endOfTrialSignal);
            
            obj.timeIntervarlsAbsolute = repmat([-startPointEst:obj.neuronData.binSize:obj.neuronData.endOfTrialDist],[size(alignSignTime,1),1]) + ...
                repmat(alignSignTime,[1, size([-startPointEst:obj.neuronData.binSize:obj.neuronData.endOfTrialDist],2)]);
            
            % Calculating Fairing Rate Matrix
            obj.firingRateMat = [];
            
            trialFilter = obj.startPoint:size(obj.timeIntervarlsAbsolute,1)-1; % removing one trial on begening and one at the end for the sake of binning the spike trian
            
            for cntTrial = trialFilter                
                % All spike in current trial
                allSpikes = obj.neuronData.spikeTime;
                
                trialTimeIntervarlsAbsolute = obj.timeIntervarlsAbsolute(cntTrial,:);
                
                spikeCount = histc(allSpikes,trialTimeIntervarlsAbsolute);
                spikeCount = spikeCount(1:end-1);
                
                if isempty(spikeCount)
                    spikeCount = zeros(length(trialTimeIntervarlsAbsolute)-1,1);
                end
                
                if size(spikeCount,1)>1
                    spikeCount = spikeCount';
                end
                
                % make the spike count equal to nan, if the bin is in
                % previous trial
                spikeCount((trialTimeIntervarlsAbsolute+obj.neuronData.binSize) < alignSignTime(cntTrial-1)+obj.neuronData.endOfTrialDist) = nan;
                
                % storing firing rate
                obj.firingRateMat = [obj.firingRateMat; spikeCount/(obj.neuronData.binSize/1000)];                
            end                     
            
            % triming absolute times based on "trialFilter"
            obj.timeIntervarlsAbsolute = obj.timeIntervarlsAbsolute(trialFilter,2:end);
            
            % triming signal value based on "trialFilter"
            signalValueNames = fieldnames(obj.neuronData.signalValue);
            for cntSignal = 1:length(signalValueNames)
                obj.neuronData.signalValue.(signalValueNames{cntSignal}) = obj.neuronData.signalValue.(signalValueNames{cntSignal})(trialFilter);
            end
            
            % triming signal times based on "trialFilter"
            signalTimeNames = fieldnames(obj.neuronData.signalTime);
            for cntSignal = 1:length(signalTimeNames)
                obj.neuronData.signalTime.(signalTimeNames{cntSignal}) = obj.neuronData.signalTime.(signalTimeNames{cntSignal})(trialFilter);
            end
            
            % calculate mean firing pattern in profile
            obj.firingRateMatMean = nanmean(obj.firingRateMat,1);                        

            % calculating boudries
            obj.calcTrialBinBoudries();
            
            % calculating model components
            obj.initMemShort();
            obj.initMemSeason();
            obj.initMemExo();
            obj.initExoSignal();
            
%             % Setup parameter set and boundries
%             obj.calcParamSet();

            obj.displayMsg('Initiation process is done!');
        end
               
        % Setup fitting procedure information
        function setFittingSet(obj, iterNum, crossPer, crossIter, modelSelection)            
            obj.displayMsg('Setup fitting parameters...');
            
            if nargin <2
                iterNum = 20;
                crossPer = 0.2;
                crossIter = 10;               
                modelSelection = 1;
            elseif nargin <3                
                crossPer = 0.2;
                crossIter = 10;    
                modelSelection = 1;
            elseif nargin <4                                
                crossIter = 10;    
                modelSelection = 1;
            end
            obj.fittingSet.iterNum = iterNum;
            obj.fittingSet.crossPer = crossPer;
            obj.fittingSet.crossIter = crossIter;            
            obj.fittingSet.modelSelection = modelSelection;           
            
            obj.fittingSet.selectedComponents.seqOrderName = {}; % a dictonary for converting seq order num to name of components
            obj.fittingSet.selectedComponents.totalPermutationNumber = []; % total number of model permutations
            obj.fittingSet.selectedComponents.currentSelectedSeq = []; % The current selected sequence of components (binary seq)
            
            % Selection of ExoSignals
            obj.fittingSet.selectedComponents.exoSignal.selected = ~isempty(fieldnames(obj.component.exoSignal));
            obj.fittingSet.selectedComponents.seqOrderName{end+1} = 'exoSignal';
            obj.fittingSet.selectedComponents.exoSignal.seqNumber = length(obj.fittingSet.selectedComponents.seqOrderName);
            
            % Selection of short and seasonal memories
            obj.fittingSet.selectedComponents.memShort.selected = obj.component.memShort.defined;
            obj.fittingSet.selectedComponents.seqOrderName{end+1} = 'memShort';
            obj.fittingSet.selectedComponents.memShort.seqNumber = length(obj.fittingSet.selectedComponents.seqOrderName);
            
            obj.fittingSet.selectedComponents.memSeason.selected = obj.component.memSeason.defined;
            obj.fittingSet.selectedComponents.seqOrderName{end+1} = 'memSeason';
            obj.fittingSet.selectedComponents.memSeason.seqNumber = length(obj.fittingSet.selectedComponents.seqOrderName);
            

            % Calculating exogenous memory selection
            if ~isempty(fieldnames(obj.component.memExo))
                exoMemNames = fieldnames(obj.fittingModelData.memExo);
                
                for cntExoMem = 1:length(exoMemNames)
                    exoMemNameTemp = exoMemNames{cntExoMem};
                    obj.fittingSet.selectedComponents.memExo.(exoMemNameTemp).selected = 1;
                    obj.fittingSet.selectedComponents.seqOrderName{end+1} = ['memExo.',exoMemNameTemp];
                    obj.fittingSet.selectedComponents.memExo.(exoMemNameTemp).seqNumber = length(obj.fittingSet.selectedComponents.seqOrderName);
                end
            end
                        
            obj.fittingSet.selectedComponents.totalPermutationNumber = 2.^length(obj.fittingSet.selectedComponents.seqOrderName);
            obj.fittingSet.selectedComponents.currentSelectedSeq = ones(1,length(obj.fittingSet.selectedComponents.seqOrderName));
        end
        
        % Fit the model
        function fit(obj)        
           obj.saveModel();
            
            if obj.fittingSet.modelSelection
                obj.displayMsg('Start Model Selection Process...');
                
                for cntModel = 0:obj.fittingSet.selectedComponents.totalPermutationNumber-1
                    obj.displayMsg(['Running model [',num2str(cntModel),']...']);
                                        
                    % Set current selection Mask
                    obj.fittingSet.selectedComponents.currentSelectedSeq = dec2bin(cntModel,length(obj.fittingSet.selectedComponents.currentSelectedSeq));
                    obj.fittingSet.selectedComponents.currentSelectedSeq = obj.fittingSet.selectedComponents.currentSelectedSeq(end:-1:1);
                    
                    obj.fittingSet.currentModel = cntModel+1;                    
                    
                    % Setup parameter set and boundries
                    obj.calcParamSet();
            
                    % Fit a single model
                    fit_singleModel(obj);
                    
                    obj.displayMsg(['End of model [',num2str(cntModel),']...']);
                end
                
                obj.chooseBestModel();                                
                
            else
                obj.displayMsg('Running General Model...');
                cntModel = 0;
                obj.fittingSet.currentModel = cntModel+1;
                
                % Set current selection Mask
                obj.fittingSet.selectedComponents.currentSelectedSeq = dec2bin(obj.fittingSet.selectedComponents.totalPermutationNumber-1,length(obj.fittingSet.selectedComponents.currentSelectedSeq));
                obj.fittingSet.selectedComponents.currentSelectedSeq = obj.fittingSet.selectedComponents.currentSelectedSeq(end:-1:1);
                
                % Setup parameter set and boundries
                obj.calcParamSet();
                
                % Fit a single model
                fit_singleModel(obj);
                obj.displayMsg('End of General Model...');
            end
             
            obj.saveModelSnapshot(1);
            obj.saveModel(1);
        end
        
        % Fit a single model
        function fit_singleModel(obj)           
            dataMask = [];
            dataMask.mask = logical(ones(size(obj.fittingModelData.output)))';
            minFunc = @(v,d) obj.costfunction(v,d.mask);
            positiveFlag = obj.paramSet.limit.sorted(:,1) >=0;
            limitSet = {};
            for cntV = 1:size(obj.paramSet.limit.sorted,1)
                limitSet{cntV} = obj.paramSet.limit.sorted(cntV,:);
            end
            
            [vars, stat, ret] = Fitting_Function(minFunc, dataMask, size(obj.paramSet.limit.sorted,1),...
                obj.fittingSet.iterNum, obj.fittingSet.crossPer, obj.fittingSet.crossIter,...
                positiveFlag, limitSet, {'mask'}, 1);           
            
            obj.fittingResults.models{obj.fittingSet.currentModel}.res.vars = vars;
            obj.fittingResults.models{obj.fittingSet.currentModel}.res.stat = stat;
            obj.fittingResults.models{obj.fittingSet.currentModel}.res.ret = ret;
            
            obj.fittingResults.models{obj.fittingSet.currentModel}.paramSet = obj.paramSet;
            
            obj.saveModelSnapshot();
        end
    end
    
    methods(Access = private)
        % Display a message
        function displayMsg(obj, msg)
            if obj.display
                display(msg);
            end
        end
        
        % Compute InterTrial Intervals
        function compITI(obj)
            %   Computing ITI for all trials
            obj.displayMsg('Computing ITI for all trials...');
            
            for cntTemp= 2:obj.neuronData.trialRange(end)
                ITITemp(cntTemp) = obj.neuronData.signalTime.(obj.neuronData.beginOfTrialSignal)(cntTemp) - obj.neuronData.signalTime.(obj.neuronData.endOfTrialSignal)(cntTemp-1); % All the ITIs
            end
            obj.ITI = ITITemp;
        end                
        
        % Compute Short term memory
        function initMemShort(obj)
            if (obj.component.memShort.defined)
                %   Setup the neural data for Short memory
                obj.displayMsg('Setup the neural data for Short memory...');
                
                dataPointIdx = 0;
                obj.fittingModelData.memShort.Rate = zeros(obj.totalDataPoint, obj.component.memShort.order);
                obj.fittingModelData.memShort.Time = zeros(obj.totalDataPoint, obj.component.memShort.order);
                
                for cntTrials = obj.boundries.trial
                    for cntBins = obj.boundries.bin
                        % clc
                        if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                            
                            % calculating the stating bin number of non nan
                            % data
                            startPointNoNAN = sum(isnan(obj.firingRateMat(cntTrials,:)))+1;
                            
                            % calculating the sliding windeo; negative values
                            % indicate the value from previouse trial
                            binsSlidingBox = [cntBins-obj.component.memShort.order:cntBins-1];
                            binsSlidingBox = [binsSlidingBox(binsSlidingBox<startPointNoNAN)-startPointNoNAN+1 binsSlidingBox(binsSlidingBox>=startPointNoNAN)];
                            
                            %
                            %                         yTemp = spikeCountsGeneral(cntTrials,cntBins);
                            %                         yMeanTemp = spikeCountsGeneralMean(cntBins);
                            %
                            %                         yTempZS = spikeCountsGeneralZS(cntTrials,cntBins);
                            %                         yMeanTempZS = spikeCountsGeneralMeanZS(cntBins);
                            
                            maxBinLen = size(obj.firingRateMat,2);
                            
                            % Initializing data storages
                            ARSTemp = zeros(1,obj.component.memShort.order);
                            timeS = zeros(1,obj.component.memShort.order);
                            
                            cntId = 1;
                            for cntSlidingBox = binsSlidingBox
                                if cntSlidingBox>0
                                    ARSTemp(cntId) = obj.firingRateMat(cntTrials,cntSlidingBox);
                                    timeS(cntId) = obj.timeIntervarlsAbsolute(cntTrials,cntSlidingBox);
                                    cntId = cntId + 1;
                                else
                                    ARSTemp(cntId) = obj.firingRateMat(cntTrials-1,maxBinLen+cntSlidingBox);
                                    timeS(cntId) = obj.timeIntervarlsAbsolute(cntTrials-1,maxBinLen+cntSlidingBox);
                                    cntId = cntId + 1;
                                end
                            end
                            
                            for cntBinTimeS = 1:length(timeS)-1
                                timeS(cntBinTimeS) = timeS(cntBinTimeS+1)-timeS(cntBinTimeS);
                            end
                            timeS(end) = obj.neuronData.binSize;
                            
                            % Inverse the vectors
                            ARSTemp = ARSTemp(end:-1:1);
                            timeS = timeS(end:-1:1);
                            
                            dataPointIdx = dataPointIdx + 1;
                            obj.fittingModelData.memShort.Rate(dataPointIdx,:) = ARSTemp;
                            obj.fittingModelData.memShort.Time(dataPointIdx,:) = timeS;
                            
                        end
                        
                    end
                end
                
                obj.fittingModelData.memShort.mask = ones(obj.totalDataPoint, obj.component.memShort.order);
                if ~isempty(obj.component.memShort.effectiveTimeWin)
                    obj.fittingModelData.memShort.mask = obj.clcTimWin(obj.component.memShort.effectiveTimeWin, obj.fittingModelData.memShort.mask);
                end
            end
        end       
        
        % Compute Seasonal memory
        function initMemSeason(obj)
            if (obj.component.memSeason.defined)
                %   Setup the neural data for Seasonal memory
                obj.displayMsg('Setup the neural data for Seasonal memory...');
                
                dataPointIdx = 0;
                obj.fittingModelData.memSeason.Rate = zeros(obj.totalDataPoint, obj.component.memSeason.order);
                obj.fittingModelData.memSeason.Time = zeros(obj.totalDataPoint, obj.component.memSeason.order);
                for cntTrials = obj.boundries.trial
                    for cntBins = obj.boundries.bin
                        % clc
                        if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                            
                            ARLTemp = obj.firingRateMat(cntTrials-obj.component.memSeason.order:cntTrials-1,cntBins);
                            timeL = obj.timeIntervarlsAbsolute(cntTrials-obj.component.memSeason.order:cntTrials-1,cntBins);
                            
                            % Inverse the vectors
                            ARLTemp = ARLTemp(end:-1:1);
                            timeL = timeL(end:-1:1);
                            
                            dataPointIdx = dataPointIdx + 1;
                            obj.fittingModelData.memSeason.Rate(dataPointIdx,:) = ARLTemp;
                            obj.fittingModelData.memSeason.Time(dataPointIdx,:) = timeL;
                            
                        end
                        
                    end
                end
                
                obj.fittingModelData.memSeason.mask = ones(obj.totalDataPoint, obj.component.memSeason.order);
                if ~isempty(obj.component.memSeason.effectiveTimeWin)
                    obj.fittingModelData.memSeason.mask = obj.clcTimWin(obj.component.memSeason.effectiveTimeWin, obj.fittingModelData.memSeason.mask);
                end
            end
        end
        
        % Compute exogenous memories
        function initMemExo(obj)
            if ~isempty(fieldnames(obj.component.memExo))
                %   Setup the neural data for exo memories
                obj.displayMsg('Setup the neural data for exo memories...');
                
                exoMemNames = fieldnames(obj.component.memExo);
                
                for cntExoMem = 1:length(exoMemNames)
                    exoMemNameTemp = exoMemNames{cntExoMem};
                    exoMemTemp = obj.component.memExo.(exoMemNameTemp);
                    
                    if ~exoMemTemp.interaction
                        if ~exoMemTemp.exoSelEff
                            
                            dataPointIdx = 0;
                            obj.fittingModelData.memExo.(exoMemNameTemp).Signal = zeros(obj.totalDataPoint, exoMemTemp.memOrder);
                            obj.fittingModelData.memExo.(exoMemNameTemp).Time = zeros(obj.totalDataPoint, exoMemTemp.memOrder);
                            for cntTrials = obj.boundries.trial
                                for cntBins = obj.boundries.bin
                                    % clc
                                    if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                                        
                                        % clc
                                        % current time bigger than the signal
                                        if obj.timeIntervarlsAbsolute(cntTrials,cntBins) > (obj.neuronData.signalTime.(exoMemTemp.signalTime)(cntTrials))
                                            exoSigTemp = obj.neuronData.signalValue.(exoMemTemp.signalName)(cntTrials-exoMemTemp.memOrder+1:cntTrials);
                                            exoTimeTemp = obj.timeIntervarlsAbsolute(cntTrials,cntBins) -...
                                                obj.neuronData.signalTime.(exoMemTemp.signalTime)(cntTrials-exoMemTemp.memOrder+1:cntTrials);
                                            
                                        else % current time smaller than the signal
                                            exoSigTemp = obj.neuronData.signalValue.(exoMemTemp.signalName)(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                            exoTimeTemp = obj.timeIntervarlsAbsolute(cntTrials,cntBins) -...
                                                obj.neuronData.signalTime.(exoMemTemp.signalTime)(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                        end
                                        
                                        % Inverse the vectors
                                        exoSigTemp = exoSigTemp(end:-1:1);
                                        exoTimeTemp = exoTimeTemp(end:-1:1);
                                        
                                        dataPointIdx = dataPointIdx + 1;
                                        obj.fittingModelData.memExo.(exoMemNameTemp).Signal(dataPointIdx,:) = exoSigTemp;
                                        obj.fittingModelData.memExo.(exoMemNameTemp).Time(dataPointIdx,:) = exoTimeTemp;
                                        
                                    end
                                end
                            end
                            
                        else
                            
                            dataPointIdx = 0;
                            obj.fittingModelData.memExo.(exoMemNameTemp).Signal = zeros(obj.totalDataPoint, exoMemTemp.memOrder);
                            obj.fittingModelData.memExo.(exoMemNameTemp).Time = zeros(obj.totalDataPoint, exoMemTemp.memOrder);
                            for cntTrials = obj.boundries.trial
                                for cntBins = obj.boundries.bin
                                    % clc
                                    if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                                        
                                        exoSigTemp = obj.neuronData.signalValue.(exoMemTemp.signalName)(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                        exoTimeTemp = obj.timeIntervarlsAbsolute(cntTrials,cntBins) -...
                                            obj.neuronData.signalTime.(exoMemTemp.signalTime)(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                        
                                        % Inverse the vectors
                                        exoSigTemp = exoSigTemp(end:-1:1);
                                        exoTimeTemp = exoTimeTemp(end:-1:1);
                                        
                                        dataPointIdx = dataPointIdx + 1;
                                        obj.fittingModelData.memExo.(exoMemNameTemp).Signal(dataPointIdx,:) = exoSigTemp;
                                        obj.fittingModelData.memExo.(exoMemNameTemp).Time(dataPointIdx,:) = exoTimeTemp;
                                        
                                    end
                                end
                            end
                            
                        end
                    else
                        if ~exoMemTemp.exoSelEff
                            
                            dataPointIdx = 0;
                            obj.fittingModelData.memExo.(exoMemNameTemp).Signal = zeros(obj.totalDataPoint, exoMemTemp.memOrder);
                            obj.fittingModelData.memExo.(exoMemNameTemp).Time = zeros(obj.totalDataPoint, exoMemTemp.memOrder);
                            for cntTrials = obj.boundries.trial
                                for cntBins = obj.boundries.bin
                                    % clc
                                    if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                                        
                                        % clc
                                        % current time bigger than the signal
                                        if obj.timeIntervarlsAbsolute(cntTrials,cntBins) > (obj.neuronData.signalTime.(exoMemTemp.signalTime)(cntTrials))
                                            for cntSignals = 1:length(exoMemTemp.signalName)
                                                exoSigTemp(:,cntSignals) = obj.neuronData.signalValue.(exoMemTemp.signalName{cntSignals})(cntTrials-exoMemTemp.memOrder+1:cntTrials);
                                            end
                                            
                                            exoSigTemp = obj.interactionCalc(exoSigTemp, exoMemTemp.interactionType);
                                            exoTimeTemp = obj.timeIntervarlsAbsolute(cntTrials,cntBins) -...
                                                obj.neuronData.signalTime.(exoMemTemp.signalTime)(cntTrials-exoMemTemp.memOrder+1:cntTrials);
                                            
                                        else % current time smaller than the signal
                                            for cntSignals = 1:length(exoMemTemp.signalName)
                                                exoSigTemp(:,cntSignals) = obj.neuronData.signalValue.(exoMemTemp.signalName{cntSignals})(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                            end
                                            
                                            exoSigTemp = obj.interactionCalc(exoSigTemp, exoMemTemp.interactionType);
                                            exoTimeTemp = obj.timeIntervarlsAbsolute(cntTrials,cntBins) -...
                                                obj.neuronData.signalTime.(exoMemTemp.signalTime)(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                        end
                                        
                                        % Inverse the vectors
                                        exoSigTemp = exoSigTemp(end:-1:1);
                                        exoTimeTemp = exoTimeTemp(end:-1:1);
                                        
                                        dataPointIdx = dataPointIdx + 1;
                                        obj.fittingModelData.memExo.(exoMemNameTemp).Signal(dataPointIdx,:) = exoSigTemp;
                                        obj.fittingModelData.memExo.(exoMemNameTemp).Time(dataPointIdx,:) = exoTimeTemp;
                                        
                                    end
                                end
                            end
                            
                        else
                            
                            dataPointIdx = 0;
                            obj.fittingModelData.memExo.(exoMemNameTemp).Signal = zeros(obj.totalDataPoint, exoMemTemp.memOrder);
                            obj.fittingModelData.memExo.(exoMemNameTemp).Time = zeros(obj.totalDataPoint, exoMemTemp.memOrder);
                            for cntTrials = obj.boundries.trial
                                for cntBins = obj.boundries.bin
                                    % clc
                                    if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                                        
                                        for cntSignals = 1:length(exoMemTemp.signalName)
                                            exoSigTemp(:,cntSignals) = obj.neuronData.signalValue.(exoMemTemp.signalName{cntSignals})(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                        end
                                        
                                        exoSigTemp = obj.interactionCalc(exoSigTemp, exoMemTemp.interactionType);
                                        exoTimeTemp = obj.timeIntervarlsAbsolute(cntTrials,cntBins) -...
                                            obj.neuronData.signalTime.(exoMemTemp.signalTime)(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                        
                                        % Inverse the vectors
                                        exoSigTemp = exoSigTemp(end:-1:1);
                                        exoTimeTemp = exoTimeTemp(end:-1:1);
                                        
                                        dataPointIdx = dataPointIdx + 1;
                                        obj.fittingModelData.memExo.(exoMemNameTemp).Signal(dataPointIdx,:) = exoSigTemp;
                                        obj.fittingModelData.memExo.(exoMemNameTemp).Time(dataPointIdx,:) = exoTimeTemp;
                                        
                                    end
                                end
                            end
                            
                        end
                    end
                    
                    obj.fittingModelData.memExo.(exoMemNameTemp).meanData = obj.firingRateMatMeanArranged;
                    
                    obj.fittingModelData.memExo.(exoMemNameTemp).mask = ones(obj.totalDataPoint, obj.component.memExo.(exoMemNameTemp).memOrder);
                    if ~isempty(obj.component.memExo.(exoMemNameTemp).effectiveTimeWin)
                        obj.fittingModelData.memExo.(exoMemNameTemp).mask = obj.clcTimWin(obj.component.memExo.(exoMemNameTemp).effectiveTimeWin, obj.fittingModelData.memExo.(exoMemNameTemp).mask);
                    end
                    
                end
            end
        end
        
        % Compute exogenous signals
        function initExoSignal(obj)
            if ~isempty(fieldnames(obj.component.exoSignal))
                %   Setup the neural data for exo signals
                obj.displayMsg('Setup the neural data for exo signals...');
                
                exoSignalNames = fieldnames(obj.component.exoSignal);
                
                for cntExoSignal = 1:length(exoSignalNames)
                    exoSignalNameTemp = exoSignalNames{cntExoSignal};
                    exoSignalTemp = obj.component.exoSignal.(exoSignalNameTemp);
                    
                    if ~exoSignalTemp.interaction
                        dataPointIdx = 0;
                        obj.fittingModelData.exoSignal.(exoSignalNameTemp).Signal = zeros(obj.totalDataPoint, 1);
                        
                        for cntTrials = obj.boundries.trial
                            for cntBins = obj.boundries.bin
                                % clc
                                if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                                    
                                    % clc
                                    % current time bigger than the signal
                                    if obj.timeIntervarlsAbsolute(cntTrials,cntBins) > (obj.neuronData.signalTime.(exoSignalTemp.signalName)(cntTrials))
                                        exoSigTemp = obj.neuronData.signalValue.(exoSignalTemp.signalName)(cntTrials-1+1:cntTrials);
                                    else % current time smaller than the signal
                                        exoSigTemp = obj.neuronData.signalValue.(exoSignalTemp.signalName)(cntTrials-1:cntTrials-1);
                                    end
                                    
                                    % Inverse the vectors
                                    exoSigTemp = exoSigTemp(end:-1:1);
                                    dataPointIdx = dataPointIdx + 1;
                                    obj.fittingModelData.exoSignal.(exoSignalNameTemp).Signal(dataPointIdx,:) = exoSigTemp;
                                end
                            end
                        end
                        
                    else
                        
                        dataPointIdx = 0;
                        obj.fittingModelData.exoSignal.(exoSignalNameTemp).Signal = zeros(obj.totalDataPoint, 1);
                        
                        for cntTrials = obj.boundries.trial
                            for cntBins = obj.boundries.bin
                                % clc
                                if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                                    
                                    % clc
                                    for cntSignals = 1:length(exoSignalTemp.signalName)
                                        % current time bigger than the signal
                                        if obj.timeIntervarlsAbsolute(cntTrials,cntBins) > (obj.neuronData.signalTime.(exoSignalTemp.signalName{cntSignals})(cntTrials))
                                            exoSigTemp(:,cntSignals) = obj.neuronData.signalValue.(exoSignalTemp.signalName{cntSignals})(cntTrials-1+1:cntTrials);
                                        else % current time smaller than the signal
                                            exoSigTemp(:,cntSignals) = obj.neuronData.signalValue.(exoSignalTemp.signalName{cntSignals})(cntTrials-1:cntTrials-1);
                                        end
                                    end
                                    
                                    exoSigTemp = obj.interactionCalc(exoSigTemp, exoSignalTemp.interactionType);
                                    
                                    % Inverse the vectors
                                    exoSigTemp = exoSigTemp(end:-1:1);
                                    dataPointIdx = dataPointIdx + 1;
                                    obj.fittingModelData.exoSignal.(exoSignalNameTemp).Signal(dataPointIdx,:) = exoSigTemp;
                                end
                            end
                        end
                        
                    end
                    
                    obj.fittingModelData.exoSignal.(exoSignalNameTemp).mask = ones(obj.totalDataPoint, 1);
                    if ~isempty(obj.component.exoSignal.(exoSignalNameTemp).effectiveTimeWin)
                        obj.fittingModelData.exoSignal.(exoSignalNameTemp).mask = obj.clcTimWin(obj.component.exoSignal.(exoSignalNameTemp).effectiveTimeWin, obj.fittingModelData.exoSignal.(exoSignalNameTemp).mask);
                    end
                end
            end
        end       
        
        % Compute trial boundries
        function calcTrialBinBoudries(obj)
            % calculate maximum order
            allOrders = obj.component.memSeason.order;
            exoMemNames = fieldnames(obj.component.memExo);
            for cntExoMems = 1:length(exoMemNames)
                allOrders(end+1) = obj.component.memExo.(exoMemNames{cntExoMems}).memOrder;
            end
            
            obj.maxOrder = max(allOrders);
            
            % calculating boudries
            obj.boundries.trial = obj.maxOrder+2:size(obj.firingRateMat,1);
            obj.boundries.bin = 1:size(obj.firingRateMat,2);
            
            % Calculating total datapoints
            obj.totalDataPoint = 0;
            obj.fittingModelData.output = [];
            obj.firingRateMatMeanArranged = [];
            for cntTrials = obj.boundries.trial
                for cntBins = obj.boundries.bin
                    
                    if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                        obj.totalDataPoint = obj.totalDataPoint + 1;
                        obj.fittingModelData.output(end+1) = obj.firingRateMat(cntTrials,cntBins);
                        obj.firingRateMatMeanArranged(end+1) = obj.firingRateMatMean(cntBins);
                    end
                end
            end
        end
        
        % Compute time windows for time mask
        function tempMask = clcTimWin(obj, timeWin, tempMask)
            % initiate a temp mask
            dataPointIdx = 0;
%             tempMask = ones(obj.totalDataPoint, obj.component.memShort.order);
            
            if timeWin.type == 1
                for cntTrials = obj.boundries.trial
                    for cntBins = obj.boundries.bin
                        if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                            dataPointIdx = dataPointIdx + 1;
                            
                            if ~ismember(cntBins,timeWin.binWindow) 
                                tempMask(dataPointIdx,:) = tempMask(dataPointIdx,:) * 0;
                            end
                        end
                    end
                end
            else
                
                for cntTrials = obj.boundries.trial
                    for cntBins = obj.boundries.bin
                        % clc
                        if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                            dataPointIdx = dataPointIdx + 1;
                            % clc
                            % current time smaller than the mask begin time
                            if obj.timeIntervarlsAbsolute(cntTrials,cntBins) < (obj.neuronData.signalTime.(timeWin.beginSig_Name)(cntTrials)+timeWin.beginSig_TimeDist)
                                tempMask(dataPointIdx,:) = tempMask(dataPointIdx,:) * 0;
                            % current time biger than the mask end time
                            elseif obj.timeIntervarlsAbsolute(cntTrials,cntBins) > (obj.neuronData.signalTime.(timeWin.endSig_Name)(cntTrials)+timeWin.endSig_TimeDist) 
                                tempMask(dataPointIdx,:) = tempMask(dataPointIdx,:) * 0;
                            end
                            
                        end
                    end
                end
                
            end
            
        end
        
        % Compute interaction data
        function data = interactionCalc(obj, data, interactionType)
            % Assumption is that size(data,2) shows the number of
            % parameters involved in interaction and values can be either 1
            % or -1.
            %
            % Function retuns a unique value correspond to the binary value of
            % input between -(paramNum^2)/2 and (paramNum^2)/2, excluding 0
           
            if interactionType == 1
                data = prod(data,2);      
            elseif interactionType == 2
                paramNum = size(data,2);
                dataPos = (data+1)/2;
                tempOutput = bin2dec(num2str(dataPos))+1;
                if tempOutput <= (paramNum^2)/2
                    tempOutput = tempOutput - ((paramNum^2)/2 + 1);
                else
                    tempOutput = tempOutput - ((paramNum^2)/2);
                end
            end
            data = tempOutput;
        end
        
        % Setup parameter set and boundries
        function calcParamSet(obj)
            
            obj.displayMsg('Setup parameter settings...');
            
            curly = @(x, varargin) x{varargin{:}};
            
            obj.paramSet.size.bias          = 1;
            obj.paramSet.size.exoSignal     = 1;
            obj.paramSet.size.memShort      = obj.component.memShort.order;
            obj.paramSet.size.memSeason     = obj.component.memSeason.order;
            obj.paramSet.size.memExo        = 2;
            obj.paramSet.size.mean          = 1;
            
            
            obj.paramSet.flag.bias          = [];
            obj.paramSet.flag.exoSignal     = [];
            obj.paramSet.flag.memShort      = [];
            obj.paramSet.flag.memSeason     = [];
            obj.paramSet.flag.memExo        = [];
            obj.paramSet.flag.mean          = [];
            
            obj.paramSet.effect.bias        = 0;
            obj.paramSet.effect.exoSignal   = 0;
            obj.paramSet.effect.memShort    = 0;
            obj.paramSet.effect.memSeason   = 0;
            obj.paramSet.effect.memExo      = 0;
            obj.paramSet.effect.mean        = 0;
                       
            
            obj.paramSet.limit.bias         = [-inf inf];
            obj.paramSet.limit.exoSignal    = [-inf inf];
            obj.paramSet.limit.memShort     = [-inf inf];
            obj.paramSet.limit.memSeason    = [-inf inf];
            obj.paramSet.limit.memExo       = [-inf inf];            
            obj.paramSet.limit.mean         = [-inf inf];
            
            
            obj.paramSet.flag.bias = 1;
            currFlag = 1;
            obj.paramSet.effect.bias = 1;
            
            obj.paramSet.limit.sorted       = [];
            obj.paramSet.limit.sorted(1,:)  = [-inf inf];

            
            if ~isempty(fieldnames(obj.component.exoSignal)) && ...
                    str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.exoSignal.seqNumber))
                obj.paramSet.flag.exoSignal = [currFlag+1:currFlag+obj.paramSet.size.exoSignal*(length(fieldnames(obj.component.exoSignal)))];
                currFlag = obj.paramSet.flag.exoSignal(end);
                obj.paramSet.effect.exoSignal = 1;
                
                for cntN = 1:length(fieldnames(obj.component.exoSignal))                    
                    obj.paramSet.limit.exoSignal((cntN-1)*obj.paramSet.size.exoSignal+1:cntN*obj.paramSet.size.exoSignal,:) = ...
                        obj.component.exoSignal.(curly(fieldnames(obj.component.exoSignal),cntN)).paramLim;
                    obj.paramSet.limit.sorted(end+1:end+obj.paramSet.size.exoSignal,:)  = ...
                        obj.component.exoSignal.(curly(fieldnames(obj.component.exoSignal),cntN)).paramLim;
                end
            end                        
            
            if (obj.component.memShort.defined) && ...
                    str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.memShort.seqNumber))
                obj.paramSet.flag.memShort = [currFlag+1:currFlag+obj.paramSet.size.memShort];
                currFlag = obj.paramSet.flag.memShort(end);
                obj.paramSet.effect.memShort = 1;
                obj.paramSet.limit.memShort(1:obj.paramSet.size.memShort,:) = repmat(obj.component.memShort.paramLim,[obj.paramSet.size.memShort,1]);
                obj.paramSet.limit.sorted(end+1:end+obj.paramSet.size.memShort,:)  = repmat(obj.component.memShort.paramLim,[obj.paramSet.size.memShort,1]);
            end
            
            if (obj.component.memSeason.defined) && ...
                    str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.memSeason.seqNumber))
                obj.paramSet.flag.memSeason = [currFlag+1:currFlag+obj.paramSet.size.memSeason];
                currFlag = obj.paramSet.flag.memSeason(end);
                obj.paramSet.effect.memSeason = 1;
                obj.paramSet.limit.memSeason(1:obj.paramSet.size.memSeason,:) = repmat(obj.component.memSeason.paramLim,[obj.paramSet.size.memSeason,1]);
                obj.paramSet.limit.sorted(end+1:end+obj.paramSet.size.memSeason,:)  = repmat(obj.component.memSeason.paramLim,[obj.paramSet.size.memSeason,1]);
            end
            
            if ~isempty(fieldnames(obj.component.memExo))
                exoMemNames = fieldnames(obj.fittingModelData.memExo);
                exoEffeciveNumTemp = 0; % store the numbers of effective exo memory components
                for cntN = 1:length(exoMemNames)
                    exoMemNameTemp = exoMemNames{cntN};
                    if str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.memExo.(exoMemNameTemp).seqNumber))
                        obj.paramSet.limit.memExo((cntN-1)*obj.paramSet.size.memExo+1:cntN*obj.paramSet.size.memExo,:) = ...
                            obj.component.memExo.(curly(fieldnames(obj.component.memExo),cntN)).paramLim;
                        obj.paramSet.limit.sorted(end+1:end+obj.paramSet.size.memExo,:)  = ...
                            obj.component.memExo.(curly(fieldnames(obj.component.memExo),cntN)).paramLim;
                        
                        exoEffeciveNumTemp = exoEffeciveNumTemp + 1;
                    end
                end
                
                if exoEffeciveNumTemp>0
                    obj.paramSet.flag.memExo = [currFlag+1:currFlag+obj.paramSet.size.memExo*exoEffeciveNumTemp];
                    currFlag = obj.paramSet.flag.memExo(end);
                    obj.paramSet.effect.memExo = 1;
                end
                
                
                % Add mean flag
                obj.paramSet.flag.mean = [currFlag+1];
                currFlag = obj.paramSet.flag.mean;
                obj.paramSet.effect.mean = 1;
                obj.paramSet.limit.sorted(end+1,:) = obj.paramSet.limit.mean;
            end

        end                
        
        % Compute output of the model
        function ret = modelOutput(obj, params, dataMask)    
            if nargin < 3
                dataMask = logical(ones(size(obj.fittingModelData.output)))';
            end
            
            ex              = @(x ,t) x(1).*(exp(-t./x(2)));
            
            anonymTools;
            
            % Validate the parameters
            rejectParams = 0;
            
            interupt = 0;
            
            if size(params,2) ~= size(obj.paramSet.limit.sorted,1)
                interupt = 1;
            end
                        
            for cntP = 1:size(obj.paramSet.limit.sorted,1)
                if (params(cntP) < obj.paramSet.limit.sorted(cntP,1)) | (params(cntP) > obj.paramSet.limit.sorted(cntP,2))
                    interupt = 1;
                end
            end                       
            
            % Run the model
            if interupt
                totOutputTemp = ones(size(obj.fittingModelData.output(dataMask'))) * 10^45;
            else
%                 ex              = @(x ,t) x(1).*(exp(-t./x(2))+x(3));
%                 ex_rew          = @(x ,t) x(1).*(exp(-t./x(2)));% + (abs(x(1))>4 | x(2)>(4*20)) * 10^45;
%                 ex_choice       = @(x ,t) x(1).*(exp(-t./x(2)));
%                 %     ex_rew          = @(x ,t) x(1).*(exp(-t./x(2)));% + (abs(x(1))>4 | x(2)>(4*20)) * 10^45;
%                 ex_l            = @(x ,t) x(1).*(exp(-t./x(2))+x(3));% + (abs(x(1))>4 | x(2)>(4*20)) * 10^45;
                
                totInput        =   min(obj.totalDataPoint, sum(dataMask));
                
                biasPart        =   0;
                exoSignalPart   =   zeros(totInput,1);
                memShortPart    =   zeros(totInput,1);
                memSeasonPart   =   zeros(totInput,1);
                memExoPart      =   zeros(totInput,1);
                ARMeanPart      =   zeros(totInput,1);
                
                % zscoring if needed
                if obj.zscoring == 1
                    obj.zscoreData(dataMask);
                end
                
                % Calculating bias component
                if obj.paramSet.effect.bias
                    biasPart    =   params(obj.paramSet.flag.bias);
                end
                
                % Calculating exogenous signal(s) component
                if obj.paramSet.effect.exoSignal
                    exoSignalNames = fieldnames(obj.fittingModelData.exoSignal);
                    
                    for cntExoSignal = 1:length(exoSignalNames)
                        exoSignalNameTemp = exoSignalNames{cntExoSignal};
%                         exoSignalTemp = obj.fittingModelData.exoSignal.(exoSignalNameTemp);
                        exoSignalPart(:,cntExoSignal)   =   (obj.fittingModelData.snapShot.exoSignal.(exoSignalNameTemp).mask .* obj.fittingModelData.snapShot.exoSignal.(exoSignalNameTemp).Signal * params(obj.paramSet.flag.exoSignal(cntExoSignal)));
                    end  
                    
                    exoSignalPart   =   sum(exoSignalPart,2);
                end
                
                % Calculating short memory component
                if obj.paramSet.effect.memShort                                       
                    obj.fittingModelData.snapShot.memShort.Rate(isnan(obj.fittingModelData.snapShot.memShort.Rate))             = 0;
                    memShortPart     =   sum(obj.fittingModelData.snapShot.memShort.mask .* obj.fittingModelData.snapShot.memShort.Rate .*  repmat(params(obj.paramSet.flag.memShort),[totInput,1]),2);
                end
                
                % Calculating long memory component
                if obj.paramSet.effect.memSeason
                    obj.fittingModelData.snapShot.memSeason.Rate(isnan(obj.fittingModelData.snapShot.memSeason.Rate))             = 0;
                    memSeasonPart     =   sum(obj.fittingModelData.snapShot.memSeason.mask .* obj.fittingModelData.snapShot.memSeason.Rate .*  repmat(params(obj.paramSet.flag.memSeason),[totInput,1]),2);
                end
                
                % Calculating exogenous memory component
                if obj.paramSet.effect.memExo
                    exoMemNames = fieldnames(obj.fittingModelData.memExo);
                    cntExoMemTemp = 0;
                    for cntExoMem = 1:length(exoMemNames)
                        exoMemNameTemp = exoMemNames{cntExoMem};
                        if str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.memExo.(exoMemNameTemp).seqNumber))
                            exoMemTemp = obj.fittingModelData.snapShot.memExo.(exoMemNameTemp);
                            cntExoMemTemp = cntExoMemTemp + 1;
                            
                            memExoPart(:,cntExoMemTemp) = sum(repmat(exoMemTemp.meanData',[1, size(exoMemTemp.Signal,2)]) .*...
                                obj.component.memExo.(exoMemNameTemp).memoryFunction(params(obj.paramSet.flag.memExo((cntExoMemTemp-1)*2+1:(cntExoMemTemp-1)*2+2)),...
                                exoMemTemp.Time/1000) .*...
                                exoMemTemp.Signal .* ...
                                exoMemTemp.mask,2);
                            
                            if params(obj.paramSet.flag.memExo((cntExoMemTemp-1)*2+2))<0
                                rejectParams = 1;
                            end
                        end
                    end                                       
                    
                    memExoPart = sum(memExoPart,2);
                    
                    ARMeanPart   =   exoMemTemp.meanData' * params(obj.paramSet.flag.mean(end));
                end                                              
                
                memShortPart(isnan(memShortPart))       = 0;
                memSeasonPart(isnan(memSeasonPart))     = 0;
                memExoPart(isnan(memExoPart))           = 0;
                
                totOutputTemp = (biasPart + ARMeanPart + exoSignalPart + memShortPart + memSeasonPart + memExoPart)';
                
                if rejectParams
                    totOutputTemp = ones(size(totOutputTemp))*10^45;
                end
                ret = totOutputTemp;
                % toc
            end

        end
        
        % Compute the costfunction
        function ret = costfunction(obj, params, dataMask)
            if nargin < 3
                dataMask = logical(ones(size(obj.fittingModelData.output)))';
            end
            out_temp = obj.modelOutput(params, dataMask);            
            ret = obj.fittingModelData.output(dataMask') - out_temp;
        end  
        
        % ZScoreData
        function zscoreData(obj, dataMask)
            if nargin < 2
                dataMask = logical(ones(size(obj.fittingModelData.output)))';
            end
            % Exo Signals
            exoSignalNames = fieldnames(obj.fittingModelData.exoSignal);
            
            for cntExoSignal = 1:length(exoSignalNames)
                exoSignalNameTemp = exoSignalNames{cntExoSignal};
                obj.fittingModelData.snapShot.exoSignal.(exoSignalNameTemp).Signal = nanzscore(obj.fittingModelData.exoSignal.(exoSignalNameTemp).Signal(dataMask));                
                obj.fittingModelData.snapShot.exoSignal.(exoSignalNameTemp).mask = obj.fittingModelData.exoSignal.(exoSignalNameTemp).mask(dataMask);
            end
            
            % Short Memory
            obj.fittingModelData.snapShot.memShort.Rate = nanzscore(obj.fittingModelData.memShort.Rate(dataMask,:));
            obj.fittingModelData.snapShot.memShort.mask = (obj.fittingModelData.memShort.mask(dataMask,:));
            obj.fittingModelData.snapShot.memShort.Time = (obj.fittingModelData.memShort.Time(dataMask,:));
            
            % Seasonal Memory
            obj.fittingModelData.snapShot.memSeason.Rate = nanzscore(obj.fittingModelData.memSeason.Rate(dataMask,:));
            obj.fittingModelData.snapShot.memSeason.mask = (obj.fittingModelData.memSeason.mask(dataMask,:));
            obj.fittingModelData.snapShot.memSeason.Time = (obj.fittingModelData.memSeason.Time(dataMask,:));
            
            % Exogenous Memory
            exoMemNames = fieldnames(obj.fittingModelData.memExo);
            
            for cntExoMem = 1:length(exoMemNames)
                exoMemNameTemp = exoMemNames{cntExoMem};
                obj.fittingModelData.snapShot.memExo.(exoMemNameTemp).Signal = nanzscore(obj.fittingModelData.memExo.(exoMemNameTemp).Signal(dataMask,:));
                obj.fittingModelData.snapShot.memExo.(exoMemNameTemp).meanData = nanzscore(obj.fittingModelData.memExo.(exoMemNameTemp).meanData(dataMask));
                obj.fittingModelData.snapShot.memExo.(exoMemNameTemp).Time = (obj.fittingModelData.memExo.(exoMemNameTemp).Time(dataMask,:));
                obj.fittingModelData.snapShot.memExo.(exoMemNameTemp).mask = (obj.fittingModelData.memExo.(exoMemNameTemp).mask(dataMask,:));
            end
        end
        
        function saveModelSnapshot(obj,noDateTime)
            if nargin<2
                noDateTime = 0;
            end
            [~,~,~] = mkdir(obj.saveAddress);
            fittingResults = obj.fittingResults;
            if noDateTime
                save([obj.saveAddress, obj.saveNamePrefix,'_snapShot_Model_[',num2str(obj.fittingSet.currentModel),'].mat'],'fittingResults')
            else
                save([obj.saveAddress, obj.saveNamePrefix,'_snapShot_Model_[',num2str(obj.fittingSet.currentModel),']_',datestr(datetime),'.mat'],'fittingResults')
            end
        end
        
        function saveModel(obj,noDateTime)
            if nargin<2
                noDateTime = 0;
            end
            [~,~,~] = mkdir(obj.saveAddress);
            if noDateTime
                save([obj.saveAddress, obj.saveNamePrefix,'_all.mat'],'obj')
            else
                save([obj.saveAddress, obj.saveNamePrefix,'_all_',datestr(datetime),'.mat'],'obj')
            end
        end
        
        function chooseBestModel(obj)
            comparingFactor = 'RSS';
            comparingFunction = @max;
            comparingSet = [];
            for modelCnt = 1:length(obj.fittingResults.models)
                comparingSet(end+1) = obj.fittingResults.models{modelCnt}.res.stat.test{1}.(comparingFactor);
            end
            [V I] = comparingFunction(comparingSet);
            
            obj.fittingResults.bestModel = obj.fittingResults.models{I};
            obj.fittingResults.comparingSet = comparingSet;
            obj.fittingResults.bestModelNumber = I;
            
            % Convert best model into general format
            obj.fittingResults.bestModelGeneralFormat.components = dec2bin(I-1,length(obj.fittingSet.selectedComponents.currentSelectedSeq));
            obj.fittingResults.bestModelGeneralFormat.components = obj.fittingResults.bestModelGeneralFormat.components(end:-1:1);
            obj.fittingResults.bestModelGeneralFormat.params = nan(1,obj.fittingResults.models{end}.paramSet.flag.mean);
            obj.fittingResults.bestModelGeneralFormat.P_test = nan(1,obj.fittingResults.models{end}.paramSet.flag.mean);
            
            cellNames = fieldnames(obj.fittingResults.bestModel.paramSet.flag);
            for cntCell = [1:4 6]
                cellNameTemp = cellNames{cntCell};
                if ~isempty(obj.fittingResults.bestModel.paramSet.flag.(cellNameTemp))
                    obj.fittingResults.bestModelGeneralFormat.params(obj.fittingResults.models{end}.paramSet.flag.(cellNameTemp)) =...
                        obj.fittingResults.bestModel.res.vars(obj.fittingResults.bestModel.paramSet.flag.(cellNameTemp));
                    
                    obj.fittingResults.bestModelGeneralFormat.P_test(obj.fittingResults.models{end}.paramSet.flag.(cellNameTemp)) =...
                        obj.fittingResults.bestModel.res.stat.all.P_test(obj.fittingResults.bestModel.paramSet.flag.(cellNameTemp));
                end
            end
            
            cellNameTemp = 'memExo';
            exoMemNames = fieldnames(obj.fittingModelData.memExo);
            exoEffeciveNumTemp = 1; % store the numbers of effective exo memory components
            
            for cntN = 1:length(exoMemNames)
                exoMemNameTemp = exoMemNames{cntN};
                if str2num(obj.fittingResults.bestModelGeneralFormat.components(obj.fittingSet.selectedComponents.memExo.(exoMemNameTemp).seqNumber))
                    memExoTempIdx_Geneal = (cntN-1)*obj.paramSet.size.memExo+1:cntN*obj.paramSet.size.memExo;                    
                    
                    memExoTempIdx_Special = (exoEffeciveNumTemp-1)*obj.paramSet.size.memExo+1:exoEffeciveNumTemp*obj.paramSet.size.memExo;
                    
                    obj.fittingResults.bestModelGeneralFormat.params(obj.fittingResults.models{end}.paramSet.flag.(cellNameTemp)(memExoTempIdx_Geneal)) = ...
                        obj.fittingResults.bestModel.res.vars(obj.fittingResults.bestModel.paramSet.flag.(cellNameTemp)(memExoTempIdx_Special));
                    
                    obj.fittingResults.bestModelGeneralFormat.P_test(obj.fittingResults.models{end}.paramSet.flag.(cellNameTemp)(memExoTempIdx_Geneal)) = ...
                        obj.fittingResults.bestModel.res.stat.all.P_test(obj.fittingResults.bestModel.paramSet.flag.(cellNameTemp)(memExoTempIdx_Special));
                    
                    exoEffeciveNumTemp = exoEffeciveNumTemp + 1;
                end                                
            end
                       
        end
    end
    
    methods(Access = public)
        
        function tm = timeMaskBin(obj, binWindow)
            %   creates and returns timeMask (tm) object based on bin info
            
            tm.type = 1;  % 1: bin; 2: signal dependent            
            tm.binWindow = binWindow;
        end
        
        function tm = timeMaskSignal(obj, beginSig_Name, beginSig_TimeDist, endSig_Name, endSig_TimeDist)
            %   creates and returns timeMask (tm) object based on signals
            
            tm.type = 2;  % 1: bin; 2: signal dependent            
            tm.beginSig_Name = beginSig_Name;
            tm.beginSig_TimeDist = beginSig_TimeDist;
            
            tm.endSig_Name = endSig_Name;
            tm.endSig_TimeDist = endSig_TimeDist;
        end
        
        function addMemShort(obj, order, PACFTrans, effWin, paramLim)
            %   Initialize data structure for short memory
            %       .order:                 Indicates the memory order
            %
            %       .effectiveTimeWin:      Time window that this memory is
            %                               effective (refer to timeMask function)
            %
            %       .PACF_transfer:         Whether transfer AR Coeffs into
            %                               PACF (Partial AutoCorrelation Function) domain (Boolean)  
            %
            %       .paramLim:              Parameter fitting boundries
            
            if nargin < 2
                warning('Not enough arguments!');
                return
            elseif nargin < 3
                effWin = [];
                PACFTrans = false;
                paramLim = [-inf inf];
            elseif nargin < 4
                effWin = [];
                paramLim = [-inf inf];
            elseif nargin < 5                
                paramLim = [-inf inf];
            end
            
            proceedFlag = false;
            
            if obj.component.memShort.defined
                str = input('You already have a short memory component in your model.\n Do you want to rewrite the short memory? y/n [y]: ','s');
                if isempty(str)
                    str = 'y';
                end
                
                if strcmp(str,'y')
                    proceedFlag = true;
                end
            else
                proceedFlag = true;
            end
            
            if proceedFlag
                obj.component.memShort.defined = true;
                obj.component.memShort.order = order;
                obj.component.memShort.effectiveTimeWin = effWin;
                obj.component.memShort.PACF_transfer = PACFTrans;
                obj.component.memShort.paramLim = paramLim;
            end
        end
        
        function addMemSeason(obj, order, PACFTrans, effWin, paramLim)
            %   Initialize data structure for seasonal memory
            %       .order:                 Indicates the memory order
            %
            %       .effectiveTimeWin:      Time window that this memory is
            %                               effective (refer to timeMask function)
            %
            %       .PACF_transfer:         Whether transfer AR Coeffs into
            %                               PACF (Partial AutoCorrelation Function) domain (Boolean)       
            %
            %       .paramLim:              Parameter fitting boundries
            
            if nargin < 2
                warning('Not enough arguments!');
                return
            elseif nargin < 3
                effWin = [];
                PACFTrans = false;
                paramLim = [-inf inf];
            elseif nargin < 4
                effWin = [];
                paramLim = [-inf inf];
            elseif nargin < 5                
                paramLim = [-inf inf];
            end
            
            proceedFlag = false;
            
            if obj.component.memSeason.defined
                str = input('You already have a seasonal memory component in your model.\n Do you want to rewrite the seasonal memory? y/n [y]: ','s');
                if isempty(str)
                    str = 'y';
                end
                
                if strcmp(str,'y')
                    proceedFlag = true;
                end
            else
                proceedFlag = true;
            end
            
            if proceedFlag
                obj.component.memSeason.defined = true;
                obj.component.memSeason.order = order;
                obj.component.memSeason.effectiveTimeWin = effWin;
                obj.component.memSeason.PACF_transfer = PACFTrans;
                obj.component.memSeason.paramLim = paramLim;
            end
        end
                        
        function addMemExo(obj, name, signalNames, interaction, interactionType, signalTime, expo, expoOrder, memOrder, exoSelEff, effectiveTimeWin, memoryFunction, paramLim)
            %   Add and initialize data structure for exo memory
            %
            %           .name:                  Name of the exogenous memory
            %
            %           .signalName:        A set contains name(s) of the
            %                               signals involve in this memory       
            %
            %           .interaction:       Whether memory is based on
            %                               interaction of signals
            %
            %           .interactionType:   Type of interactio
            %                               1) simple product
            %                               2) seperate values
            %
            %           .signalTime:        Name of the 'signalTime'
            %                               corresponding for this exo memory                                      
            %
            %           .expo:              Whether memory function is a exponential
            %                               function
            %
            %           .expoOrder:         Order of exponential function
            %
            %           .memOrder:          Indicates the memory order                
            %
            %           .exoSelEff:         Whether consider the effect of exo
            %                               signal (Boolean)
            %                               We should have similar name exo
            %                               signal in "exoSignal" component
            %
            %           .effectiveTimeWin:  Time window that this memory is
            %                               effective (refer to timeMask function)
            % 
            %           .paramLim:          Parameter fitting boundries

            if nargin < 8
                warning('Not enough arguments!');
                return
            elseif nargin < 10
                effectiveTimeWin = [];
                exoSelEff = true;
                memoryFunction = @(x ,t) x(1).*(exp(-t./x(2)));
                paramLim = [-inf inf; -inf inf];
            elseif nargin < 11
                effectiveTimeWin = [];
                memoryFunction = @(x ,t) x(1).*(exp(-t./x(2)));
                paramLim = [-inf inf; -inf inf];
            elseif nargin < 12
                memoryFunction = @(x ,t) x(1).*(exp(-t./x(2)));
                paramLim = [-inf inf; -inf inf];
            elseif nargin < 13
                paramLim = [-inf inf; -inf inf];
            end
            
            if isempty(memoryFunction)
                memoryFunction = @(x ,t) x(1).*(exp(-t./x(2)));
            end
            
            proceedFlag = false;                       
                        
            if sum(strcmp(fieldnames(obj.component.memExo),name)) > 0
                str = input(['You already have a exo memory component with the name of "',name ,'" in your model.\n Do you want to rewrite this exo memory? y/n [y]: '],'s');
                if isempty(str)
                    str = 'y';
                end
                
                if strcmp(str,'y')
                    proceedFlag = true;
                end
            else
                proceedFlag = true;
            end
            
            if proceedFlag
                obj.component.memExo.(name).signalName = signalNames;
                obj.component.memExo.(name).interaction = interaction;
                obj.component.memExo.(name).interactionType = interactionType;
                obj.component.memExo.(name).signalTime = signalTime;
                obj.component.memExo.(name).expo = expo;
                obj.component.memExo.(name).expoOrder = expoOrder;
                obj.component.memExo.(name).memOrder = memOrder;
                obj.component.memExo.(name).exoSelEff = exoSelEff;
                obj.component.memExo.(name).effectiveTimeWin = effectiveTimeWin;
                obj.component.memExo.(name).memoryFunction = memoryFunction;
                obj.component.memExo.(name).paramLim = paramLim;
            end
        end
                        
        function addExoSignal(obj, name, signalNames, interaction, interactionType, effectiveTimeWin, paramLim)
            %   Add and initialize data structure for exo signal
            %       .name:                  Name of the exogenous memory
            %
            %           .signalName:        A set contains name(s) of the
            %                               signals involve in this memory
            %
            %           .interaction:       Whether memory is based on
            %                               interaction of signals
            %
            %           .interactionType:   Type of interactio
            %                               1) simple product
            %                               2) seperate values
            %
            %           .effectiveTimeWin:  Time window that this memory is
            %                               effective (refer to timeMask function)
            %
            %           .paramLim:          Parameter fitting boundries
            
            if nargin < 4
                warning('Not enough arguments!');
                return
            elseif nargin < 5
                interactionType = 2;
                effectiveTimeWin = [];
                paramLim = [-inf inf];
            elseif nargin < 6
                effectiveTimeWin = [];
                paramLim = [-inf inf];
            elseif nargin < 7
                paramLim = [-inf inf];
            end
            
            proceedFlag = false;                       
                        
            if sum(strcmp(fieldnames(obj.component.exoSignal),name)) > 0
                str = input(['You already have a exo signal component with the name of "',name ,'" in your model.\n Do you want to rewrite this exo signal? y/n [y]: '],'s');
                if isempty(str)
                    str = 'y';
                end
                
                if strcmp(str,'y')
                    proceedFlag = true;
                end
            else
                proceedFlag = true;
            end
            
            if proceedFlag
                obj.component.exoSignal.(name).signalName = signalNames;
                obj.component.exoSignal.(name).interaction = interaction;
                obj.component.exoSignal.(name).interactionType = interactionType;
                obj.component.exoSignal.(name).effectiveTimeWin = effectiveTimeWin;
                obj.component.exoSignal.(name).paramLim = paramLim;
            end
        end
    end
end

