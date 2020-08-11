% This file shows how to use the ARMAXNeuro framework
%
% You need to fetch the data for neuron and rearrengment it in proper format
% After initiaing ARMAXNeuro model, you can feed the neural data for
% further analysis.
%
% Copyright 2019 Mehran M. Spitmaan (mehran.m.spitman@gmail.com).
% 
% #########################################################################
% ######################### SAMPLE NEURAL DATA FOR SINGLE NEURON ##########
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

clc
clear all

anonymTools;


%% Loading Sample neural data
load('neuronDataSample.mat')

%% Setting up the sample neuralData Structure
neuronData                             = [];

neuronData.spikeTime                   = neuronDataSample.spikeTime;

neuronData.signalTime.feedback         = neuronDataSample.feedbackStart;
neuronData.signalTime.reward           = neuronDataSample.rewardStart;
neuronData.signalTime.choice           = neuronDataSample.choiceStart;
neuronData.signalTime.target           = neuronDataSample.targetStart;
neuronData.signalTime.magnitude        = neuronDataSample.magStart;
neuronData.signalTime.gaze             = neuronDataSample.gazeStart;
neuronData.signalTime.fixation         = neuronDataSample.fixationStart;

neuronData.signalValue.reward          = neuronDataSample.trialData(:,6)*2-1;
neuronData.signalValue.choice          = neuronDataSample.trialData(:,4)*2-1;
neuronData.signalValue.choicecolor     = neuronDataSample.trialData(:,5)*2-1;

neuronData.endOfTrialSignal            = 'feedback';
neuronData.endOfTrialDist              = 500;
neuronData.maxTrialLen                 = 6500;
neuronData.beginOfTrialSignal          = 'fixation';
neuronData.binSize                     = 50;
neuronData.trialRange                  = neuronDataSample.trialRange;

%% Initiating an ARMAXNeuro model

% Adding Input data
ARMAXNeuroModel = ARMAXNeuro(neuronData, 1, './MyFolder/', 'MyNeuralModel');

% Create a time mask for choice from 0 to 500 msec around choice signal:
% - BeginSignalName: "choice";
% - beginSignalTimeDistance: 0;
% - EndSignalName: "choice";
% - EndSignalTimeDistance: 500.
tm_choice1 = ARMAXNeuroModel.timeMaskSignal('choice', 0, 'choice', 500);

% Create a time mask for reward choice interaction from 0 to 500 msec around feedback signal:
% - BeginSignalName: "feedback";
% - beginSignalTimeDistance: 0;
% - EndSignalName: "feedback";
% - EndSignalTimeDistance: 500.
tm_reward_x_choice = ARMAXNeuroModel.timeMaskSignal('feedback', 0, 'feedback', 500);

% Adding a short memory component with the depth of 5
ARMAXNeuroModel.addMemShort(5,1,tm_choice1);

% Adding a seasonal memory component with the depth of 5
ARMAXNeuroModel.addMemSeason(5,1,tm_choice1);

% Adding an exogenous memory component with information below:
% - Name: "reward"; 
% - No interactios of signals; 
% - SignalName: "reward";
% - SignalTime: "feedback"; 
% - Single exponential memory; 
% - Depth of 5;
% - Considering the effect of exoSignal with the same name as "reward";
% - No time window limit.
ARMAXNeuroModel.addMemExo('reward', 'reward', 0, [], 'feedback', 1, 1, 5, 1, [], [], [-4 4; -eps 4*20+eps]);

% Adding an exogenous memory component with information below:
% - Name: "choice1"; 
% - No interactios of signals; 
% - SignalName: "choice";
% - SignalTime: "choice"; 
% - Single exponential memory; 
% - Depth of 5;
% - Considering the effect of exoSignal with the same name as "choice";
% - Timewindow of [0:500] around choice, defined by "tm_choice1".
ARMAXNeuroModel.addMemExo('choice1', 'choice', 0, [], 'choice', 1, 1, 5, 1, tm_choice1, [], [-4 4; -eps 4*20+eps]);

% Adding an exogenous signal component with information below:
% - Name: "reward"; 
% - No interactios of signals; 
% - SignalName: "reward";
% - No time window limit.
ARMAXNeuroModel.addExoSignal('reward', 'reward', 0, [], []);

% Adding an exogenous signal component with information below:
% - Name: "choice1"; 
% - No interactios of signals; 
% - SignalName: "reward";
% - Timewindow of [0:500] around choice, defined by "tm_choice1".
ARMAXNeuroModel.addExoSignal('choice1', 'choice', 0, [], tm_choice1);

% Adding an exogenous signal component with information below:
% - Name: "reward_x_choice"; 
% - Interaction of "reward" and "choice" signals;
% - Interaction type (2) seperate values
% - Timewindow of [0:500] around feedback, defined by "tm_reward_x_choice".
ARMAXNeuroModel.addExoSignal('reward_x_choice', {'reward','choice'}, 1, 2, tm_reward_x_choice);

% Initializing the data model
tic
ARMAXNeuroModel.initData();
toc

%% Fitting the Model

%  Setting fitting options
iterNum = 3; % number of initial points || 30
crossPer = 0.1; % Cross-validation test data percentage
crossIter = 5; % Cross-validation instances || 50
modelSelection = 0; % Whether do model selection or not

ARMAXNeuroModel.setFittingSet(iterNum, crossPer, crossIter, modelSelection);

% Fitting the model
ARMAXNeuroModel.fit();
