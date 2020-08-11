# NeuroARMAX_LeeLabData
Here, you can find neural data and computer codes required for model fitting and data analyses presented in Spitmaan et al.

This a ARMAX framework plus exponential memory trace components, 
in order to estimate short, seasonal, and exogenous (reward and choice) memories 
plus exogenous selctivity of neurons, based on the activity of neurons in any task.

## Initiation and Loading Data
In order to initiate an instance of the model you need feed it with neural data.
```
ARMAXNeuroModel = ARMAXNeuro(neuronData);
```

## Neural Data Structure
The neural data set (in this case neuronData data structure) should follow a specific structure of data.
To see an example you can refer to "Setting up the sample neuralData Structure" in [RunMe.m](https://github.com/Spitmaan/ARMAXNeuro/blob/master/RunMe.m)

Here is the list of fields in neuronData structure:
neuronData (Struct)
- .spikeTime:         All the spike time for entire experiment
- .signalTime:        Structure contain absolute time of all signal for each trial. Name of each field represents name of each signal.
- .signalValue:       Structure contain values of all the possible signal for each trial. There might be signals with signalTime but without signalValue. Name of each field represents name of each signal.
- .endOfTrialSignal:  Name of the signal that indicates the end of the trial. There might be a time distance before/after "endOfTrialSignal" to the real end of the trial which can be set by "endOfTrialDist". (String).
- .endOfTrialDist:    Time distanse from "endOfTrialSignal" that indicated the real end of the trial. (msec)
- .maxTrialLen:       Maximum time for each trial. Useful for arranging data into bins. (msec) The real value for trial length would be the caculated by 
```
"(signalTime.(endOfTrialSignal)(currTrial) +
endOfTrialDist) - ...
(signalTime.(endOfTrialSignal)(currTrial-1) +
endOfTrialDist);
```
- .beginOfTrialSignal:    Name of the signal that indicates the begin of the trial. (String)
- .binSize:           Size of each bin (msec)
- .trialRange:        Range of trials; Total number of trials.


## Adding Components into the Model
Here you can find the example codes in order to add desired model components into your defined instanse (ARMAXNeuroModel):
### Time Mask for Exogenous Signals
Create a time mask for choice from 0 to 500 msec around choice signal:
- BeginSignalName: "choice";
- beginSignalTimeDistance: 0;
- EndSignalName: "choice";
- EndSignalTimeDistance: 500.
```
tm_choice1 = ARMAXNeuroModel.timeMaskSignal('choice', 0, 'choice', 500);
```
Create a time mask for reward choice interaction from 0 to 500 msec around feedback signal:
- BeginSignalName: "feedback";
- beginSignalTimeDistance: 0;
- EndSignalName: "feedback";
- EndSignalTimeDistance: 500.
```
tm_reward_x_choice = ARMAXNeuroModel.timeMaskSignal('feedback', 0, 'feedback', 500);
```
### Intrinsic Memory
Adding an intrinsic memory component with the depth of 5
```
ARMAXNeuroModel.addMemShort(5);
```
### Seasonal Memory
Adding a seasonal memory component with the depth of 5
```
ARMAXNeuroModel.addMemSeason(5);
```
### Exogenous Memories
Adding an exogenous memory component with information below:
- Name: "reward"; 
- No interactios of signals; 
- SignalName: "reward";
- SignalTime: "feedback"; 
- Single exponential memory; 
- Depth of 5;
- Considering the effect of exoSignal with the same name as "reward";
- No time window limit.
```
ARMAXNeuroModel.addMemExo('reward', 'reward', 0, 'feedback', 1, 1, 5, 1, []);
```
Adding an exogenous memory component with information below:
- Name: "choice1"; 
- No interactios of signals; 
- SignalName: "choice";
- SignalTime: "choice"; 
- Single exponential memory; 
- Depth of 5;
- Considering the effect of exoSignal with the same name as "choice";
- Timewindow of [0:500] around choice, defined by "tm_choice1".
```
ARMAXNeuroModel.addMemExo('choice1', 'choice', 0, 'choice', 1, 1, 5, 1, tm_choice1);
```
### Exogenous Signals
Adding an exogenous signal component with information below:
- Name: "reward"; 
- No interactios of signals; 
- SignalName: "reward";
- No time window limit.
```
ARMAXNeuroModel.addExoSignal('reward', 'reward', 0, []);
```
Adding an exogenous signal component with information below:
- Name: "choice1"; 
- No interactios of signals; 
- SignalName: "reward";
- Timewindow of [0:500] around choice, defined by "tm_choice1".
```
ARMAXNeuroModel.addExoSignal('choice1', 'choice', 0, tm_choice1);
```
Adding an exogenous signal component with information below:
- Name: "reward_x_choice"; 
- Interaction of "reward" and "choice" signals;
- Timewindow of [0:500] around feedback, defined by "tm_reward_x_choice".
```
ARMAXNeuroModel.addExoSignal('reward_x_choice', {'reward','choice'}, 1, tm_reward_x_choice);
```

## Fitting the Model
Here you can find the example codes in order to prepear, setup and fit the predefined model instanse (ARMAXNeuroModel):
### Initializing Model Data
Initialize the data model based on introduced components and dataset:
```
ARMAXNeuroModel.initData();
```
### Setting Fitting Options
Setting a fitting procedure with information below:
- iterNum: Number of initial points
- crossPer: Cross-validation test data percentage
- crossIter: Cross-validation instances
```
ARMAXNeuroModel.setFittingSet(iterNum, crossPer, crossIter);
```
### Fitting the Model
Fitting the model based on introduced (or default) options.
```
ARMAXNeuroModel.fit();
```

## Misc.
CCNL, Dartmouth College, (c) 2020
