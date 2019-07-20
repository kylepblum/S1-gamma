clear, clc
filename = 'C:\Users\kpb8927\data\td-library\td-spindle\Han_20171106_TRT_5ms.mat';
load(filename)

modelNames = {'muscleLen','muscleVel','spindle','emg'};
nFolds = 5;




%% Remove buffer for spindle data & process spikes
splitParams.split_idx_name = 'idx_startTime';
splitParams.extra_bins = [20 20];
trial_data = splitTD(trial_data,splitParams);

smoothParams.signals = {'S1_spikes'};
smoothParams.calc_rate = true;
smoothParams.width = 0.050;
trial_data = smoothSignals(trial_data,smoothParams);

smoothParams.signals = {'spindle','emgNorm','muscle_len','muscle_vel'};
smoothParams.calc_rate = false;
trial_data = smoothSignals(trial_data,smoothParams);

%Get rid of unsorted units / apply lags for different variables
sorted_idx = find(trial_data(1).S1_unit_guide(:,2)~=0);

for trial = 1:numel(trial_data)
    trial_data(trial).S1_spikes = trial_data(trial).S1_spikes(:,sorted_idx);
    trial_data(trial).S1_unit_guide = trial_data(trial).S1_unit_guide(sorted_idx,:);
    
    trial_data(trial).musVelRel = trial_data(trial).musVelRel(21:10:end-21,:);
    trial_data(trial).musLenRel = trial_data(trial).musLenRel(21:10:end-21,:);
    trial_data(trial).muscle_vel = trial_data(trial).muscle_vel(21:10:end-21,:);
    trial_data(trial).muscle_len = trial_data(trial).muscle_len(21:10:end-21,:);
    trial_data(trial).spindle = trial_data(trial).spindle(21:10:end-21,:); 
    trial_data(trial).emgNorm = trial_data(trial).emgNorm(21:10:end-21,:);
    trial_data(trial).S1_spikes = trial_data(trial).S1_spikes(21:10:end-21,:);
end
%% Fit GLMs to neural data and cross-validate
clc
nTrials = numel(trial_data);
trials = 1:nTrials;
randTrials = randperm(nTrials);
blockSize = ceil(nTrials/nFolds);

temp = trial_data;
test_out = [];
for fold = 1:nFolds
    disp(['Now permorming fold ' num2str(fold) ' out of ' num2str(nFolds)])
    blockIdx = (fold-1)*blockSize < trials & fold*blockSize >= trials;
    train_idx = randTrials(~blockIdx);
    test_idx = randTrials(blockIdx);
    
    td_train = trial_data(train_idx);
    td_test{fold} = trial_data(test_idx);
    
    for model = 1:numel(modelNames)
        modelParams.model_name = modelNames{model};
        modelParams.model_type = 'glm';
        modelParams.out_signals = {'S1_spikes','all'};
        switch model
            case 1
                modelParams.in_signals = {'musLenRel','all'};
            case 2
                modelParams.in_signals = {'musVelRel','all'};
            case 3
                modelParams.in_signals = {'spindle','all'};
            case 4 
                modelParams.in_signals = {'emgNorm','all'};
            case 5
        end
        
        [~, model_info] = getModel(td_train,modelParams);
        td_test{fold} = getModel(td_test{fold}, model_info);

        %Evaluate models
        modelParams.eval_metric = 'pr2';
        modelParams.model_name = modelNames{model};
        test(model,fold,:) = evalModel(td_test{fold},modelParams);

    end
end


%% Evaluate models
clc

modelParams.eval_metric = 'vaf';

for model = 1:numel(modelNames)
    modelParams.model_name = modelNames{model};
    switch model
        case 1
            modelParams.in_signals = {'musLenRel','all'};
        case 2
            modelParams.in_signals = {'spindle','all'};
        case 3
            modelParams.in_signals = {'emgNorm','all'};
        case 4
            modelParams.in_signals = {'emgNormAFF','all'};
    end
    test(model,fold,:) = evalModel(td_test{fold},modelParams);
end


