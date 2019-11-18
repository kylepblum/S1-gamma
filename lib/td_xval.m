function trial_data = td_xval(trial_data,params)

randTrials = randperm(numel(trial_data));
trials = 1:numel(trial_data);
blockSize = ceil(numel(trials)/params.nFolds);

for thisFold = 1:params.nFolds
    disp(['Now performing fold ' num2str(thisFold) ' out of ' num2str(params.nFolds)])
    blockIdx = (thisFold-1)*blockSize < trials & thisFold*blockSize >= trials;
    train_idx = randTrials(~blockIdx);
    test_idx = randTrials(blockIdx);
    
    td_train{thisFold} = trial_data(train_idx);
    td_test{thisFold} = trial_data(test_idx);
    
    [~, model_info] = getModel(td_train{thisFold},params);
    td_test{thisFold} = getModel(td_test{thisFold}, model_info);
end

td_test = vertcat(td_test{:});
trial_data = td_test(randTrials);

end