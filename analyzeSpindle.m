%% load & process TD file
clear, clc
pathname = '~/LimbLab/Projects/S1-gamma/data';
% filename = 'Duncan_20190524_CObumpmove_10ms.mat';
filename = 'Han_20170203_COactpas_5ms.mat';
load([pathname filesep filename]);

% Smooth spikes, get split TD, get movement onsets, get rid of unsorted units

%Smooth spikes
smoothParams.signals = {'S1_spikes'};
smoothParams.width = 0.05;
smoothParams.calc_rate = true;
% trial_data.S1_spikes_bins = trial_data.S1_spikes;
td = smoothSignals(trial_data,smoothParams);

%Get speed
for i = 1:numel(td)
    td(i).speed = sqrt(td(i).vel(:,1).^2 + td(i).vel(:,2).^2);
end



% Smooth kinematic variables
% smoothParams.signals = {'spindle'};
% smoothParams.width = 0.01;
% smoothParams.calc_rate = false;
% td = smoothSignals(td,smoothParams);

%Get rid of unsorted units
sorted_idx = find(td(1).S1_unit_guide(:,2)~=0);

for i = 1:numel(td)
    td(i).S1_spikes = td(i).S1_spikes(:,sorted_idx);
    td(i).S1_unit_guide = td(i).S1_unit_guide(sorted_idx,:);
end

%Split TD
if numel(td) == 1
    splitParams.split_idx_name = 'idx_startTime';
    splitParams.linked_fields = {'result','bumpDir','tgtDir'};
    td = splitTD(td,splitParams);
elseif numel(td(1).tgtDir) > 1
    for i = 1:numel(td)
        td(i).tgtDir = td(i).tgtDir(i);
        td(i).result = td(i).result(i);
        td(i).bumpDir = td(i).bumpDir(i);
    end
end

%Get movement onset
moveOnsetParams.start_idx = 'idx_startTime';
moveOnsetParams.end_idx = 'idx_trial_end';
moveOnsetParams.which_field = 'speed';
td = getMoveOnsetAndPeak(td,moveOnsetParams);

% Separate TD into bump, act, and comb structures

% Get bump-move trials


% Get bump-only trials
td_temp = [];
for i = 1:numel(td)
    if ~isnan(td(i).idx_movement_on) && ~isnan(td(i).bumpDir) && td(i).idx_bumpTime < td(i).idx_goCueTime
        if strcmp(td(i).result,'R') || strcmp(td(i).result,'F')
            if isempty(td_temp)
                td_temp = td(i);
            else
                td_temp = [td_temp; td(i)];
            end
        end
    end
end

td_bump = td_temp; %td bump

% Get non-bump trials
td_temp = [];
for i = 1:numel(td)
    if ~isnan(td(i).idx_movement_on) && isnan(td(i).bumpDir)
        if strcmp(td(i).result,'R') && ~isnan(td(i).idx_goCueTime)
            if isempty(td_temp)
                td_temp = td(i);
            else
                td_temp = [td_temp; td(i)];
            end
        end
    end
end
td_act = td_temp; %td non-bump




%% Prep TDs for encoding/decoding models

srParams.powers = [0.5];
srParams.signals = {'vel'};
td_act = getPowerTD(td_act,srParams);


for i = 1:numel(td_act)
    idx = td_act(i).idx_movement_on:td_act(i).idx_movement_on + 130;
    td_act(i).pos = td_act(i).pos(idx,:);
    td_act(i).vel_0_5 = td_act(i).vel_0_5(idx,:);
    td_act(i).vel = td_act(i).vel(idx,:);
    td_act(i).speed = td_act(i).speed(idx,:);
    td_act(i).joint_vel = td_act(i).joint_vel(idx,:);
    td_act(i).musVelRel = td_act(i).musVelRel(idx,:);
    td_act(i).force = td_act(i).force(idx,:);
    td_act(i).spindle = td_act(i).spindle(idx,:);
    td_act(i).S1_spikes = td_act(i).S1_spikes(idx,:);
    
end

td_bump = getPowerTD(td_bump,srParams);

for i = 1:numel(td_bump)
    idx = td_bump(i).idx_bumpTime:td_bump(i).idx_bumpTime + 130;
    td_bump(i).pos = td_bump(i).pos(idx,:);
    td_bump(i).vel_0_5 = td_bump(i).vel_0_5(idx,:);
    td_bump(i).vel = td_bump(i).vel(idx,:);
    td_bump(i).speed = td_bump(i).speed(idx,:);
    td_bump(i).joint_vel = td_bump(i).joint_vel(idx,:);
    td_bump(i).musVelRel = td_bump(i).musVelRel(idx,:);
    td_bump(i).force = td_bump(i).force(idx,:);
    td_bump(i).spindle = td_bump(i).spindle(idx,:);
    td_bump(i).S1_spikes = td_bump(i).S1_spikes(idx,:);
end





%%

clearvars -except td_act td_bump

modelParams.nFolds = 10;
modelParams.model_name = 'spindle';
modelParams.out_signals = {'vel','all'};
modelParams.in_signals = {'S1_spikes','all'};
modelParams.model_type = 'linmodel';
modelParams.eval_metric = 'vaf';


td_bump = td_xval(td_bump,modelParams);
[~, decoder_bump] = getModel(td_bump,modelParams); %Get model trained on bumps
test_bump2bump_spindle = evalModel(td_bump,modelParams); 

td_act = getModel(td_act,decoder_bump);
test_bump2act_spindle = evalModel(td_act,modelParams);



modelParams.in_signals = {'musVelRel','all'};
modelParams.model_name = 'muscle_vel';

td_bump = td_xval(td_bump,modelParams);
[~, decoder_bump] = getModel(td_bump,modelParams); %Get model trained on bumps
test_bump2bump_musVel = evalModel(td_bump,modelParams); 

td_act = getModel(td_act,decoder_bump);
test_bump2act_musVel = evalModel(td_act,modelParams);



plot(test_bump2act_spindle,test_bump2act_musVel,'.')
axis([-1 1 -1 1])




