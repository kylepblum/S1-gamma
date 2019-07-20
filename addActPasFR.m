%% load & process TD file
pathname = 'C:\Users\kpb8927\data\td-library';
filename = 'Duncan_20190524_CObumpmove_10ms.mat';
load([pathname filesep filename]);

% Smooth spikes, get split TD, get movement onsets, get rid of unsorted units

%Smooth spikes
smoothParams.signals = {'S1_spikes'};
smoothParams.width = 0.05;
smoothParams.calc_rate = true;
td = smoothSignals(trial_data,smoothParams);

%Get speed
td.speed = sqrt(td.vel(:,1).^2 + td.vel(:,2).^2);

%Remove offset
td.pos(:,1) = td.pos(:,1)+0;
td.pos(:,2) = td.pos(:,2)+32;
%Get rid of unsorted units
sorted_idx = find(td.S1_unit_guide(:,2)~=0);
td.S1_spikes = td.S1_spikes(:,sorted_idx);
td.S1_unit_guide = td.S1_unit_guide(sorted_idx,:);

%Split TD
splitParams.split_idx_name = 'idx_startTime';
splitParams.linked_fields = {'result','bumpDir','tgtDir'};
td = splitTD(td,splitParams);

%Get movement onset
moveOnsetParams.start_idx = 'idx_startTime';
moveOnsetParams.end_idx = 'idx_endTime';
td = getMoveOnsetAndPeak(td,moveOnsetParams);

td_temp = [];
