%% load TD file
pathname = 'C:\Users\kpb8927\data\td-library';
filename = 'Duncan_20190515_CObumpmove_10ms.mat';
load([pathname filesep filename]);

%% Smooth spikes, get split TD, get movement onsets, get rid of unsorted units

%Smooth spikes
smoothParams.signals = {'S1_spikes'};
smoothParams.width = 0.05;
smoothParams.calc_rate = true;
td = smoothSignals(trial_data,smoothParams);

%Get speed
td.speed = sqrt(td.vel(:,1).^2 + td.vel(:,2).^2);

%Get rid of unsorted units
sorted_idx = find(td.S1_unit_guide(:,2)~=0);
td.S1_spikes = td.S1_spikes(:,sorted_idx);
td.S1_unit_guide = td.S1_unit_guide(sorted_idx,:);

%Split TD
splitParams.split_idx_name = 'idx_startTime';
splitParams.linked_fields = {'result','bumpDir'};
td = splitTD(td,splitParams);

%Get movement onset
moveOnsetParams.start_idx = 'idx_startTime';
moveOnsetParams.end_idx = 'idx_endTime';
td = getMoveOnsetAndPeak(td,moveOnsetParams);

td_temp = [];

for i = 1:numel(td)
    if ~isnan(td(i).idx_movement_on) && ~isnan(td(i).bumpDir)
        if strcmp(td(i).result,'R') || strcmp(td(i).result,'F')
            if isempty(td_temp)
                td_temp = td(i);
            else
                td_temp = [td_temp; td(i)];
            end
        end
    end
end

td = td_temp;

%% Plot trial averaged firing rates for each bump direction
params.savepath = '~/LimbLab/Projects/EMGanalysis/Figures/Han/Han_20171207_COactpas/';
params.filetype = '.pdf';
params.savefig = 0;
params.trialType = 'bump';
params.spikeArray = 'S1_spikes';


bumpDirs = -135:45:180;
for unit = 51:60
    
    
    
    
    params.frToPlot = unit;
    plotFR(td,params)
    
    
    
    
end


%% Get tuning curves for passive bumps
bumpStruct = struct('S1_spikes',[],'bumpDir',[],'vel',[],'monkey',[],'task',[],'date',[]);
for i = 1:numel(td)
    if ~isnan(td(i).bumpDir)
        bumpIdx = td(i).idx_bumpTime:(td(i).idx_bumpTime+12);
        if i ==1
            bumpStruct(end).S1_spikes = td(i).S1_spikes(bumpIdx,:);
        else
            bumpStruct(end+1).S1_spikes = td(i).S1_spikes(bumpIdx,:);
        end
        bumpStruct(end).bumpDir = td(i).bumpDir;
        bumpStruct(end).vel = td(i).vel(bumpIdx,:);
        bumpStruct(end).monkey = td(i).monkey;
        bumpStruct(end).task = td(i).task;
        bumpStruct(end).date = td(i).date_time(1:9);
    end
end

PDparams.out_signals = 'S1_spikes';
PDparams.in_signals = 'vel';
PDparams.numBoots = 1;

pdTable =  getTDPDs(bumpStruct,PDparams);


%% Plot PDs
velPDs = pdTable.velPD(:);
velMDs = pdTable.velModdepth(:);
figure; polarhistogram(velPDs,30);
figure; subplot(3,3,[2 3 5 6]); hold on;
plot(velPDs,velMDs,'.');

subplot(3,3,[1 4]); hold on;
[counts,bins] = hist(velMDs);
barh(bins,counts);

subplot(3,3,[8 9]); hold on;
[counts,bins] = hist(velPDs);
bar(bins,counts);

