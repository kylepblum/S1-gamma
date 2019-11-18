%% load & process TD file
pathname = '~/LimbLab/Projects/S1-gamma/data';

filename = 'Duncan_20191111_OOR_10ms.mat';

load([pathname filesep filename]);


% Smooth spikes, get split TD, get movement onsets, get rid of unsorted units

%Smooth spikes
smoothParams.signals = {'S1_spikes'};
smoothParams.width = 0.03;
smoothParams.calc_rate = true;
trial_data.S1_spikes_bins = trial_data.S1_spikes;
td = smoothSignals(trial_data,smoothParams);


%Get speed
td.speed = sqrt(td.vel(:,1).^2 + td.vel(:,2).^2);

%Get accel
td.acc = diff(td.vel)./td.bin_size;
td.acc(end+1,:) = td.acc(end,:);

%Remove offset
td.pos(:,1) = td.pos(:,1)+0;
td.pos(:,2) = td.pos(:,2)+32;

% Smooth kinematic variables
smoothParams.signals = {'pos','vel','acc','force'};
smoothParams.width = 0.03;
smoothParams.calc_rate = false;
td = smoothSignals(td,smoothParams);

%Get rid of unsorted units
sorted_idx = find(td.S1_unit_guide(:,2)~=0);
td.S1_spikes = td.S1_spikes(:,sorted_idx);
td.S1_spikes_bins = td.S1_spikes_bins(:,sorted_idx);
td.S1_unit_guide = td.S1_unit_guide(sorted_idx,:);

%Split TD
splitParams.split_idx_name = 'idx_startTime';
splitParams.linked_fields = {'result','forceDir','tgtDir'};
td = splitTD(td,splitParams);

%Get movement onset
moveOnsetParams.start_idx = 'idx_startTime';
moveOnsetParams.end_idx = 'idx_endTime';
td = getMoveOnsetAndPeak(td,moveOnsetParams);

for i = 1:numel(td)
    if td(i).forceDir == 360
        td(i).forceDir = td(i).forceDir - 360;
    end
end

td_temp = [];
for i = 1:numel(td)
    if ~isnan(td(i).idx_movement_on) && abs(td(i).forceDir-td(i).tgtDir)==180
        if strcmp(td(i).result,'R') || strcmp(td(i).result,'F')
            if isempty(td_temp)
                td_temp = td(i);
            else
                td_temp = [td_temp; td(i)];
            end
        end
    end
end

td_resist = td_temp;


td_temp = [];
for i = 1:numel(td)
    if ~isnan(td(i).idx_movement_on) && abs(td(i).forceDir-td(i).tgtDir)==0
        if strcmp(td(i).result,'R') || strcmp(td(i).result,'F')
            if isempty(td_temp)
                td_temp = td(i);
            else
                td_temp = [td_temp; td(i)];
            end
        end
    end
end

td_assist = td_temp;






%% Prep TDs for encoding/decoding models
modelNames = {'act','pas','comb'};

srParams.powers = [0.5];
srParams.signals = {'vel'};

for i = 1:numel(td_resist)
    idx = td_resist(i).idx_movement_on+0:td_resist(i).idx_movement_on+50;
    td_resist(i).pos = td_resist(i).pos(idx,:);
    td_resist(i).vel = td_resist(i).vel(idx,:);
    td_resist(i).acc = td_resist(i).acc(idx,:);
    td_resist(i).speed = td_resist(i).speed(idx,:);
    td_resist(i).force = td_resist(i).force(idx,:);
    td_resist(i).S1_spikes = td_resist(i).S1_spikes(idx,:);
    td_resist(i).meanTraj = wrapTo2Pi(atan2(mean(td_resist(i).vel(:,2)),mean(td_resist(i).vel(:,1))));
    td_resist(i).forceTraj = wrapTo2Pi(atan2(mean(td_resist(i).force(:,2)),mean(td_resist(i).force(:,1))));
    td_resist(i).accTraj = wrapTo2Pi(atan2(mean(td_resist(i).acc(:,2)),mean(td_resist(i).acc(:,1))));
end


for i = 1:numel(td_assist)
    idx = td_assist(i).idx_movement_on+0:td_assist(i).idx_movement_on+50;
    td_assist(i).pos = td_assist(i).pos(idx,:);
    td_assist(i).vel = td_assist(i).vel(idx,:);
    td_assist(i).acc = td_assist(i).acc(idx,:);
    td_assist(i).speed = td_assist(i).speed(idx,:);
    td_assist(i).force = td_assist(i).force(idx,:);
    td_assist(i).S1_spikes = td_assist(i).S1_spikes(idx,:);
    td_assist(i).meanTraj = wrapTo2Pi(atan2(mean(td_assist(i).vel(:,2)),mean(td_assist(i).vel(:,1))));
    td_assist(i).forceTraj = wrapTo2Pi(atan2(mean(td_assist(i).force(:,2)),mean(td_assist(i).force(:,1))));
    td_assist(i).accTraj = wrapTo2Pi(atan2(mean(td_assist(i).acc(:,2)),mean(td_assist(i).acc(:,1))));
end

%% Plot S1 FRs
clearvars -except td_assist td_resist
clc


plotTuning = 1;
meanFRs_RES = [];
meanFRs_ASS = [];

for context = 1:2
    clear meanFRs_1 bumpDir_1 smoothFRs1 smoothSTD1 meanSens
    if context == 1
        td_1 = td_resist;
        color1 = [0.9 0.2 0.2]; % RED
    elseif context == 2
        td_1 = td_assist;
        color1 = [0.2 0.2 0.7]; % BLUE
    end
    
    
    nCell = size(td_1(1).S1_spikes,2);
    for cell = 1:nCell
        for i = 1:numel(td_1)
            meanFRs_1(cell,i) = mean(td_1(i).S1_spikes(:,cell));
            meanSens(cell,i) = mean(td_1(i).S1_spikes(:,cell)./td_1(i).speed);
            
            if context == 1
                bumpDir_1(cell,i) = wrapTo2Pi(td_1(i).meanTraj);
            elseif context == 2
                bumpDir_1(cell,i) = wrapTo2Pi(td_1(i).meanTraj);
            end
        end
        [bumpDir_1(cell,:),I1] = sort(bumpDir_1(cell,:));
        meanFRs_1(cell,:) = meanFRs_1(cell,I1);
        smoothFRs1(cell,:) = smooth(meanFRs_1(cell,:),10);
        smoothSTD1(cell,:) = movstd(meanFRs_1(cell,:),10);
    end
    smoothFRs1(:,end+1) = smoothFRs1(:,1);
    smoothSTD1(:,end+1) = smoothSTD1(:,1);
    meanFRs_1(:,end+1) = meanFRs_1(:,1);
    meanSens(:,end+1) = meanSens(:,1);
    bumpDir_1(:,end+1) = bumpDir_1(:,1);
%     
    for cell = 1:nCell
        if context == 1 && plotTuning == 1
            hfig(cell) = figure(cell); set(gcf,'Color','white');
            hfig(cell).Position(3) = hfig(cell).Position(3)*1.5;
            hfig(cell).Position = hfig(cell).Position/3;

            a1(cell) = polaraxes;
            set(a1(cell),'FontName','Helvetica','FontSize',9,...
                'clipping','off'); hold on
            a1(cell).Position(1) = a1(cell).Position(1)-0.25;
            a2(cell) = polaraxes;
            set(a2(cell),'FontName','Helvetica','FontSize',9,...
                'clipping','off'); hold on
            a2(cell).Position(1) = a2(cell).Position(1)+0.25;

        end
        edges = 0:pi/4:2*pi;
        Y = discretize(bumpDir_1(cell,:),edges);
        for i = 1:numel(edges)
            FRs_bins(i) = mean(meanFRs_1(cell,Y==i));
            FRs_ct(i) = sum(Y==i);
            FRs_sem(i) = std(meanFRs_1(cell,Y==i));
            
            Sens_bins(i) = mean(meanSens(cell,Y==i));
            Sens_sem(i) = std(meanSens(cell,Y==i));            
        end
        
        FRs_ct(isnan(FRs_bins)) = [];
        FRs_sem(isnan(FRs_bins)) = [];
        Sens_bins(isnan(FRs_bins)) = [];
        Sens_sem(isnan(FRs_bins)) = [];
        edges(isnan(FRs_bins)) = [];
        FRs_bins(isnan(FRs_bins)) = [];
        
        edges(end+1) = edges(1);
        FRs_bins(end+1) = FRs_bins(1);
        FRs_sem(end+1) = FRs_sem(1);
        FRs_ct(end+1) = FRs_ct(1);
        Sens_bins(end+1) = Sens_bins(1);
        Sens_sem(end+1) = Sens_sem(1);
        
        if plotTuning == 1
            polarplot(a1(cell),edges,FRs_bins,'color',color1,'LineWidth',2)
            polarplot(a1(cell),edges,FRs_bins-FRs_sem./FRs_ct,'color',color1)
            polarplot(a1(cell),edges,FRs_bins+FRs_sem./FRs_ct,'color',color1)
            
            [maxFR,maxI] = max(FRs_bins);
            polarplot(a1(cell),[edges(maxI) edges(maxI)],[0 30],'color',color1,'LineWidth',1.5)
            
            polarplot(a2(cell),edges,Sens_bins,'color',color1,'LineWidth',2)
            polarplot(a2(cell),edges,Sens_bins-Sens_sem./FRs_ct,'color',color1)
            polarplot(a2(cell),edges,Sens_bins+Sens_sem./FRs_ct,'color',color1)
        end
        
        switch context
            case 1
                meanFRs_RES(end+1,:) = Sens_bins;
            case 2
                meanFRs_ASS(end+1,:) = Sens_bins;
        end
        
        clear FRs_bins FRs_sem clear FRs_ct
    end
end