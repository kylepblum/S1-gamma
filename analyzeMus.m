%% load & process TD file
clear, clc
pathname = '~/LimbLab/Projects/S1-gamma/data';
% filename = 'Duncan_20190524_CObumpmove_10ms.mat';
filename = 'Han_20171101_TRT_5ms.mat';
load([pathname filesep filename]);

% Smooth spikes, get split TD, get movement onsets, get rid of unsorted units

%Smooth spikes
smoothParams.signals = {'S1_spikes'};
smoothParams.width = 0.03;
smoothParams.calc_rate = true;
td = smoothSignals(trial_data,smoothParams);


%Get speed
td.speed = sqrt(td.vel(:,1).^2 + td.vel(:,2).^2);

% %Get accel
% td.acc = diff(td.vel)./td.bin_size;
% td.acc(end+1,:) = td.acc(end,:);

%Remove offset
td.pos(:,1) = td.pos(:,1)+0;
td.pos(:,2) = td.pos(:,2)+32;

%Get norm muscle signals
td = getNormEMG(td,[]);


% Smooth kinematic variables
smoothParams.signals = {'pos','vel','emg','emgNorm',...
    'muscle_len','muscle_vel'};
smoothParams.width = 0.01;
smoothParams.calc_rate = false;
td = smoothSignals(td,smoothParams);

%Get rid of unsorted units
sorted_idx = find(td.S1_unit_guide(:,2)~=0);
td.S1_spikes = td.S1_spikes(:,sorted_idx);
% td.S1_spikes_bins = td.S1_spikes_bins(:,sorted_idx);
td.S1_unit_guide = td.S1_unit_guide(sorted_idx,:);

%Split TD
splitParams.split_idx_name = 'idx_startTime';
splitParams.linked_fields = {'result','bumpDir','spaceNum','trialID'};
td = splitTD(td,splitParams);

%Get movement onset
moveOnsetParams.start_idx = 'idx_startTime';
moveOnsetParams.end_idx = 'idx_endTime';
td = getMoveOnsetAndPeak(td,moveOnsetParams);

% Separate TD into bump, act, and comb structures


% Get bump-only trials
td_temp = [];
for i = 1:numel(td)
    if ~isnan(td(i).idx_movement_on) && ~isnan(td(i).bumpDir) && td(i).idx_bumpTime < td(i).idx_goCueTime(1)
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
        if strcmp(td(i).result,'R') && ~isnan(td(i).idx_goCueTime(1))
            if isempty(td_temp)
                td_temp = td(i);
            else
                td_temp = [td_temp; td(i)];
            end
        end
    end
end
td_act = td_temp; %td non-bump




% Truncate bump trials
td_bumpSLR = td_bump;
for i = 1:numel(td_bump)
    idx = td_bump(i).idx_bumpTime:td_bump(i).idx_bumpTime+10;
    td_bumpSLR(i).pos = td_bump(i).pos(idx,:);
    td_bumpSLR(i).vel = td_bump(i).vel(idx,:);
    td_bumpSLR(i).muscle_len = td_bump(i).muscle_len(idx,:);
    td_bumpSLR(i).muscle_vel = td_bump(i).muscle_vel(idx,:);
    td_bumpSLR(i).emgNorm = td_bump(i).emgNorm(idx,:);
    td_bumpSLR(i).joint_ang = td_bump(i).joint_ang(idx,:);
    td_bumpSLR(i).joint_vel = td_bump(i).joint_vel(idx,:);
    td_bumpSLR(i).speed = td_bump(i).speed(idx,:);
    td_bumpSLR(i).S1_spikes = td_bump(i).S1_spikes(idx,:);
end

td_bumpLLR = td_bump;
for i = 1:numel(td_bump)
    idx = td_bump(i).idx_bumpTime+10:td_bump(i).idx_bumpTime+25;
    td_bumpLLR(i).pos = td_bump(i).pos(idx,:);
    td_bumpLLR(i).vel = td_bump(i).vel(idx,:);
    td_bumpLLR(i).muscle_len = td_bump(i).muscle_len(idx,:);
    td_bumpLLR(i).muscle_vel = td_bump(i).muscle_vel(idx,:);
    td_bumpLLR(i).emgNorm = td_bump(i).emgNorm(idx,:);
    td_bumpLLR(i).joint_ang = td_bump(i).joint_ang(idx,:);
    td_bumpLLR(i).joint_vel = td_bump(i).joint_vel(idx,:);
    td_bumpLLR(i).speed = td_bump(i).speed(idx,:);
    td_bumpLLR(i).S1_spikes = td_bump(i).S1_spikes(idx,:);
end


td_bumpS1 = td_bump;
for i = 1:numel(td_bump)
    idx = td_bump(i).idx_bumpTime:td_bump(i).idx_bumpTime+60;
    td_bumpS1(i).pos = td_bump(i).pos(idx,:);
    td_bumpS1(i).vel = td_bump(i).vel(idx,:);
    td_bumpS1(i).muscle_len = td_bump(i).muscle_len(idx,:);
    td_bumpS1(i).muscle_vel = td_bump(i).muscle_vel(idx,:);
    td_bumpS1(i).emgNorm = td_bump(i).emgNorm(idx,:);
    td_bumpS1(i).joint_ang = td_bump(i).joint_ang(idx,:);
    td_bumpS1(i).joint_vel = td_bump(i).joint_vel(idx,:);
    td_bumpS1(i).speed = td_bump(i).speed(idx,:);
    td_bumpS1(i).S1_spikes = td_bump(i).S1_spikes(idx,:);
end
%% Analyze bump EMG PDs
clearvars -except td_bumpSLR td_bumpLLR td_bump

params.out_signals = 'emgNorm';
params.in_signals = 'vel';
params.num_boots = 10;

workspaces = vertcat(td_bump.spaceNum);
td_1 = td_bump(workspaces==1);
td_2 = td_bump(workspaces==2);

emg_pdTable_1 = getTDPDs(td_1,params);
emg_pdTable_2 = getTDPDs(td_2,params);

PDs1 = rad2deg(emg_pdTable_1.velPD);
PDs2 = rad2deg(emg_pdTable_2.velPD);

MDs1 = rad2deg(emg_pdTable_1.velModdepth);
MDs2 = rad2deg(emg_pdTable_2.velModdepth);



plot(PDs1,PDs2,'.')

% figure; plot(MDs1,MDs2,'.')



%% Plot EMGs

for reflex = 1:2
    clearvars -except td_bumpSLR td_bumpLLR td_bumpS1 reflex

    bumpDirs = 0:pi/4:7*pi/4;

    if reflex == 1
        td_bump = td_bumpSLR;
        color1 = [0.9 0.2 0.2];
        color2 = [0.2 0.2 0.9];
    elseif reflex == 2
        td_bump = td_bumpLLR;
        color1 = [0.7 0.2 0.2];
        color2 = [0.2 0.2 0.7];

    end
    
    
    workspaces = vertcat(td_bump.spaceNum);
    td_1 = td_bump(workspaces==1);
    td_2 = td_bump(workspaces==2);
    
    nMus = 22;
    for mus = 1:nMus
        for i = 1:numel(td_1)
            meanEMG_1(mus,i) = mean(td_1(i).emgNorm(:,mus));
            bumpDir_1(mus,i) = td_1(i).bumpDir;
        end
        [bumpDir_1(mus,:),I1] = sort(bumpDir_1(mus,:));
        meanEMG_1(mus,:) = meanEMG_1(I1);
        smoothEMG1(mus,:) = smooth(meanEMG_1(mus,:),10);
        smoothSTD1(mus,:) = movstd(meanEMG_1(mus,:),10);
    end
    
    
    
    for mus = 1:nMus
        for i = 1:numel(td_2)
            meanEMG_2(mus,i) = mean(td_2(i).emgNorm(:,mus));
            bumpDir_2(mus,i) = td_2(i).bumpDir;
        end
        [bumpDir_2(mus,:),I2] = sort(bumpDir_2(mus,:));
        meanEMG_2(mus,:) = meanEMG_2(I2);
        smoothFRs2(mus,:) = smooth(meanEMG_2(mus,:),10);
        smoothSTD2(mus,:) = movstd(meanEMG_2(mus,:),10);
    end
    
    smoothEMG1(:,end+1) = smoothEMG1(:,1);
    smoothSTD1(:,end+1) = smoothSTD1(:,1);
    bumpDir_1(:,end+1) = bumpDir_1(:,1);
    smoothFRs2(:,end+1) = smoothFRs2(:,1);
    smoothSTD2(:,end+1) = smoothSTD2(:,1);
    bumpDir_2(:,end+1) = bumpDir_2(:,2);
    
    
    
    for mus = 1:nMus
        figure(mus); set(gcf,'Color','white'); 
        if reflex == 1
            a(mus) = polaraxes;
            set(a(mus),'FontName','Helvetica','FontSize',9,'RLim',[0 0.2]); hold on
        end
        %    polar(bumpDir_1(mus,:),meanEMG_1(mus,:),'r.'); hold on;
        h1(mus) = polarplot(bumpDir_1(mus,:),smoothEMG1(mus,:),'color',color1); hold on;
        h1(mus).LineWidth = 2;
        hs1_pos = polarplot(bumpDir_1(mus,:),smoothEMG1(mus,:) + smoothSTD1(mus,:)/10,'color',color1);
        hs1_neg = polarplot(bumpDir_1(mus,:),smoothEMG1(mus,:) - smoothSTD1(mus,:)/10,'color',color1);
        
        
        %    polar(bumpDir_2(mus,:),meanEMG_2(mus,:),'b.');
        h2(mus) = polarplot(bumpDir_2(mus,:),smoothFRs2(mus,:),'color',color2);
        h2(mus).LineWidth = 2;
        hs2_pos = polarplot(bumpDir_2(mus,:),smoothFRs2(mus,:) + smoothSTD2(mus,:)/10,'color',color2);
        hs2_neg = polarplot(bumpDir_2(mus,:),smoothFRs2(mus,:) - smoothSTD2(mus,:)/10,'color',color2);
        
    end
    
end


%% Plot S1 FRs

for reflex = 1
    clearvars -except td_bumpSLR td_bumpLLR td_bumpS1 reflex

    bumpDirs = 0:pi/4:7*pi/4;

    if reflex == 1
        td_bump = td_bumpS1;
        color1 = [0.9 0.2 0.2];
        color2 = [0.2 0.2 0.9];
    elseif reflex == 2
        td_bump = td_bumpLLR;
        color1 = [0.7 0.2 0.2];
        color2 = [0.2 0.2 0.7];

    end
    
    
    workspaces = vertcat(td_bump.spaceNum);
    td_1 = td_bump(workspaces==1);
    td_2 = td_bump(workspaces==2);
    
    nCell = size(td_bumpS1(1).S1_unit_guide,1);
    for cell = 1:nCell
        for i = 1:numel(td_1)
            meanFRs_1(cell,i) = mean(td_1(i).S1_spikes(:,cell));
            bumpDir_1(cell,i) = td_1(i).bumpDir;
        end
        [bumpDir_1(cell,:),I1] = sort(bumpDir_1(cell,:));
        meanFRs_1(cell,:) = meanFRs_1(I1);
        smoothFRs1(cell,:) = smooth(meanFRs_1(cell,:),30);
        smoothSTD1(cell,:) = movstd(meanFRs_1(cell,:),30);
    end
    
    
    
    for cell = 1:nCell
        for i = 1:numel(td_2)
            meanFRs_2(cell,i) = mean(td_2(i).S1_spikes(:,cell));
            bumpDir_2(cell,i) = td_2(i).bumpDir;
        end
        [bumpDir_2(cell,:),I2] = sort(bumpDir_2(cell,:));
        meanFRs_2(cell,:) = meanFRs_2(I2);
        smoothFRs2(cell,:) = smooth(meanFRs_2(cell,:),30);
        smoothSTD2(cell,:) = movstd(meanFRs_2(cell,:),30);
    end
    
    smoothFRs1(:,end+1) = smoothFRs1(:,1);
    smoothSTD1(:,end+1) = smoothSTD1(:,1);
    bumpDir_1(:,end+1) = bumpDir_1(:,1);
    smoothFRs2(:,end+1) = smoothFRs2(:,1);
    smoothSTD2(:,end+1) = smoothSTD2(:,1);
    bumpDir_2(:,end+1) = bumpDir_2(:,2);
    
    
    
    for cell = 1:nCell
        figure; set(gcf,'Color','white'); 
        if reflex == 1
            a(cell) = polaraxes;
            set(a(cell),'FontName','Helvetica','FontSize',9,'RLim',[0 50]); hold on
        end
        %    polar(bumpDir_1(cell,:),meanFRs_1(cell,:),'r.'); hold on;
        h1(cell) = polarplot(bumpDir_1(cell,:),smoothFRs1(cell,:),'color',color1); hold on;
        h1(cell).LineWidth = 2;
        hs1_pos = polarplot(bumpDir_1(cell,:),smoothFRs1(cell,:) + smoothSTD1(cell,:)/10,'color',color1);
        hs1_neg = polarplot(bumpDir_1(cell,:),smoothFRs1(cell,:) - smoothSTD1(cell,:)/10,'color',color1);
        
        
        %    polar(bumpDir_2(cell,:),meanEMG_2(cell,:),'b.');
        h2(cell) = polarplot(bumpDir_2(cell,:),smoothFRs2(cell,:),'color',color2);
        h2(cell).LineWidth = 2;
        hs2_pos = polarplot(bumpDir_2(cell,:),smoothFRs2(cell,:) + smoothSTD2(cell,:)/10,'color',color2);
        hs2_neg = polarplot(bumpDir_2(cell,:),smoothFRs2(cell,:) - smoothSTD2(cell,:)/10,'color',color2);
        
    end
    
end
