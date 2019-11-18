%% load & process TD file
clear; clc; 

pathname = '~/LimbLab/Projects/S1-gamma/data';
% filename = 'Duncan_20190524_CObumpmove_10ms.mat';
filename = 'Duncan_20190710_CObumpmove_10ms.mat';
% filename = 'Han_20191010_CObumpmove_10ms.mat';
% filename = 'Han_20191101_CObumpmove_10ms.mat';
% filename = 'Duncan_20191106_CObumpmove_10ms.mat';

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

%Get rectified velocity
td.vel_rect = [td.vel(:,1) td.vel(:,2) -td.vel(:,1) -td.vel(:,2)];
td.vel_rect(td.vel_rect < 0) = 0;

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
splitParams.linked_fields = {'result','bumpDir','tgtDir'};
td = splitTD(td,splitParams);

%Get movement onset
td(isnan([td.idx_goCueTime])) = [];

moveOnsetParams.start_idx = 'idx_goCueTime';
moveOnsetParams.end_idx = 'idx_endTime';
td = getMoveOnsetAndPeak(td,moveOnsetParams);

% Separate TD into bump, act, and comb structures

% Get bump-move trials
td_temp = [];
for i = 1:numel(td)
    if ~isnan(td(i).idx_movement_on) && ~isnan(td(i).bumpDir) && td(i).idx_bumpTime > td(i).idx_goCueTime
        if strcmp(td(i).result,'R') || strcmp(td(i).result,'F')
            if isempty(td_temp)
                td_temp = td(i);
            else
                td_temp = [td_temp; td(i)];
            end
        end
    end
end

td_comb = td_temp;

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



td_bl = td;


%% Plot move-bump trajectories
tgtDirs = 0:45:315;
relBumpDirs = 180;
numTgt = numel(tgtDirs);
numBump = numel(relBumpDirs);
numCond = numTgt*numBump;


kinfig = figure;
set(gcf,'Color','White')
i = 0;
j = 0;
for rd = tgtDirs
    i = i+1;
    j = 0;
    for rbd = relBumpDirs
        j = j+1;
        for trial = 1:numel(td_comb)
            plotStartIdx = td_comb(trial).idx_bumpTime;
            plotStopIdx = td_comb(trial).idx_bumpTime + 13;
            absBumpDir = td_comb(trial).bumpDir;
            relBumpDir = absBumpDir - td_comb(trial).tgtDir;
            if relBumpDir < 0
                relBumpDir = relBumpDir + 360;
            end
            if td_comb(trial).tgtDir == rd && relBumpDir == rbd
                set(0,'CurrentFigure',kinfig)
                set(gcf,'visible','on')
                subplot(numTgt,numBump,(i-1)*(numBump)+j); hold on;
                axis([-10 10 -10 10])
                [targX,targY] = pol2cart(deg2rad(rd),9);
                plot(targX,targY,'ko','MarkerFaceColor','k','markersize',5);
                % Get movement directions and colors set up
                velX = td_comb(trial).vel(plotStartIdx:plotStopIdx,1);
                velY = td_comb(trial).vel(plotStartIdx:plotStopIdx,2);
                traj = atan2(velY,velX);
                traj = rad2deg(traj);
                trajColorMap = [0.9 0.25 0; 0.3 0.2 1; 0.9 0.25 0];
                x = linspace(-180,179.9999,3);
                xq = td_comb(trial).tgtDir - traj;
                xq(xq<-180) = xq(xq<-180) + 360;
                xq(xq>180) = xq(xq>180) - 360;
                v = trajColorMap;
                vq = interp1(x,v,xq);
                idxStart = td_comb(trial).idx_goCueTime;
                plot(td_comb(trial).pos(idxStart:plotStartIdx,1),td_comb(trial).pos(idxStart:plotStartIdx,2),'.','markersize',6,'color',[0.9 0.9 0.9])
                for t = 1:numel(traj)
                    plot(td_comb(trial).pos(plotStartIdx+t,1),td_comb(trial).pos(plotStartIdx+t,2),'.','markersize',6,'color',vq(t,:)), hold on
                end
                % Formating
                set(gca,'TickDir','out','clipping','off')
                if j ~= 1
                    set(gca,'yticklabel',[])
                else
                    ylabel(num2str(rd))
                end
                if i ~= numTgt
                    set(gca,'xticklabel',[])
                end
                if i == 1
                    title(num2str(rbd))
                end
            end
        end
    end
end



%% Plot bump only trajectories

kinfig = figure;
i = 1;
j = 0;

tgtDirs = 0:45:315;
relBumpDirs = 0:45:315;
numTgt = numel(tgtDirs);
numBump = numel(relBumpDirs);
numCond = numTgt*numBump;

for rbd = relBumpDirs
    j = j+1;
    for trial = 1:numel(td_bump)
        plotStartIdx = td_bump(trial).idx_bumpTime;
        plotStopIdx = td_bump(trial).idx_bumpTime+13;
        absBumpDir = td_bump(trial).bumpDir;
        relBumpDir = absBumpDir;
        
        if relBumpDir < 0
            relBumpDir = relBumpDir + 360;
        end
        if relBumpDir >= 360
            relBumpDir = relBumpDir - 360;
        end
        
        if relBumpDir == rbd
            set(0,'CurrentFigure',kinfig)
            set(gcf,'visible','on')
            subplot(numBump,numTgt,(i-1)*(numBump)+j); hold on;
            axis([-3 3 -3 3])
            % Get movement directions and colors set up
            velX = td_bump(trial).vel(plotStartIdx:plotStopIdx,1);
            velY = td_bump(trial).vel(plotStartIdx:plotStopIdx,2);
            traj = atan2(velY,velX);
            traj = rad2deg(traj);
            trajColorMap = [1 0 0; 0.0 0.5 1; 1 0 0];
            x = linspace(-180,179.9999,3);
            xq = relBumpDir - traj;
            xq(xq<-180) = xq(xq<-180) + 360;
            xq(xq>180) = xq(xq>180) - 360;
            v = trajColorMap;
            vq = interp1(x,v,xq);
            idxStart = td_bump(trial).idx_goCueTime;
            plot(td_bump(trial).pos(idxStart:plotStartIdx,1),td_bump(trial).pos(idxStart:plotStartIdx,2),'.','markersize',6,'color',[0.9 0.9 0.9])
            for t = 1:numel(traj)
                plot(td_bump(trial).pos(plotStartIdx+t,1),td_bump(trial).pos(plotStartIdx+t,2),'.','markersize',6,'color',vq(t,:)), hold on
            end
            % Formating
            set(gca,'TickDir','out','clipping','off')
            if j ~= 1
                set(gca,'yticklabel',[])
            else
                ylabel(num2str(rbd))
            end
            if i ~= numTgt
                set(gca,'xticklabel',[])
            end
            if i == 1
                title(num2str(rbd))
            end
        end
    end
end


%% Plot active trajectories
kinfig = figure;
i = 0;
j = 1;
tgtDirs = 0:45:315;


for rd = tgtDirs
    i = i+1;
    for trial = 1:numel(td_act)
        plotStartIdx = td_act(trial).idx_movement_on;
        plotStopIdx = td_act(trial).idx_movement_on+13;
        tgtDir = td_act(trial).tgtDir;
        
        
        
        
        if tgtDir == rd
            set(0,'CurrentFigure',kinfig)
            set(gcf,'visible','on')
            subplot(numTgt,numBump,(i-1)*(numBump)+j); hold on;
            axis([-10 10 -10 10])
            [targX,targY] = pol2cart(deg2rad(rd),9);
            plot(targX,targY,'ko','MarkerFaceColor','k','markersize',5);
            % Get movement directions and colors set up
            velX = td_act(trial).vel(plotStartIdx:plotStopIdx,1);
            velY = td_act(trial).vel(plotStartIdx:plotStopIdx,2);
            traj = atan2(velY,velX);
            traj = rad2deg(traj);
            trajColorMap = [1 0 0; 0.0 0.5 1; 1 0 0];
            x = linspace(-180,179.9999,3);
            xq = td_act(trial).tgtDir - traj;
            xq(xq<-180) = xq(xq<-180) + 360;
            xq(xq>180) = xq(xq>180) - 360;
            v = trajColorMap;
            vq = interp1(x,v,xq);
            idxStart = td_act(trial).idx_goCueTime;
            plot(td_act(trial).pos(idxStart:plotStartIdx,1),td_act(trial).pos(idxStart:plotStartIdx,2),'.','markersize',6,'color',[0.9 0.9 0.9])
            for t = 1:numel(traj)
                plot(td_act(trial).pos(plotStartIdx+t,1),td_act(trial).pos(plotStartIdx+t,2),'.','markersize',6,'color',vq(t,:)), hold on
            end
            % Formating
            set(gca,'TickDir','out','clipping','off')
            if j ~= 1
                set(gca,'yticklabel',[])
            else
                ylabel(num2str(rd))
            end
            if i ~= numTgt
                set(gca,'xticklabel',[])
            end
            if i == 1
                title(num2str(rd))
            end
        end
    end
end

%% Prep TDs for encoding/decoding models
modelNames = {'act','pas','comb'};

srParams.powers = [0.5];
srParams.signals = {'vel'};

td_act_prep = td_act;

for i = 1:numel(td_act)
    idx = td_act(i).idx_movement_on:td_act(i).idx_movement_on+13;
    td_act(i).pos = td_act(i).pos(idx,:);
    td_act(i).vel = td_act(i).vel(idx,:);
    td_act(i).vel_rect = td_act(i).vel_rect(idx,:);
    td_act(i).acc = td_act(i).acc(idx,:);
    td_act(i).speed = td_act(i).speed(idx,:);
    td_act(i).force = td_act(i).force(idx,:);
    td_act(i).S1_spikes = td_act(i).S1_spikes(idx,:);
    td_act(i).meanTraj = rad2deg(atan2(mean(td_act(i).vel(:,2)),mean(td_act(i).vel(:,1))));
    td_act(i).forceTraj = rad2deg(atan2(mean(td_act(i).force(:,2)),mean(td_act(i).force(:,1))));
    td_act(i).accTraj = rad2deg(atan2(mean(td_act(i).acc(:,2)),mean(td_act(i).acc(:,1))));
end

for i = 1:numel(td_act_prep)
    idx = td_act_prep(i).idx_goCueTime-13:td_act_prep(i).idx_goCueTime;
    td_act_prep(i).pos = td_act_prep(i).pos(idx,:);
    td_act_prep(i).vel = td_act_prep(i).vel(idx,:);
    td_act_prep(i).vel_rect = td_act_prep(i).vel_rect(idx,:);
    td_act_prep(i).acc = td_act_prep(i).acc(idx,:);
    td_act_prep(i).speed = td_act_prep(i).speed(idx,:);
    td_act_prep(i).force = td_act_prep(i).force(idx,:);
    td_act_prep(i).S1_spikes = td_act_prep(i).S1_spikes(idx,:);
    td_act_prep(i).meanTraj = rad2deg(atan2(mean(td_act_prep(i).vel(:,2)),mean(td_act_prep(i).vel(:,1))));
    td_act_prep(i).forceTraj = rad2deg(atan2(mean(td_act_prep(i).force(:,2)),mean(td_act_prep(i).force(:,1))));
    td_act_prep(i).accTraj = rad2deg(atan2(mean(td_act_prep(i).acc(:,2)),mean(td_act_prep(i).acc(:,1))));
end


for i = 1:numel(td_bump)
    idx = td_bump(i).idx_bumpTime:td_bump(i).idx_bumpTime+13;
    td_bump(i).pos = td_bump(i).pos(idx,:);
    td_bump(i).vel = td_bump(i).vel(idx,:);
    td_bump(i).vel_rect = td_bump(i).vel_rect(idx,:);
    td_bump(i).acc = td_bump(i).acc(idx,:);
    td_bump(i).speed = td_bump(i).speed(idx,:);
    td_bump(i).force = td_bump(i).force(idx,:);
    td_bump(i).S1_spikes = td_bump(i).S1_spikes(idx,:);
    td_bump(i).meanTraj = rad2deg(atan2(mean(td_bump(i).vel(:,2)),mean(td_bump(i).vel(:,1))));
    td_bump(i).forceTraj = rad2deg(atan2(mean(td_bump(i).force(:,2)),mean(td_bump(i).force(:,1))));
    td_bump(i).accTraj = rad2deg(atan2(mean(td_bump(i).acc(:,2)),mean(td_bump(i).acc(:,1))));
end

for i = 1:numel(td_comb)
    idx = td_comb(i).idx_bumpTime:td_comb(i).idx_bumpTime+13;
    td_comb(i).pos = td_comb(i).pos(idx,:);
    td_comb(i).vel = td_comb(i).vel(idx,:);
    td_comb(i).vel_rect = td_comb(i).vel_rect(idx,:);
    td_comb(i).acc = td_comb(i).acc(idx,:);
    td_comb(i).speed = td_comb(i).speed(idx,:);
    td_comb(i).force = td_comb(i).force(idx,:);
    td_comb(i).S1_spikes = td_comb(i).S1_spikes(idx,:);
    td_comb(i).S1_spikes_bins = td_comb(i).S1_spikes_bins(idx,:);
    td_comb(i).meanTraj = rad2deg(atan2(mean(td_comb(i).vel(:,2)),mean(td_comb(i).vel(:,1))));
    td_comb(i).forceTraj = rad2deg(atan2(mean(td_comb(i).force(:,2)),mean(td_comb(i).force(:,1))));
    td_comb(i).accTraj = rad2deg(atan2(mean(td_comb(i).acc(:,2)),mean(td_comb(i).acc(:,1))));
    td_comb(i).trajDiff = abs(td_comb(i).tgtDir - td_comb(i).meanTraj);
end

for i = 1:numel(td_bl)
    idx = td_bl(i).idx_startTime+10:td_bl(i).idx_startTime+23;
    td_bl(i).pos = td_bl(i).pos(idx,:);
    td_bl(i).vel = td_bl(i).vel(idx,:);
    td_bl(i).acc = td_bl(i).acc(idx,:);
    td_bl(i).speed = td_bl(i).speed(idx,:);
    td_bl(i).force = td_bl(i).force(idx,:);
    td_bl(i).S1_spikes = td_bl(i).S1_spikes(idx,:);
    td_bl(i).S1_spikes_bins = td_bl(i).S1_spikes_bins(idx,:);
    td_bl(i).meanTraj = rad2deg(atan2(mean(td_bl(i).vel(:,2)),mean(td_bl(i).vel(:,1))));
    td_bl(i).forceTraj = rad2deg(atan2(mean(td_bl(i).force(:,2)),mean(td_bl(i).force(:,1))));
    td_bl(i).accTraj = rad2deg(atan2(mean(td_bl(i).acc(:,2)),mean(td_bl(i).acc(:,1))));
    td_bl(i).trajDiff = abs(td_bl(i).tgtDir - td_bl(i).meanTraj);
end
%% PD Calculations

params.out_signals = 'S1_spikes';
params.in_signals = 'vel';
params.num_boots = 10;

tgtDirs_comb = vertcat(td_comb.tgtDir);
bumpDirs_comb = vertcat(td_comb.bumpDir);
totalDirs_comb = bumpDirs_comb - tgtDirs_comb;


while sum(totalDirs_comb >= 360) > 0
    totalDirs_comb(totalDirs_comb>=360) = totalDirs_comb(totalDirs_comb>=360) - 360;
end

% td_assist = td_comb(totalDirs_comb==0);
td_resist = td_comb(totalDirs_comb==180);
% td_perpen = td_comb(totalDirs_comb==90 | totalDirs_comb==270);


% td_resist = td_resist_late;

act_pdTable = getTDPDs(td_act,params);
bump_pdTable = getTDPDs(td_bump,params);
comb_pdTable = getTDPDs(td_comb,params);
prep_pdTable = getTDPDs(td_act_prep,params);
res_pdTable = getTDPDs(td_resist,params);
% ass_pdTable = getTDPDs(td_assist,params);

if strfind(params.in_signals,'vel')
    actPDs = rad2deg(act_pdTable.velPD);
    actPDs(actPDs<0) = actPDs(actPDs<0)+360;
    actMDs = act_pdTable.velModdepth;
    
    pasPDs = rad2deg(bump_pdTable.velPD);
    pasPDs(pasPDs<0) = pasPDs(pasPDs<0)+360;
    pasMDs = bump_pdTable.velModdepth;
    
    combPDs = rad2deg(comb_pdTable.velPD);
    combPDs(combPDs<0) = combPDs(combPDs<0)+360;
    combMDs = comb_pdTable.velModdepth;
    
    prepPDs = rad2deg(prep_pdTable.velPD);
    prepPDs = wrapTo360(prepPDs);
    prepMDs = prep_pdTable.velModdepth;
    
    resPDs = rad2deg(res_pdTable.velPD);
    resPDs(resPDs<0) = resPDs(resPDs<0)+360;
    resMDs = res_pdTable.velModdepth;
%     
%     assPDs = rad2deg(ass_pdTable.velPD);
%     assPDs(assPDs<0) = assPDs(assPDs<0)+360;
%     assMDs = ass_pdTable.velModdepth;
    
elseif strfind(params.in_signals,'acc')
    actPDs = rad2deg(act_pdTable.accPD);
    actPDs(actPDs<0) = actPDs(actPDs<0)+360;
    actMDs = act_pdTable.accModdepth;
    
    pasPDs = rad2deg(bump_pdTable.accPD);
    pasPDs(pasPDs<0) = pasPDs(pasPDs<0)+360;
    pasMDs = bump_pdTable.accModdepth;
    
    combPDs = rad2deg(comb_pdTable.accPD);
    combPDs(combPDs<0) = combPDs(combPDs<0)+360;
    combMDs = comb_pdTable.accModdepth;
    
    resPDs = rad2deg(res_pdTable.accPD);
    resPDs(resPDs<0) = resPDs(resPDs<0)+360;
    resMDs = res_pdTable.accModdepth;
%     
%     assPDs = rad2deg(ass_pdTable.accPD);
%     assPDs(assPDs<0) = assPDs(assPDs<0)+360;
%     assMDs = ass_pdTable.accModdepth;
end



%% Make figures a la London and Miller 2013
actCI = rad2deg(wrapTo2Pi(act_pdTable.velPDCI));
pasCI = rad2deg(wrapTo2Pi(bump_pdTable.velPDCI));
combCI = rad2deg(wrapTo2Pi(comb_pdTable.velPDCI));
prepCI = rad2deg(wrapTo2Pi(prep_pdTable.velPDCI));
resCI = rad2deg(wrapTo2Pi(res_pdTable.velPDCI));
% assCI = rad2deg(wrapTo2Pi(ass_pdTable.velPDCI));

%Pas2act PD scatter
figure; hold on; set(gcf,'Color','White','Units','Normalized','Position',[0.25 0.25 0.5 0.5])
subplot(3,4,1), hold on;
plot(pasPDs,actPDs,'k.','MarkerSize',10);
errorbar(pasPDs,actPDs,actCI(:,1)-actPDs,actCI(:,2)-actPDs,pasCI(:,1)-pasPDs,pasCI(:,2)-pasPDs,...
    'LineStyle','none','Color','k');
set(gca,'FontName','Helvetica','FontSize',12), axis([0 360 0 360]) 
xlabel('Center-hold bump PDs'), ylabel('Active reach PDs')

shadeX = rad2deg([0 2*pi pi 0]);
shadeY = rad2deg([0 2*pi 2*pi pi]);
fill(shadeX,shadeY,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
shadeX2 = rad2deg([pi 2*pi 2*pi pi]);
shadeY2 = rad2deg([0 0 pi 0]);
fill(shadeX2,shadeY2,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)

%Pas2act differences
p2aPDdiff = pasPDs - actPDs;
p2aPDdiff = wrapTo180(p2aPDdiff);
subplot(3,4,5), hold on;
edges = -0:5:180;
histogram(abs(p2aPDdiff),edges,'FaceColor','k')
xlabel('PD difference (pas - act)'), ylabel('Number of neurons')
set(gca,'FontName','Helvetica','FontSize',12),axis([0 180 0 20])

subplot(3,4,9), hold on;
plot(pasMDs,actMDs,'.'); 
line([0;0.5],[0.0;0.5],'LineWidth',2,'Color','k'); 
axis([0 0.3 0 0.3])



%pas2res PD scatter
subplot(3,4,2), hold on;
plot(pasPDs,resPDs,'k.','MarkerSize',10);
errorbar(pasPDs,resPDs,resCI(:,1)-resPDs,resCI(:,2)-resPDs,pasCI(:,1)-pasPDs,pasCI(:,2)-pasPDs,...
    'LineStyle','none','Color','k');
set(gca,'FontName','Helvetica','FontSize',12), axis([0 360 0 360]) 
xlabel('Center-hold bump PDs'),ylabel('Perturbations during reach PDs')

shadeX = rad2deg([0 2*pi pi 0]);
shadeY = rad2deg([0 2*pi 2*pi pi]);
fill(shadeX,shadeY,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
shadeX2 = rad2deg([pi 2*pi 2*pi pi]);
shadeY2 = rad2deg([0 0 pi 0]);
fill(shadeX2,shadeY2,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)

%Pas2res differences
p2cPDdiff = pasPDs - resPDs;
p2cPDdiff = wrapTo180(p2cPDdiff);
subplot(3,4,6), hold on;
histogram(abs(p2cPDdiff),edges,'FaceColor','k')
xlabel('PD difference (pas - res)'), ylabel('Number of neurons')
set(gca,'FontName','Helvetica','FontSize',12),axis([0 180 0 20])

subplot(3,4,10), hold on;
plot(pasMDs,resMDs,'.'); 
line([0;0.5],[0.0;0.5],'LineWidth',2,'Color','k'); 
axis([0 0.3 0 0.3])



%res2act PD scatter
subplot(3,4,3), hold on;
res2act_idx = abs(resCI(:,2) - resCI(:,1))<=45 & abs(actCI(:,2) - actCI(:,1))<=45;

plot(resPDs(res2act_idx),actPDs(res2act_idx),'k.','MarkerSize',10);
errorbar(resPDs(res2act_idx),actPDs(res2act_idx),...
    actCI(res2act_idx,1)-actPDs(res2act_idx),...
    actCI(res2act_idx,2)-actPDs(res2act_idx),...
    resCI(res2act_idx,1)-resPDs(res2act_idx),...
    resCI(res2act_idx,2)-resPDs(res2act_idx),...
    'LineStyle','none','Color','k');
set(gca,'FontName','Helvetica','FontSize',12), axis([0 360 0 360]) 
xlabel('Perturbations during reach PDs'), ylabel('Active reach PDs')

shadeX = rad2deg([0 2*pi pi 0]);
shadeY = rad2deg([0 2*pi 2*pi pi]);
fill(shadeX,shadeY,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
shadeX2 = rad2deg([pi 2*pi 2*pi pi]);
shadeY2 = rad2deg([0 0 pi 0]);
fill(shadeX2,shadeY2,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)

%res2act differences
r2aPDdiff = resPDs(res2act_idx) - actPDs(res2act_idx);
r2aPDdiff = wrapTo180(r2aPDdiff);
subplot(3,4,7), hold on;
histogram(abs(r2aPDdiff),edges,'FaceColor','k')
xlabel('PD difference (res - act)'), ylabel('Number of neurons')
set(gca,'FontName','Helvetica','FontSize',12),axis([0 180 0 20])

subplot(3,4,11), hold on;
plot(resMDs(res2act_idx),actMDs(res2act_idx),'.'); 
line([0;0.5],[0.0;0.5],'LineWidth',2,'Color','k'); 
axis([0 0.3 0 0.3])

%prep2act PD scatter
subplot(3,4,4), hold on;
prep2act_idx = abs(prepCI(:,2) - prepCI(:,1))<=45 & abs(actCI(:,2) - actCI(:,1))<=45;

plot(actPDs(prep2act_idx),prepPDs(prep2act_idx),'k.','MarkerSize',10);
errorbar(actPDs(prep2act_idx),prepPDs(prep2act_idx),...
    prepCI((prep2act_idx),1)-prepPDs(prep2act_idx),...
    prepCI((prep2act_idx),2)-prepPDs(prep2act_idx),...
    actCI((prep2act_idx),1)-actPDs(prep2act_idx),...
    actCI((prep2act_idx),2)-actPDs(prep2act_idx),...
    'LineStyle','none','Color','k');
set(gca,'FontName','Helvetica','FontSize',12), axis([0 360 0 360]) 
xlabel('Preparatory'), ylabel('Active reach PDs')

shadeX = rad2deg([0 2*pi pi 0]);
shadeY = rad2deg([0 2*pi 2*pi pi]);
fill(shadeX,shadeY,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
shadeX2 = rad2deg([pi 2*pi 2*pi pi]);
shadeY2 = rad2deg([0 0 pi 0]);
fill(shadeX2,shadeY2,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)

%prep2act differences
r2aPDdiff = prepPDs(prep2act_idx) - actPDs(prep2act_idx);
r2aPDdiff(r2aPDdiff > 180) = r2aPDdiff(r2aPDdiff > 180) - 180;
r2aPDdiff(r2aPDdiff < -180) = r2aPDdiff(r2aPDdiff < -180) + 180;
subplot(3,4,8), hold on;
histogram(abs(r2aPDdiff),edges,'FaceColor','k')
xlabel('PD difference (prep - act)'), ylabel('Number of neurons')
set(gca,'FontName','Helvetica','FontSize',12),axis([0 180 0 20])

subplot(3,4,12), hold on;
plot(prepMDs(prep2act_idx),actMDs(prep2act_idx),'.'); 
line([0;0.5],[0.0;0.5],'LineWidth',2,'Color','k'); 
axis([0 0.3 0 0.3])



%% Decoders
% Train GLM decoder on comb data, try to decode act/pas data

tgtDirs_comb = vertcat(td_comb.tgtDir);
bumpDirs_comb = vertcat(td_comb.bumpDir);
totalDirs_comb = bumpDirs_comb - tgtDirs_comb;


while sum(totalDirs_comb >= 360) > 0
    totalDirs_comb(totalDirs_comb>=360) = totalDirs_comb(totalDirs_comb>=360) - 360;
end

td_assist = td_comb(totalDirs_comb==0);
td_resist = td_comb(totalDirs_comb==180);
td_perpen = td_comb(totalDirs_comb==90 | totalDirs_comb==270);


td_resist = [td_resist; td_perpen];

clear modelParams;
modelParams.in_signals = {'S1_spikes','all'};
modelParams.model_type = 'linmodel';
modelParams.eval_metric = 'vaf';

%%% Kinematic decoders %%%
modelParams.out_signals = {'acc',1:2;'vel','all';'force',1:2};


% Decoder trained on movement bumps
modelParams.model_name = 'kin_decoder_comb';
[td_comb, decoder_comb_pos] = getModel(td_comb,modelParams); %Get model trained on comb data
test_comb2comb_decoder_pos = evalModel(td_comb,modelParams); %Test model trained on comb data on comb data

td_bump = getModel(td_bump,decoder_comb_pos);
test_comb2pas_decoder_pos = evalModel(td_bump,modelParams);%Test model trained on comb data on pas data

td_act = getModel(td_act,decoder_comb_pos);
test_comb2act_decoder_pos = evalModel(td_act,modelParams);%Test model trained on comb data on act data


% Decoder trained on active reaches
modelParams.model_name = 'kin_decoder_act';
modelParams.nFolds = 5;
td_act = td_xval(td_act,modelParams);
[~, decoder_act_pos] = getModel(td_act,modelParams); %Get model trained on act data
test_act2act_decoder_pos = evalModel(td_act,modelParams); %Test model trained on act data on act data

td_bump = getModel(td_bump,decoder_act_pos);
test_act2pas_decoder_pos = evalModel(td_bump,modelParams);%Test model trained on act data on pas data

td_comb = getModel(td_comb,decoder_act_pos);
test_act2comb_decoder_pos = evalModel(td_comb,modelParams);%Test model trained on act data on comb data

td_resist = getModel(td_resist,decoder_act_pos);
test_act2resist_decoder_pos = evalModel(td_resist,modelParams);%Test model trained on act data on comb data

td_assist = getModel(td_assist,decoder_act_pos);
test_act2assist_decoder_pos = evalModel(td_assist,modelParams);%Test model trained on act data on comb data


% Decoder trained on center-hold bumps
modelParams.model_name = 'kin_decoder_pas';
[td_bump, decoder_pas_pos] = getModel(td_bump,modelParams); %Get model trained on pas data
test_pas2pas_decoder_pos = evalModel(td_bump,modelParams); %Test model trained on pas data on pas data

td_act = getModel(td_act,decoder_pas_pos);
test_pas2act_decoder_pos = evalModel(td_act,modelParams);%Test model trained on pas data on act data

td_comb = getModel(td_comb,decoder_pas_pos);
test_pas2comb_decoder_pos = evalModel(td_comb,modelParams);%Test model trained on act data on comb data




%% Plot decoder results
% Active trials
for i = 1:numel(td_act)
    act2act_kinPred(i,:) = mean(td_act(i).linmodel_kin_decoder_act);
    pas2act_kinPred(i,:) = mean(td_act(i).linmodel_kin_decoder_pas);
    comb2act_kinPred(i,:) = mean(td_act(i).linmodel_kin_decoder_comb);
    vel_act(i,:) = mean(td_act(i).vel);
    FR_act(i,:) = mean(td_act(i).S1_spikes);
end

for j = 1:numel(td_bump)
    pas2pas_kinPred(j,:) = mean(td_bump(j).linmodel_kin_decoder_pas);
    act2pas_kinPred(j,:) = mean(td_bump(j).linmodel_kin_decoder_act);
    comb2pas_kinPred(j,:) = mean(td_bump(j).linmodel_kin_decoder_comb);
    vel_pas(j,:) = mean(td_bump(j).vel);
%     fracVel_pas(j,:) = mean(td_bump(j).vel_0_5);
    force_pas(j,:) = mean(td_bump(j).force(1:2));
    FR_pas(j,:) = mean(td_bump(j).S1_spikes);
end

for k = 1:numel(td_comb)
    vel_comb(k,:) = mean(td_comb(k).vel);
    comb2comb_kinPred(k,:) = mean(td_comb(k).linmodel_kin_decoder_comb);
    pas2comb_kinPred(k,:) = mean(td_comb(k).linmodel_kin_decoder_pas);
    act2comb_kinPred(k,:) = mean(td_comb(k).linmodel_kin_decoder_pas);
    FR_comb(k,:) = mean(td_comb(k).S1_spikes);
end

for i = 1:numel(td_resist)
act2resist_pred(i,:) = mean(td_resist(i).linmodel_kin_decoder_act);
vel_resist(i,:) = mean(td_resist(i).vel);
end

for i = 1:numel(td_assist)
act2assist_pred(i,:) = mean(td_assist(i).linmodel_kin_decoder_act);
vel_assist(i,:) = mean(td_assist(i).vel);
end

plot(vel_pas,act2pas_kinPred(:,3:4),'.'); hold on;
line([-30 30],[-30 30],'LineWidth',2,'Color','k'); line([-30 30],[30 -30],'LineWidth',2,'Color','k');
xlabel('hand velocity (cm/s)'); ylabel('train:act; test:pas; Velocity')
axis([-40 40 -40 40])

figure;
plot(vel_act,pas2act_kinPred(:,3:4),'.'); hold on;
line([-30 30],[-30 30],'LineWidth',2,'Color','k'); line([-30 30],[30 -30],'LineWidth',2,'Color','k');
xlabel('hand velocity (cm/s)'); ylabel('train:pas; test:act; Velocity')
axis([-40 40 -40 40])

figure;
plot(vel_act,comb2act_kinPred(:,3:4),'.'); hold on;
line([-30 30],[-30 30],'LineWidth',2,'Color','k'); line([-30 30],[30 -30],'LineWidth',2,'Color','k');
xlabel('hand velocity (cm/s)'); ylabel('train:comb; test:act; Velocity')
axis([-40 40 -40 40])

figure;
% plot(vel_comb,act2comb_kinPred(:,3:4),'k.'); hold on
plot(vel_resist,act2resist_pred(:,3:4),'b.'); hold on;
plot(vel_assist,act2assist_pred(:,3:4),'r.'); hold on
line([-30 30],[-30 30],'LineWidth',2,'Color','k'); line([-30 30],[30 -30],'LineWidth',2,'Color','k');
xlabel('hand velocity (cm/s)'); ylabel('train:act; test:comb; Velocity')
axis([-40 40 -40 40])

figure;
plot(vel_act,act2act_kinPred(:,3:4),'.'); hold on
line([-30 30],[-30 30],'LineWidth',2,'Color','k'); line([-30 30],[30 -30],'LineWidth',2,'Color','k');
xlabel('hand velocity (cm/s)'); ylabel('train:act; test:act; Velocity')
axis([-40 40 -40 40])







%% Plot decoder inputs and outputs
figure; hold on;
set(gcf,'Color','white')
trialIdx = 1:20;
bs = td_comb(1).bin_size;
rates = vertcat(td_comb(trialIdx).S1_spikes_bins);

ts = cell(size(rates,2),1);
for neuron = 1:size(rates,2)
    for bin = 1:size(rates,1)
        thisTime = (bin-1)*bs;
        numSpikes = rates(bin,neuron);
        ints = 1:numSpikes;
        spaces = ints*bs/((numSpikes-1)+2);
        if numel(spaces)>1
            newSpikes = spaces(2:end-1);
        else
            newSpikes = spaces;
        end
        newSpikes = newSpikes + 0.005*randn(size(newSpikes));
        ts{neuron} = [ts{neuron} newSpikes+thisTime];
    end
end


outputs = vertcat(td_act(trialIdx).linmodel_kin_decoder_act);
outputs = outputs(:,3:4);
vels = vertcat(td_act(trialIdx).vel);

subplot(2,1,1); hold on
set(gca,'FontName','Helvetica','tickdir','out')

for i = 1:numel(ts)
   plot(ts{i},i*ones(size(ts{i})),'.','MarkerSize',0.1), hold on  
end
axis([0 3.0 0 numel(ts)])


subplot(2,1,2); hold on
set(gca,'FontName','Helvetica','tickdir','out')
plot(outputs)
plot(vels)


vaf = test_act2act_decoder_pos(3:4);
text(0,0,['vaf = ' num2str(round(vaf*100)/100)])

%% Split bumps during movements decoders

tgtDirs_comb = vertcat(td_comb.tgtDir);
bumpDirs_comb = vertcat(td_comb.bumpDir);
totalDirs_comb = bumpDirs_comb - tgtDirs_comb;


while sum(totalDirs_comb >= 360) > 0
    totalDirs_comb(totalDirs_comb>=360) = totalDirs_comb(totalDirs_comb>=360) - 360;
end

td_assist = td_comb(totalDirs_comb==0);
td_resist = td_comb(totalDirs_comb==180);
td_perpen = td_comb(totalDirs_comb==90 | totalDirs_comb==270);

td_resist = [td_resist; td_perpen];

% Train linear decoder on 

clear modelParams;
modelParams.in_signals = {'S1_spikes','all'};
modelParams.model_type = 'linmodel';
modelParams.eval_metric = 'r2';

%%% Kinematic decoders %%%
modelParams.out_signals = {'acc',1:2;'vel','all';'force',1:2};


% Decoder trained on resistive movement bumps
modelParams.model_name = 'kin_decoder_resist';
modelParams.nFolds = 10;


td_resist = td_xval(td_resist,modelParams);
test_resist2resist = evalModel(td_resist,modelParams);

[TEST, decoder_resist] = getModel(td_resist,modelParams);
TESTEVAL = evalModel(TEST,modelParams);

td_assist = getModel(td_assist,decoder_resist);
test_resist2assist = evalModel(td_assist,modelParams);

td_perpen = getModel(td_perpen,decoder_resist);
test_resist2perpen = evalModel(td_perpen,modelParams);

td_act = getModel(td_act,decoder_resist);
test_resist2act = evalModel(td_act,modelParams);

td_bump = getModel(td_bump,decoder_resist);
test_resist2bump = evalModel(td_bump,modelParams);




figure;
for i = 1:numel(td_resist)
resist2resist_pred(i,:) = mean(td_resist(i).linmodel_kin_decoder_resist);
vel_resist(i,:) = mean(td_resist(i).vel);
end
plot(vel_resist,resist2resist_pred(:,3:4),'.'); hold on;
line([-30 30],[-30 30],'LineWidth',2,'Color','k'); line([-30 30],[30 -30],'LineWidth',2,'Color','k');
xlabel('hand velocity (cm/s)'); ylabel('train:resist; test:resist; Velocity')
axis([-40 40 -40 40])

figure;
for i = 1:numel(td_assist)
resist2assist_pred(i,:) = mean(td_assist(i).linmodel_kin_decoder_resist);
vel_assist(i,:) = mean(td_assist(i).vel);
end
plot(vel_assist,resist2assist_pred(:,3:4),'.'); hold on
line([-30 30],[-30 30],'LineWidth',2,'Color','k'); line([-30 30],[30 -30],'LineWidth',2,'Color','k');
xlabel('hand velocity (cm/s)'); ylabel('train:resist; test:assist; Velocity')
axis([-40 40 -40 40])


figure;
for i = 1:numel(td_assist)
resist2perpen_pred(i,:) = mean(td_perpen(i).linmodel_kin_decoder_resist);
vel_perpen(i,:) = mean(td_perpen(i).vel);
end
plot(vel_perpen,resist2perpen_pred(:,3:4),'.'); hold on
line([-30 30],[-30 30],'LineWidth',2,'Color','k'); line([-30 30],[30 -30],'LineWidth',2,'Color','k');
xlabel('hand velocity (cm/s)'); ylabel('train:resist; test:perpen; Velocity')
axis([-40 40 -40 40])


%% Separate act TDs into direction, plot FR as function of velocity
clearvars -except td_act td_bump td_comb td_resist

% td_comb = [];
% td_comb = td_resist;
% 
params.out_signals = 'S1_spikes';
params.in_signals = 'vel';
params.num_boots = 1;
actpdTable = getTDPDs(td_act,params);

actPDs = rad2deg(actpdTable.velPD);
actPDs(actPDs<0) = actPDs(actPDs<0)+360;

tgtDirs_act = vertcat(td_act.tgtDir);

td_act000 = td_act(tgtDirs_act==0);
for i = 1:numel(td_act000)
    meanSpeed000(i,:) = mean(td_act000(i).speed);
    meanFR000(i,:) = mean(td_act000(i).S1_spikes);
end
% dir = 0;
% idx000 = abs(PDs - dir) <= 22.5;
% meanFR000 = mean(meanFR000(:,idx000),2);

td_act045 = td_act(tgtDirs_act==45);
for i = 1:numel(td_act045)
    meanSpeed045(i,:) = mean(td_act045(i).speed);
    meanFR045(i,:) = mean(td_act045(i).S1_spikes);
end
% dir = 45;
% idx045 = abs(PDs - dir) <= 22.5;
% meanFR045 = mean(meanFR045(:,idx045),2);

td_act090 = td_act(tgtDirs_act==90);
for i = 1:numel(td_act090)
    meanSpeed090(i,:) = mean(td_act090(i).speed);
    meanFR090(i,:) = mean(td_act090(i).S1_spikes);
end
% dir = 90;
% idx090 = abs(PDs - dir) <= 22.5;
% meanFR090 = mean(meanFR090(:,idx090),2);

td_act135 = td_act(tgtDirs_act==135);
for i = 1:numel(td_act135)
    meanSpeed135(i,:) = mean(td_act135(i).speed);
    meanFR135(i,:) = mean(td_act135(i).S1_spikes);
end
% dir = 135;
% idx135 = abs(PDs - dir) <= 22.5;
% meanFR135 = mean(meanFR135(:,idx135),2);

td_act180 = td_act(tgtDirs_act==180);
for i = 1:numel(td_act180)
    meanSpeed180(i,:) = mean(td_act180(i).speed);
    meanFR180(i,:) = mean(td_act180(i).S1_spikes);
end
% dir = 180;
% idx180 = abs(PDs - dir) <= 22.5;
% meanFR180 = mean(meanFR180(:,idx180),2);

td_act225 = td_act(tgtDirs_act==225);
for i = 1:numel(td_act225)
    meanSpeed225(i,:) = mean(td_act225(i).speed);
    meanFR225(i,:) = mean(td_act225(i).S1_spikes);
end
% dir = 225;
% idx225 = abs(PDs - dir) <= 22.5;
% meanFR225 = mean(meanFR225(:,idx225),2);


td_act270 = td_act(tgtDirs_act==270);
for i = 1:numel(td_act270)
    meanSpeed270(i,:) = mean(td_act270(i).speed);
    meanFR270(i,:) = mean(td_act270(i).S1_spikes);
end
% dir = 270;
% idx270 = abs(PDs - dir) <= 22.5;
% meanFR270 = mean(meanFR270(:,idx270),2);

td_act315 = td_act(tgtDirs_act==315);
for i = 1:numel(td_act315)
    meanSpeed315(i,:) = mean(td_act315(i).speed);
    meanFR315(i,:) = mean(td_act315(i).S1_spikes);
end
% dir = 315;
% idx315 = abs(PDs - dir) <= 22.5;
% meanFR315 = mean(meanFR315(:,idx315),2);

actFR = mean(cat(1,meanFR000,meanFR045,meanFR090,meanFR135,meanFR180,...
    meanFR225,meanFR270,meanFR315),2);
actSpeed = cat(1,meanSpeed000,meanSpeed045,meanSpeed090,meanSpeed135,...
    meanSpeed180,meanSpeed225,meanSpeed270,meanSpeed315);


edges = 0:3:50;

Y = discretize(actSpeed,edges);
for i = 1:numel(edges)
    actFR_bins(i) = mean(actFR(Y==i));
    actFR_ct(i) = sum(Y==i);
    actFR_std(i) = std(actFR(Y==i));

end

actFR_bins(actFR_ct < 2) = NaN;
actFR_std(actFR_ct < 2) = NaN;

set(gcf,'Color','White','Units','Normalized','position',[0.25 0.25 0.6 0.5])
set(gca,'FontName','Helvetica')
title('Mean S1 firing rates versus speed')
subplot(4,1,1:2); hold on

plot(edges,actFR_bins,'-b.','Markersize',10); hold on;
errorbar(edges,actFR_bins,-actFR_std,actFR_std,'linestyle','none','color','b')

% plot(actSpeed,actFR,'b.','MarkerSize',0.1)

subplot(4,1,3); hold on; box off;

n = histcounts(actSpeed,edges);
histogram('bincounts',n,'binedges',edges,'Normalization','probability','FaceColor',[0 0 1])






% plot(1:0.01:45*pow, (1:0.01:45*pow).^(0.35/pow)+8,'k')



clearvars -except td_act td_bump td_comb PDs edges

tgtDirs_pas = vertcat(td_bump.bumpDir);
tgtDirs_pas(tgtDirs_pas>=360) = tgtDirs_pas(tgtDirs_pas>=360) - 360;

td_pas000 = td_bump(tgtDirs_pas==0);
for i = 1:numel(td_pas000)
    meanSpeed000(i,:) = mean(td_pas000(i).speed);
    meanFR000(i,:) = mean(td_pas000(i).S1_spikes);
end
% dir = 0;
% idx000 = abs(PDs - dir) <= 22.5;
% meanFR000 = mean(meanFR000(:,idx000),2);


td_pas045 = td_bump(tgtDirs_pas==45);
for i = 1:numel(td_pas045)
    meanSpeed045(i,:) = mean(td_pas045(i).speed);
    meanFR045(i,:) = mean(td_pas045(i).S1_spikes);
end
% dir = 45;
% idx045 = abs(PDs - dir) <= 22.5;
% meanFR045 = mean(meanFR045(:,idx045),2);

td_pas090 = td_bump(tgtDirs_pas==90);
for i = 1:numel(td_pas090)
    meanSpeed090(i,:) = mean(td_pas090(i).speed);
    meanFR090(i,:) = mean(td_pas090(i).S1_spikes);
end
% dir = 90;
% idx090 = abs(PDs - dir) <= 22.5;
% meanFR090 = mean(meanFR090(:,idx090),2);

td_pas135 = td_bump(tgtDirs_pas==135);
for i = 1:numel(td_pas135)
    meanSpeed135(i,:) = mean(td_pas135(i).speed);
    meanFR135(i,:) = mean(td_pas135(i).S1_spikes);
end
% dir = 135;
% idx135 = abs(PDs - dir) <= 22.5;
% meanFR135 = mean(meanFR135(:,idx135),2);

td_pas180 = td_bump(tgtDirs_pas==180);
for i = 1:numel(td_pas180)
    meanSpeed180(i,:) = mean(td_pas180(i).speed);
    meanFR180(i,:) = mean(td_pas180(i).S1_spikes);
end
% dir = 180;
% idx180 = abs(PDs - dir) <= 22.5;
% meanFR180 = mean(meanFR180(:,idx180),2);

td_pas225 = td_bump(tgtDirs_pas==225);
for i = 1:numel(td_pas225)
    meanSpeed225(i,:) = mean(td_pas225(i).speed);
    meanFR225(i,:) = mean(td_pas225(i).S1_spikes);
end
% dir = 225;
% idx225 = abs(PDs - dir) <= 22.5;
% meanFR225 = mean(meanFR225(:,idx225),2);

td_pas270 = td_bump(tgtDirs_pas==270);
for i = 1:numel(td_pas270)
    meanSpeed270(i,:) = mean(td_pas270(i).speed);
    meanFR270(i,:) = mean(td_pas270(i).S1_spikes);
end
% dir = 270;
% idx270 = abs(PDs - dir) <= 22.5;
% meanFR270 = mean(meanFR270(:,idx270),2);

td_pas315 = td_bump(tgtDirs_pas==315);
for i = 1:numel(td_pas315)
    meanSpeed315(i,:) = mean(td_pas315(i).speed);
    meanFR315(i,:) = mean(td_pas315(i).S1_spikes);
end
% dir = 315;
% idx315 = abs(PDs - dir) <= 22.5;
% meanFR315 = mean(meanFR315(:,idx315),2);




pasFR = mean(cat(1,meanFR000,meanFR045,meanFR090,meanFR135,meanFR180,...
    meanFR225,meanFR270,meanFR315),2);
pasSpeed = cat(1,meanSpeed000,meanSpeed045,meanSpeed090,meanSpeed135,...
    meanSpeed180,meanSpeed225,meanSpeed270,meanSpeed315);
Y = discretize(pasSpeed,edges);
for i = 1:numel(edges)
    pasFR_bins(i) = mean(pasFR(Y==i));
    pasFR_ct(i) = sum(Y==i);
    pasFR_std(i) = std(pasFR(Y==i));

end

pasFR_bins(pasFR_ct < 5) = NaN;
pasFR_std(pasFR_ct < 5) = NaN;

subplot(4,1,1:2); hold on

plot(edges,pasFR_bins,'-r.','Markersize',10); hold on;
errorbar(edges,pasFR_bins,-pasFR_std,pasFR_std,'linestyle','none','Color','r')

% plot(pasSpeed,pasFR,'r.','MarkerSize',0.1)

subplot(4,1,3); hold on

n = histcounts(pasSpeed,edges);
histogram('bincounts',n,'binedges',edges,'Normalization','probability','FaceColor',[1 0 0],'FaceAlpha',0.8)


clearvars -except td_act td_bump td_comb PDs edges

tgtDirs_comb = vertcat(td_comb.tgtDir);
bumpDirs_comb = vertcat(td_comb.bumpDir);
totalDirs_comb = bumpDirs_comb - tgtDirs_comb;


while sum(totalDirs_comb >= 360) > 0
    totalDirs_comb(totalDirs_comb>=360) = totalDirs_comb(totalDirs_comb>=360) - 360;
end



td_comb000 = td_comb(totalDirs_comb==0);
for i = 1:numel(td_comb000)
    meanSpeed000(i,:) = mean(td_comb000(i).speed);
    meanFR000(i,:) = mean(td_comb000(i).S1_spikes);
end
% dir = 0;
% idx000 = abs(PDs - dir) <= 22.5;
% meanFR000 = mean(meanFR000(:,idx000),2);
% 


td_comb045 = td_comb(totalDirs_comb==45);
for i = 1:numel(td_comb045)
    meanSpeed045(i,:) = mean(td_comb045(i).speed);
    meanFR045(i,:) = mean(td_comb045(i).S1_spikes);
end
% dir = 45;
% idx045 = abs(PDs - dir) <= 22.5;
% meanFR045 = mean(meanFR045(:,idx045),2);

td_comb090 = td_comb(totalDirs_comb==90);
for i = 1:numel(td_comb090)
    meanSpeed090(i,:) = mean(td_comb090(i).speed);
    meanFR090(i,:) = mean(td_comb090(i).S1_spikes);
end
% dir = 90;
% idx090 = abs(PDs - dir) <= 22.5;
% meanFR090 = mean(meanFR090(:,idx090),2);

td_comb135 = td_comb(totalDirs_comb==135);
for i = 1:numel(td_comb135)
    meanSpeed135(i,:) = mean(td_comb135(i).speed);
    meanFR135(i,:) = mean(td_comb135(i).S1_spikes);
end
% dir = 135;
% idx135 = abs(PDs - dir) <= 22.5;
% meanFR135 = mean(meanFR135(:,idx135),2);

td_comb180 = td_comb(totalDirs_comb==180);
for i = 1:numel(td_comb180)
    meanSpeed180(i,:) = mean(td_comb180(i).speed);
    meanFR180(i,:) = mean(td_comb180(i).S1_spikes);
end
% dir = 180;
% idx180 = abs(PDs - dir) <= 22.5;
% meanFR180 = mean(meanFR180(:,idx180),2);


td_comb225 = td_comb(totalDirs_comb==225);
for i = 1:numel(td_comb225)
    meanSpeed225(i,:) = mean(td_comb225(i).speed);
    meanFR225(i,:) = mean(td_comb225(i).S1_spikes);
end
% dir = 225;
% idx225 = abs(PDs - dir) <= 22.5;
% meanFR225 = mean(meanFR225(:,idx225),2);

td_comb270 = td_comb(totalDirs_comb==270);
for i = 1:numel(td_comb270)
    meanSpeed270(i,:) = mean(td_comb270(i).speed);
    meanFR270(i,:) = mean(td_comb270(i).S1_spikes);
end
% dir = 270;
% idx270 = abs(PDs - dir) <= 22.5;
% meanFR270 = mean(meanFR270(:,idx270),2);

td_comb315 = td_comb(totalDirs_comb==0315);
for i = 1:numel(td_comb315)
    meanSpeed315(i,:) = mean(td_comb315(i).speed);
    meanFR315(i,:) = mean(td_comb315(i).S1_spikes);
end
% dir = 315;
% idx315 = abs(PDs - dir) <= 22.5;
% meanFR315 = mean(meanFR315(:,idx315),2);



% figure; hold on;

% combFR = mean(cat(1,meanFR000,meanFR045,meanFR090,meanFR135,meanFR180,...
%     meanFR225,meanFR270,meanFR315),2);
% combSpeed = cat(1,meanSpeed000,meanSpeed045,meanSpeed090,meanSpeed135,...
%     meanSpeed180,meanSpeed225,meanSpeed270,meanSpeed315);

% combFR = mean(cat(1,meanFR000,meanFR090,meanFR180,meanFR270),2);
% combSpeed = cat(1,meanSpeed000,meanSpeed090,meanSpeed180,meanSpeed270);
combFR = mean(cat(1,meanFR180),2);
combSpeed = cat(1,meanSpeed180);
% combFR = mean(meanFR180,2);
% combSpeed = meanSpeed180;


%USE THESE IF DOING RELATIVE DIRECTIONS
% combFR_assist = mean(meanFR000,2);
% combSpeed_assist = meanSpeed000;

combFR_resist = mean(meanFR180,2);
combSpeed_resist = meanSpeed180;

% combFR_perpen = mean(cat(1,meanFR090,meanFR270),2);
% combSpeed_perpen = cat(1,meanSpeed090,meanSpeed270);



Y = discretize(combSpeed,edges);
for i = 1:numel(edges)
    combFR_bins(i) = mean(combFR(Y==i));
    combFR_ct(i) = sum(Y==i);
    combFR_sem(i) = std(combFR(Y==i));

end

combFR_bins(combFR_ct < 5) = NaN;
combFR_std(combFR_ct < 5) = NaN;


combColor = [0.3 0.8 0.4];
subplot(4,1,1:2); hold on
set(gca,'tickdir','out')
plot(edges,combFR_bins,'-o','Color','k','MarkerFaceColor','k','Markersize',5); hold on;
errorbar(edges,combFR_bins,-combFR_sem,combFR_sem,'linestyle','none','Color','k')
set(gca,'ylim',[5 15])
% plot(combSpeed_assist,combFR_assist,'.','Color',[0.8 0.2 0.8],'MarkerSize',5)
plot(combSpeed_resist,combFR_resist,'.','Color',[0.3 0.7 0.3],'MarkerSize',5)
% plot(combSpeed_perpen,combFR_perpen,'.','Color',[0.5 0.5 0.5],'MarkerSize',0.1)

n = histcounts(combSpeed,edges);
subplot(4,1,3); hold on
set(gca,'xtick',[])
histogram('bincounts',n,'binedges',edges,'Normalization','probability','FaceColor',[0 0 0])


subplot(4,1,4); hold on
set(gca,'ylim',[0 0.3],'ytick',[0 0.3],'xtick',[],'box','off','tickdir','out')
histogram(combSpeed_assist,edges,'Normalization','probability','FaceColor',[0.8 0.2 0.8])
histogram(combSpeed_resist,edges,'Normalization','probability','FaceColor',[0.3 0.7 0.3])
% histogram(combSpeed_perpen,edges,'Normalization','probability','FaceColor',[0.5 0.5 0.5])
% histogram(combSpeed,edges,'Normalization','probability','FaceColor',[0 0 0])



%% Plot firing rates  for move-bumps

tgtDirs = 0:45:315;
relBumpDirs = 180;
numTgt = numel(tgtDirs);
numBump = numel(relBumpDirs);
numCond = numTgt*numBump;

meanFR_comb = [];
predmeanFR_comb = [];
meanFR_pred = [];
meanFR_glm = [];

plotData = 1;

% neuron2plot = 46

for neuron2plot = 15%[33]%1:size(td_comb(1).S1_spikes,2) %36 40
    clear FR predFR
    FR(numel(tgtDirs),numel(relBumpDirs)) = struct();
    FR(1,1).rate = [];
    
    predFR(numel(tgtDirs),numel(relBumpDirs)) = struct();
    predFR(1,1).rate = [];
    
    if plotData == 1
        ratefig = figure;
        set(gcf,'Color','White')
    end
    i = 0;
    j = 0;
    for rd = tgtDirs
        i = i+1;
        j = 0;
        for rbd = relBumpDirs
            j = j+1;
            for trial = 1:numel(td_comb)
%                 plotStartIdx = td_comb(trial).idx_movement_on - 50;
%                 plotStopIdx = td_comb(trial).idx_movement_on + 50;
                plotStartIdx = 1;
                plotStopIdx = 14;
                absBumpDir = td_comb(trial).bumpDir;
                relBumpDir = td_comb(trial).bumpDir - td_comb(trial).tgtDir;
                
                if relBumpDir < 0
                    relBumpDir = relBumpDir + 360;
                end
                if td_comb(trial).tgtDir == rd && relBumpDir == rbd
                    if plotData == 1
                        set(0,'CurrentFigure',ratefig)
                        set(gcf,'visible','off')
                        subplot(numTgt,numBump,(i-1)*(numBump)+j); hold on;
                        axis([0 0.13 0 100])
                    end
                    FR(i,j).rate(end+1,:) = td_comb(trial).S1_spikes(plotStartIdx:plotStopIdx,neuron2plot);
%                     FR(i,j).rate(end,:) = FR(i,j).rate(end,:) - mean(FR(i,j).rate(end,:));
                    if plotData == 1; plot(-0:td_comb(trial).bin_size:0.13,squeeze(FR(i,j).rate(end,:)),'color',[1 0.7 0.9]); end
                    
%                     predFR(i,j).rate(end+1,:) = td_comb(trial).glm_act2comb(plotStartIdx:plotStopIdx,neuron2plot) + td_comb(trial).glm_bump2comb(plotStartIdx:plotStopIdx,neuron2plot);
%                     predFR(i,j).rate(end,:) = predFR(i,j).rate(end,:) - mean(predFR(i,j).rate(end,:));
%                     if plotData == 1; plot(-2:12,squeeze(predFR(i,j).rate(end,:)),'color',[1 0.9 0.9]); end
                    
                    
                    % Formatting
                    if plotData == 1
                        set(gca,'TickDir','out','clipping','off')
                        if j ~= 1
                            set(gca,'yticklabel',[])
                        else
                            ylabel(num2str(rd))
                        end
                        if i ~= numTgt
                            set(gca,'xticklabel',[])
                        end
                        if i == 1
                            title(num2str(rbd))
                        end
                    end
                end
            end
            meanFR = [];
            meanFR = mean(FR(i,j).rate);
            
%             predmeanFR = [];
%             predmeanFR = mean(predFR(i,j).rate);
            
            
            if plotData == 1
                set(0,'CurrentFigure',ratefig)
                plot(-0:td_comb(1).bin_size:0.13,squeeze(meanFR),'linewidth',2,'color',[0.9 0.3 0.6])
%                 plot(-2:12,squeeze(predmeanFR),'linewidth',2,'color',[1 0 0])
            end
            
            meanFR_comb(neuron2plot,i,j,:) = meanFR;
%             predmeanFR_comb(neuron2plot,i,j,:) = predmeanFR;
        end
    end
    if plotData == 1
        set(ratefig,'visible','on')
    end
end

%% Plot bump-only trials

relBumpDirs = 0:45:315;
tgtDirs = 1;
numTgt = numel(tgtDirs);
numBump = numel(relBumpDirs);
numCond = numTgt*numBump;

meanFR_bump = [];

plotData = 1;

for neuron2plot = 15%[33]%1:size(td_bump(1).S1_spikes,2) %36 40
    clear FR
    FR(numel(tgtDirs),numel(relBumpDirs)) = struct();
    FR(1,1).rate = [];
    if plotData == 1
        ratefig = figure;
    end
    i = 1;
    j = 0;
    for rbd = relBumpDirs
        j = j+1;
        for trial = 1:numel(td_bump)
            plotStartIdx = 1;
            plotStopIdx = 14;
            absBumpDir = td_bump(trial).bumpDir;
            relBumpDir = td_bump(trial).bumpDir - 180; %Don't care about the target dir for passive
            
            if relBumpDir < 0
                relBumpDir = relBumpDir + 360;
            end
            if relBumpDir >= 360
                relBumpDir = relBumpDir - 360;
            end
            
            if relBumpDir == rbd
                if plotData == 1
                    
                    set(0,'CurrentFigure',ratefig)
                    set(gcf,'visible','off','Color','white')
                    subplot(numBump,numTgt,(i-1)*(numBump)+j); hold on;
                    axis([0 0.13 0 150])
                end
                FR(i,j).rate(end+1,:) = td_bump(trial).S1_spikes(plotStartIdx:plotStopIdx,neuron2plot);
                
                if plotData == 1
                    
                    plot(0:td_bump(trial).bin_size:0.13,squeeze(FR(i,j).rate(end,:)),'color',[0.6 1 0.8])
                    % Formatting
                    set(gca,'TickDir','out','clipping','off')
                    if j ~= 1
                        set(gca,'yticklabel',[])
                    else
                        ylabel('Bumps')
                    end
                    if i ~= numTgt
                        set(gca,'xticklabel',[])
                    end

                end
            end
        end
        meanFR = [];
        meanFR = mean(FR(i,j).rate);
        if plotData == 1
            
            set(0,'CurrentFigure',ratefig)
            plot(0:td_bump(trial).bin_size:0.13,squeeze(meanFR),'linewidth',2,'color',[0.3 0.7 0.5])
        end
        meanFR_bump(neuron2plot,i,j,:) = meanFR;
    end
    if plotData == 1
        
        set(ratefig,'visible','on')
    end
end


%% Plot Active trials
relBumpDirs = 1;
tgtDirs = 0:45:315;
numTgt = numel(tgtDirs);
numBump = numel(relBumpDirs);
numCond = numTgt*numBump;

meanFR_act = [];

plotData = 1;

for neuron2plot = 15%[33]%= 1:size(td_act(1).S1_spikes,2) %36 40
    clear FR
    FR(numel(tgtDirs),numel(relBumpDirs)) = struct();
    FR(1,1).rate = [];
    if plotData == 1
        ratefig = figure;
    end
    i = 0;
    j = 1;
    for rd = tgtDirs
        i = i+1;
        for trial = 1:numel(td_act)
%             plotStartIdx = td_act(trial).idx_movement_on - 50;
%             plotStopIdx = td_act(trial).idx_movement_on + 50;
            plotStartIdx = 1;
            plotStopIdx = 14;
            absBumpDir = td_act(trial).bumpDir;
            tgtDir = td_act(trial).tgtDir; %Don't care about the target dir for passive
            
            if tgtDir < 0
                tgtDir = tgtDir + 360;
            end
            if tgtDir >= 360
                tgtDir = tgtDir - 360;
            end
            
            if tgtDir == rd
                
                if plotData == 1
                    set(0,'CurrentFigure',ratefig)
                    set(gcf,'visible','off','Color','white')
                    subplot(numTgt,numBump,(i-1)*(numBump)+j); hold on;
                    axis([-0 0.13 0 100])
                    
                end
                
                FR(i,j).rate(end+1,:) = td_act(trial).S1_spikes(plotStartIdx:plotStopIdx,neuron2plot);
                
                
                if plotData == 1
                    plot(0:td_act(trial).bin_size:0.13,squeeze(FR(i,j).rate(end,:)),'color',[0.9 0.9 1])
                    % Formatting
                    set(gca,'TickDir','out','clipping','off')
                    if j ~= 1
                        set(gca,'yticklabel',[])
                    else
                        ylabel('Reaches')
                    end
                    if i ~= numTgt
                        set(gca,'xticklabel',[])
                    end

                end
            end
        end
        meanFR = [];
        meanFR = mean(FR(i,j).rate);
        if plotData == 1
            set(0,'CurrentFigure',ratefig)
            plot(-0:td_act(trial).bin_size:0.13,squeeze(meanFR),'linewidth',2,'color',[0.2 0.2 0.7])
        end
        meanFR_act(neuron2plot,i,j,:) = meanFR;
    end
    if plotData == 1
        set(ratefig,'visible','on')
    end
end


%% Predict bump-move S1 firing rates from bump and act firing rates

rd = 0:45:315;
bd = 180;
for neuron = 67
    for i = 1:numel(rd)
        bumpIdx = i:2:i+6;
        bumpIdx(bumpIdx>8) = bumpIdx(bumpIdx>8) - 8;
        reorderedFRbump(neuron,i,:,:) = meanFR_bump(neuron,:,bumpIdx,:);
    end
end
reorderedFRact = cat(3,meanFR_act,meanFR_act,meanFR_act,meanFR_act);
FR_pred = (reorderedFRbump - mean(reorderedFRbump)) + (reorderedFRact - mean(reorderedFRact,3));
% FR_pred = (reorderedFRbump) + (reorderedFRact) - mean(reorderedFRbump);
FR_actual = meanFR_comb;

for neuron = 67
    figure;
    for i = 1:8
        subplot(8,1,i), hold on
        axis([0 0.13 0 100])
        plot(0:0.01:0.13,squeeze(FR_pred(neuron,i,1,:)),'color','k','LineWidth',2)
        plot(0:0.01:0.13,squeeze(FR_actual(neuron,i,1,:)),'color',[0.9 0.3 0.6],'LineWidth',2)
    end
end
        

% n = 0;
% for i = 1:8
%     for j = 1:4
%         n = n+1;
%         for neuron = 33
%             CC_FR = corrcoef(squeeze(meanFR_comb(neuron,i,j,:)), squeeze(FR_pred(neuron,i,j,:)));
%             CC_GLM = corrcoef(squeeze(meanFR_comb(neuron,i,j,:)), squeeze(predmeanFR_comb(neuron,i,j,:)));
%             R2_FR(i,j,neuron) = CC_FR(2)^2;
%             R2_GLM(i,j,neuron) = CC_GLM(2)^2;
%         end
%         subplot(8,4,n), hold on
%         notBoxPlot(squeeze(R2_FR(i,j,:))-squeeze(R2_GLM(i,j,:)),ones(1,85),0.1,'patch')
% %         notBoxPlot(squeeze(R2_GLM(i,j,:)),2*ones(1,85),0.1,'patch')
%         
%         
%         [h(i,j),p(i,j)] = ttest(squeeze(R2_FR(i,j,:)),squeeze(R2_GLM(i,j,:)));
%     end
% end
% 
% 


%% Plot S1 FRs
clearvars -except td_act td_bump td_comb td_act_prep td_bl neurons_grp1 neurons_grp2 neurons_grp3 neurons_grp4
clc
tgtDirs_comb = vertcat(td_comb.tgtDir);
bumpDirs_comb = vertcat(td_comb.bumpDir);
totalDirs_comb = bumpDirs_comb - tgtDirs_comb;

while sum(totalDirs_comb >= 360) > 0 
    totalDirs_comb(totalDirs_comb>=360) = totalDirs_comb(totalDirs_comb>=360) - 360;
end


td_assist = td_comb(totalDirs_comb==0);
td_resist = td_comb(totalDirs_comb==180);
trajDiff_res = abs(wrapTo180(vertcat(td_resist.trajDiff)));


res_peakSpeeds = vertcat(td_resist.idx_peak_speed);
res_idxBumps = vertcat(td_resist.idx_bumpTime);

td_resist_early = td_resist(res_idxBumps < res_peakSpeeds);
td_resist_late = td_resist(res_idxBumps > res_peakSpeeds);


for i = 1:numel(td_resist)
   while td_resist(i).meanTraj >=360
       td_resist(i).meanTraj = td_resist(i).meanTraj - 360;
   end
   while td_resist(i).meanTraj < 0
       td_resist(i).meanTraj = td_resist(i).meanTraj + 360;
   end
   while td_resist(i).forceTraj >=360
       td_resist(i).forceTraj = td_resist(i).forceTraj - 360;
   end
   while td_resist(i).forceTraj < 0
       td_resist(i).forceTraj = td_resist(i).forceTraj + 360;
   end
end

td_resist_slowDown = td_resist(trajDiff_res<=90);
td_resist_turnAround = td_resist(trajDiff_res>90);


% for i = 1:numel(td_resist_early)
%    while td_resist_early(i).meanTraj >=360
%        td_resist_early(i).meanTraj = td_resist_early(i).meanTraj - 360;
%    end
%    while td_resist_early(i).meanTraj < 0
%        td_resist_early(i).meanTraj = td_resist_early(i).meanTraj + 360;
%    end
%    while td_resist_early(i).forceTraj >=360
%        td_resist_early(i).forceTraj = td_resist_early(i).forceTraj - 360;
%    end
%    while td_resist_early(i).forceTraj < 0
%        td_resist_early(i).forceTraj = td_resist_early(i).forceTraj + 360;
%    end
% end
% 
% for i = 1:numel(td_resist_late)
%    while td_resist_late(i).meanTraj >=360
%        td_resist_late(i).meanTraj = td_resist_late(i).meanTraj - 360;
%    end
%    while td_resist_late(i).meanTraj < 0
%        td_resist_late(i).meanTraj = td_resist_late(i).meanTraj + 360;
%    end
%    while td_resist_late(i).forceTraj >=360
%        td_resist_late(i).forceTraj = td_resist_late(i).forceTraj - 360;
%    end
%    while td_resist_late(i).forceTraj < 0
%        td_resist_late(i).forceTraj = td_resist_late(i).forceTraj + 360;
%    end
% end

for i = 1:numel(td_assist)
   while td_assist(i).meanTraj >=360
       td_assist(i).meanTraj = td_assist(i).meanTraj - 360;
   end
   while td_assist(i).meanTraj < 0
       td_assist(i).meanTraj = td_assist(i).meanTraj + 360;
   end
   while td_assist(i).forceTraj >=360
       td_assist(i).forceTraj = td_assist(i).forceTraj - 360;
   end
   while td_assist(i).forceTraj < 0
       td_assist(i).forceTraj = td_assist(i).forceTraj + 360;
   end
end

for i = 1:numel(td_bump)
   while td_bump(i).meanTraj >=360
       td_bump(i).meanTraj = td_bump(i).meanTraj - 360;
   end
   while td_bump(i).meanTraj < 0
       td_bump(i).meanTraj = td_bump(i).meanTraj + 360;
   end
   while td_bump(i).forceTraj >=360
       td_bump(i).forceTraj = td_bump(i).forceTraj - 360;
   end
   while td_bump(i).forceTraj < 0
       td_bump(i).forceTraj = td_bump(i).forceTraj + 360;
   end
end

for i = 1:numel(td_act)
   while td_act(i).meanTraj >=360
       td_act(i).meanTraj = td_act(i).meanTraj - 360;
   end
   while td_act(i).meanTraj < 0
       td_act(i).meanTraj = td_act(i).meanTraj + 360;
   end
   while td_act(i).forceTraj >=360
       td_act(i).forceTraj = td_act(i).forceTraj - 360;
   end
   while td_act(i).forceTraj < 0
       td_act(i).forceTraj = td_act(i).forceTraj + 360;
   end
end


plotTuning = 1;
meanFRs_RES = [];
meanFRs_ACT = [];
meanFRs_PAS = [];
meanFRs_RESTA = [];

mean_blFRs = mean(vertcat(td_bl.S1_spikes));

for context = [1,4]
    clear meanFRs_1 bumpDir_1 smoothFRs1 smoothSTD1 meanSens
    if context == 1
        td_1 = td_act;
        color1 = [0.9 0.3 0.6]; % RED
    elseif context == 2
        td_1 = td_act;
        color1 = [0.2 0.2 0.7]; % BLUE
    elseif context == 3
        td_1 = td_bump;
        color1 = [0.3 0.7 0.5]; % GREEN
    elseif context == 4
        td_1 = td_act_prep; 
        color1 = [0.6 0.2 0.6]; % PURPLE
    end
    
    
    nCell = size(td_1(1).S1_spikes,2);
    for cell = 1:nCell
        for i = 1:numel(td_1)
            meanFRs_1(cell,i) = mean(td_1(i).S1_spikes(:,cell));
%             meanSens(cell,i) = mean(abs(td_1(i).S1_spikes(:,cell) - mean_blFRs(cell))./td_1(i).speed);
            meanSens(cell,i) = mean(td_1(i).S1_spikes(:,cell)./td_1(i).speed);

            if context == 1
                bumpDir_1(cell,i) = wrapTo2Pi(deg2rad(td_1(i).tgtDir));
            elseif context == 2
                bumpDir_1(cell,i) = wrapTo2Pi(deg2rad(td_1(i).tgtDir));
            elseif context == 3
                bumpDir_1(cell,i) = wrapTo2Pi(deg2rad(td_1(i).bumpDir));
            else
                bumpDir_1(cell,i) = wrapTo2Pi(deg2rad(td_1(i).tgtDir));
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
%                 polarplot(bumpDir_1(cell,:),meanFRs_1(cell,:),'.','Color',color1); hold on;
%                 h1(cell) = polarplot((bumpDir_1(cell,:)),smoothFRs1(cell,:),'color',color1); hold on;
        %         h1(cell).LineWidth = 2;
        %         hs1_pos = polarplot((bumpDir_1(cell,:)),smoothFRs1(cell,:) + smoothSTD1(cell,:)/10,'color',color1);
        %         hs1_neg = polarplot((bumpDir_1(cell,:)),smoothFRs1(cell,:) - smoothSTD1(cell,:)/10,'color',color1);
        
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
                meanFRs_ACT(end+1,:) = Sens_bins;
            case 3
                meanFRs_PAS(end+1,:) = Sens_bins;
            case 4 
                meanFRs_RESTA(end+1,:) = Sens_bins;
        end
        
        clear FRs_bins FRs_sem clear FRs_ct
    end
end

%% Separating Neurons into groups 

%% clustering

for i = 1:numel(meanFRs_ACT,1)
    [~,maxACT] = max(meanFRs_ACT(i,:));
    meanFRs_ACT(i,:) = circshift(meanFRs_ACT(i,:),maxACT);
    meanFRs_PAS(i,:) = circshift(meanFRs_PAS(i,:),maxACT);
    meanFRs_RES(i,:) = circshift(meanFRs_RES(i,:),maxACT);

end

meanFRs_TOT = normalize([meanFRs_ACT(:,1:end-1) meanFRs_PAS(:,end-1) meanFRs_RES(:,end-1)],2,'range');

[COEFF, SCORE] = pca(meanFRs_TOT);

% plot3(SCORE(:,1),SCORE(:,2),SCORE(:,3),'.')


opts = statset('Display','final');
[idx,C] = kmeans(SCORE,5,'Replicates',1000,'Options',opts);

figure;
plot3(SCORE(idx==1,1),SCORE(idx==1,2),SCORE(idx==1,3),'r.','MarkerSize',12)
hold on
plot3(SCORE(idx==2,1),SCORE(idx==2,2),SCORE(idx==2,3),'b.','MarkerSize',12)
plot3(SCORE(idx==3,1),SCORE(idx==3,2),SCORE(idx==3,3),'g.','MarkerSize',12)
plot3(SCORE(idx==4,1),SCORE(idx==4,2),SCORE(idx==4,3),'c.','MarkerSize',12)

plot3(C(:,1),C(:,2),C(:,3),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off

figure;
silhouette(SCORE,idx)

neurons = 1:size(SCORE,1);
neurons_grp1 = neurons(idx==1);
neurons_grp2 = neurons(idx==2);
neurons_grp3 = neurons(idx==3);
neurons_grp4 = neurons(idx==4);


%% Split into Early and Late bumps
res_peakSpeeds = vertcat(td_resist.idx_peak_speed);
res_idxBumps = vertcat(td_resist.idx_bumpTime);

td_resist_early = td_resist(res_idxBumps < res_peakSpeeds);
td_resist_late = td_resist(res_idxBumps > res_peakSpeeds);
        


%% LDA on move bumps
clearvars -except td_act td_bump td_comb 

tgtDirs_comb = vertcat(td_comb.tgtDir);
bumpDirs_comb = vertcat(td_comb.bumpDir);
totalDirs_comb = bumpDirs_comb - tgtDirs_comb;
while sum(totalDirs_comb >= 360) > 0
    totalDirs_comb(totalDirs_comb>=360) = totalDirs_comb(totalDirs_comb>=360) - 360;
end

td_assist = td_comb(totalDirs_comb==0);
td_resist = td_comb(totalDirs_comb==180);
td_perpen = td_comb(totalDirs_comb==90 | totalDirs_comb==270);


for i = 1:numel(td_act)
    td_act(i).class = 0;
end


for j = 1:numel(td_bump)
    td_bump(j).class = 1;
end


for k = 1:numel(td_resist)
    td_resist(k).class = 2;
end

for l = 1:numel(td_assist)
    td_assist(l).class = 1;
end

for m = 1:numel(td_perpen)
    td_perpen(m).class = 2;
end

td_cat = cat(1,td_act,td_assist);
td_cat2 = cat(1,td_act,td_bump);

signal = zeros(numel(td_cat),size(td_cat(1).S1_spikes,2));
for trial = 1:numel(td_cat)
signal(trial,:) = mean(td_cat(trial).S1_spikes);
end

actpascomb = vertcat(td_cat.class);

mdl = fitcdiscr(signal,actpascomb);
class_pred = predict(mdl,signal);
separability = sum(class_pred == actpascomb)/length(actpascomb);

signal2 = zeros(numel(td_cat2),size(td_cat2(1).S1_spikes,2));
for trial = 1:numel(td_cat2)
signal2(trial,:) = mean(td_cat2(trial).S1_spikes);
end
actpascomb2 = vertcat(td_cat2.class);
class_pred2 = predict(mdl,signal2);
separability2 = sum(class_pred2 == actpascomb2)/length(actpascomb2);


% Plot results (from Raeed)
% plot active as filled, passive as open
% bump_colors = linspecer(4);
% act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
% pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;
act_color = [114 191 111]/256;
pas_color = [88 137 176]/256;
comb_color = [176 88 137]/256;

w = mdl.Sigma\diff(mdl.Mu)';
signal_sep = signal*w;

% get basis vector orthogonal to w for plotting
null_sep = null(w');
signal_null_sep = signal*null_sep;
[~,signal_null_sep_scores] = pca(signal_null_sep);

% plot for act/pas separability
figure(1); set(gcf,'Color','White')
hold all
scatter3(signal_sep(actpascomb==1),signal_null_sep_scores(actpascomb==1,1),signal_null_sep_scores(actpascomb==1,2),50,act_color,'filled')
scatter3(signal_sep(actpascomb==0),signal_null_sep_scores(actpascomb==0,1),signal_null_sep_scores(actpascomb==0,2),50,pas_color,'filled')
scatter3(signal_sep(actpascomb==2),signal_null_sep_scores(actpascomb==2,1),signal_null_sep_scores(actpascomb==2,2),50,comb_color,'filled')
ylim = get(gca,'ylim');
zlim = get(gca,'zlim');
% plot3([0 0],ylim,[0 0],'--k','linewidth',2)
% plot3([0 0],[0 0],zlim,'--k','linewidth',2)
set(gca,'box','off','tickdir','out','xtick',[],'ztick',[],'FontName','Helvetica')
xlabel('dim 1'), zlabel('dim 2')
title(['Act-Pas separability: ' num2str(round(separability*100)/100)],...
    'Fontname','Helvetica')
view([0 0])

% plot for predicted act/comb separability

w = mdl.Sigma\diff(mdl.Mu)';
signal_sep = signal2*w;

% get basis vector orthogonal to w for plotting
null_sep = null(w');
signal_null_sep = signal2*null_sep;
[~,signal_null_sep_scores] = pca(signal_null_sep);

figure(2); set(gcf,'Color','White')
hold all
scatter3(signal_sep(actpascomb2==1),signal_null_sep_scores(actpascomb2==1,1),signal_null_sep_scores(actpascomb2==1,2),50,act_color,'filled')
scatter3(signal_sep(actpascomb2==0),signal_null_sep_scores(actpascomb2==0,1),signal_null_sep_scores(actpascomb2==0,2),50,pas_color,'filled')
scatter3(signal_sep(actpascomb2==2),signal_null_sep_scores(actpascomb2==2,1),signal_null_sep_scores(actpascomb2==2,2),50,comb_color,'filled')
ylim = get(gca,'ylim');
zlim = get(gca,'zlim');
% plot3([0 0],ylim,[0 0],'--k','linewidth',2)
% plot3([0 0],[0 0],zlim,'--k','linewidth',2)
set(gca,'box','off','tickdir','out','xtick',[],'ztick',[],'FontName','Helvetica')
xlabel('dim 1'), zlabel('dim 2')
title(['Act-Comb separability: ' num2str(round(separability2*100)/100)],...
    'Fontname','Helvetica')
view([0 0])
