function trial_data = getRelMusLen(trial_data,params)
%
% TO-DOs: 
% - Pennation angles?
% - 

if params.opensimChris % Are we post processing Chris's data?
    % Set some default values
    if ~isfield(params,'idx_opensimLen')
        params.idx_opensimLen = 15:53;
    end
    
    if ~isfield(params,'idx_opensimVel')
        params.idx_opensimVel = 54:92;
    end
    
    if ~isfield(params,'L0')
        params.L0 = 'session_mean';
    end
    
    if strcmpi(params.L0,'session_mean')
        L0 = nanmean(trial_data.opensim(:,params.idx_opensimLen));
    elseif strcmpi(params.L0,'arbitrary')
        idx = params.L0idx;
        L0 = nanmean(trial_data.opensim(idx,params.idx_opensimLen));
    elseif strcmpi(params.L0,'trial_init')
        idx = trial_data.idx_startTime; %This should be before bump
        L0 = mean(trial_data.opensim(idx:idx+4,params.idx_opensimLen));
    end
    
else
    if strcmpi(params.L0,'session_mean')
        L0 = nanmean(trial_data.muscle_len);
    elseif strcmpi(params.L0,'arbitrary')
        idx = params.L0idx;
        L0 = nanmean(trial_data.muscle_len);
    elseif strcmpi(params.L0,'trial_init')
        idx = trial_data.idx_startTime; %This should be before bump
        L0 = mean(trial_data.muscle_len);
    end

    
end

if isfield(trial_data,'muscle_len') %Raeed's data
    musLenRel = trial_data.muscle_len./L0;
    musVelRel = trial_data.muscle_vel./L0;
    musNames = trial_data.muscle_names; %This is placeholder, may not work
elseif isfield(trial_data,'opensim') %Chris's data
    musLenRel = trial_data.opensim(:,params.idx_opensimLen)./L0;
    musVelRel = trial_data.opensim(:,params.idx_opensimVel)./L0;
    musNames = trial_data.opensim_names(params.idx_opensimLen);
    for i = 1:numel(musNames)
        musNames{i} = musNames{i}(1:end-4); %Get rid of '_len'
    end
end
    
%Add to trial_data struct
trial_data.musLenRel = musLenRel;
trial_data.musVelRel = musVelRel;
trial_data.musNames = musNames;

end