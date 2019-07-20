function [td, field_names] = getPowerTD(td, params)
signals = [];
powers =  [];
keep_sign = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra parameters that aren't documented in the header


if nargin > 1, assignParams(who,params); end % overwrite defaults
   
td = check_td_quality(td);
bin_size = td(1).bin_size;
% make sure the signal input formatting is good
signals = check_signals(td(1),signals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check output field addition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
for trial = 1:length(td)
    for iSig = 1:size(signals,1)
        for j = 1:length(powers)
            field_extra{iSig} = ['_',replace(num2str(powers(j)), '.', '_')];
            data = getSig(td(trial),signals(iSig,1));
            negFlag = data <0;
            data = abs(data).^powers(j);
            if keep_sign
                negMat = ones(length(data(:,1)), length(data(1,:)));
                negMat(negFlag) = -1;
                data = data.*negMat;
            end
            td(trial).([signals{iSig,1} field_extra{iSig}]) = data;
            field_names{iSig, j} = [signals{iSig,1} field_extra{iSig}];
        end
    end
end

end