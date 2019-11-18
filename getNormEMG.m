function trial_data = getNormEMG(trial_data,params)
%
%

emgNorm = normalize(trial_data.emg,'range');

trial_data.emgNorm = emgNorm;


end