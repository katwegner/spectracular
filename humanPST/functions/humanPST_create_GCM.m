function GCM = humanPST_create_GCM(subs, run, electrodes, n_bins, model_architecture)


%% Set Parameters and path
home_dir = pwd;
n_trials = n_bins;
path = './analysis_DCM';
analysis_dir =strcat(path, '/', electrodes, '/', model_architecture, '/');
%cd(analysis_dir)
%% Create GCM
GCM = cell(0,0);%cell(n_bins*length(subs),1);
for s = subs
    if length(run) > 1
    clear GCM_session trial
    fold2 = dir(fullfile([analysis_dir 'DCM_', strcat('pp', sprintf('%02d', s)), '_bin', num2str(n_bins), '_run*','_full', '_*.mat']));
    elseif length(run) == 1
    fold2 = dir(fullfile([analysis_dir 'DCM_', strcat('pp', sprintf('%02d', s)), '_bin', num2str(n_bins), '_run', num2str(run),'_full', '_*.mat']));
    end
    GCM_session = {fold2.name}';
    GCM((length(GCM)+1):(length(GCM)+ length(GCM_session)),1) = GCM_session;
end
cd(home_dir)
end