function [data_path, analysis_dir] = monkeyPST_create_directories(m, el, electrodes, phase, spatial_ref, binned_trials, n_bins, clip, model_architecture)

home_dir = pwd;
% -------------------------------------------------------------------------
% Set paths and file names
% -------------------------------------------------------------------------
switch phase
    case {'search','repeat'}
        data_path = strcat(home_dir, '/data_search_repeat');
    case 'complete'
        data_path = strcat(home_dir, '/data_complete');
    case 'accuracy'
        data_path = strcat(home_dir, '/data_accuracy');
end
analysis_dir = monkeyPST_make_analysis_dir(m, el, electrodes, phase, spatial_ref, binned_trials, n_bins, clip, model_architecture);
mkdir(analysis_dir);

end