function  analysis_dir = monkeyPST_make_analysis_dir(m, el, electrodes, phase, spatial_ref, binned_trials, n_bins, clip, model_architecture)
%-------------------------------------------------------------------------
% Trouble shooting
%-------------------------------------------------------------------------
if length(el) == 2 && ~exist('model_architecture','var')
   error('Specify model architecture')
end
%-------------------------------------------------------------------------
% Set parameters
%-------------------------------------------------------------------------
monkey_init = {'r', 's'};

%-------------------------------------------------------------------------
% Create preinversion directory and directory for DCM results
%-------------------------------------------------------------------------
if binned_trials == 0
    path_to_analysis = [pwd, '/analysis_DCM/results_c', num2str(clip), '_', phase, '/',  spatial_ref, '_', electrodes, '/'];
elseif binned_trials == 1
    path_to_analysis = [pwd, '/analysis_DCM/results_c', num2str(clip), '_', phase, '/',  spatial_ref, '_binned_trials', num2str(n_bins), '_', electrodes, '/'];
end

if length(el) == 1
    analysis_dir = strcat(path_to_analysis,  monkey_init{m}, '_el',num2str(el),'/');      
else 
    analysis_dir = strcat(path_to_analysis,  monkey_init{m}, '_', model_architecture, '/'); 
end

end