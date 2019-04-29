function GCM_total = monkeyPST_create_GCM(m, el, days, electrodes, spatial_ref, binned_trials, n_bins, clip, phase, model_architecture)
% m                 = 1 or 2
% el                = 1 (for preSMA) or 2 (for M1)
% days              = choose days of interest
% electrodes        = 'right_motor', 'left_motor'
% spatial_ref       = 'bipolar', 'nobipolar'
% binned_trials     = 0 or 1
% n_bins            = 8, 16
% phase             = 'search', 'repeat', 'search_repeat'
% model_architecture

%% Set Parameters and path
n_trials = n_bins;
monkey_init = {'r', 's'};

%% Create GCM
GCM_total = cell(0,0);
for m_indx = m
    analysis_dir = monkeyPST_make_analysis_dir(m_indx, el, electrodes, phase, spatial_ref, binned_trials, n_bins, clip, model_architecture);

    GCM = cell(0,0);
    
    for d = days
        clear GCM_session trial
        address = ['DCM_', monkey_init{m_indx}, '_', phase, '_d_', num2str(d),'_full'];
        fold2 = dir(fullfile([analysis_dir address '_*.mat']));
        GCM_session = {fold2.name}';
        GCM((length(GCM)+1):(length(GCM)+ length(GCM_session)),1) = GCM_session;
    end
    % remove empty
    GCM = GCM(~cellfun('isempty',GCM));
    GCM_total = [GCM_total; GCM];
end

end