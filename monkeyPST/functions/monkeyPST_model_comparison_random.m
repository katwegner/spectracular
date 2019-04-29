function [F_matrix, exceedance_probabilities] = monkeyPST_model_comparison_random(m, el, days, electrodes, spatial_ref, binned_trials, n_bins, clip, phase, all_model_architecture)
home_dir = pwd;

% -------------------------------------------------------------------------
% GCM for all trials
% -------------------------------------------------------------------------

GCM = {};
for model = 1:length(all_model_architecture)
    GCM0 = {};
    model_architecture   = all_model_architecture{model};
    GCM0 = monkeyPST_create_GCM(m, el, days, electrodes, spatial_ref, binned_trials, n_bins, clip, phase, model_architecture);
    GCM = [GCM GCM0];
end

%% make GCM into F matrix

F_matrix = [];

for model = 1:length(all_model_architecture)
    model_architecture   = all_model_architecture{model};
        
    k = 0;
    for m_index = m

        analysis_dir = monkeyPST_make_analysis_dir(m_index, el, electrodes, phase, spatial_ref, binned_trials, n_bins, clip, model_architecture);
        cd(analysis_dir)
        for d = days
            for b = 1:n_bins
            k = k+1;
            fname = GCM{k,model};
            load(fname)
            F_matrix(k,model) = DCM.F;
            end
        end
        cd(home_dir)
    end
end


%% random effect model comparison
[alpha,exp_r,exceedance_probabilities,pxp,bor] = spm_BMS (F_matrix, [], 1);
% alpha   - vector of model probabilities
% exp_r   - expectation of the posterior p(r|y)
% xp      - exceedance probabilities
% pxp     - protected exceedance probabilities
% bor     - Bayes Omnibus Risk (probability that model frequencies
%           are equal)
end