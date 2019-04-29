function [F_matrix, exceedance_probabilities] = humanPST_model_comparison_random(subs, run, electrodes, n_bins, all_model_architecture)
home_dir = pwd;

% -------------------------------------------------------------------------
% GCM for all trials
% -------------------------------------------------------------------------

GCM = {};
for model = 1:length(all_model_architecture)
    GCM0 = {};
    model_architecture   = all_model_architecture{model};
    GCM0 = create_GCM_humanPST(subs, run, electrodes, n_bins, model_architecture);
    GCM = [GCM GCM0];
end

%% make GCM into F matrix

F_matrix = [];

for model = 1:length(all_model_architecture)
    model_architecture   = all_model_architecture{model};
    k = 0;
    path = './analysis_DCM';
    analysis_dir =strcat(path, '/', electrodes, '/', model_architecture, '/');
    cd(analysis_dir)
    
    
    for r = 1:size(GCM,1)
        k = k+1;
        fname = GCM{k,model};
        load(fname)
        F_matrix(k,model) = DCM.F;
    end
    cd(home_dir)
    
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