function regressors = humanPST_create_PEB_regressor(s, run,  electrodes, n_bins, spectra, model_architecture)
%% info needed
home_dir = pwd;
cd(strcat(home_dir, filesep, 'analysis_DCM', filesep, electrodes))
%% create temporal basis functions
mean_regressor = ones(n_bins,1);

% beta
try

load beta1_raw_2el
rawbeta1_temporal_basis_function = squeeze(log(real(csd_beta1_mean(1:n_bins,run,s,:)).^2));

%% center beta's
for b = 1:size(rawbeta1_temporal_basis_function,2)
    rawbeta1_temporal_basis_function_centered(:,b) = rawbeta1_temporal_basis_function(:,b) - mean(rawbeta1_temporal_basis_function(:,b));
end


end
% alpha
try
    load alpha_raw
    rawalpha_temporal_basis_function = zscore(log(real(csd_alpha_mean(1:n_trials(d,m),d,m,el)).^2));
catch
end

% monoexponential
tau =  floor(n_bins/2);

% task
task_regressor = [1;-1];


%% create regressors
X = [];
cos_index = strcmp(spectra, 'cosine');
monoexp_index = strcmp(spectra, 'monoexp');
damped_cosine_index = strcmp(spectra, 'damped_cosine');
regressors = mean_regressor;
for s = 1:length(spectra)
    switch spectra{s}
        case 'beta1'
            regressor = beta1_temporal_basis_function;
        case'rawbeta1'
            regressor = rawbeta1_temporal_basis_function_centered;
        case'beta2'
            regressor = beta2_temporal_basis_function;
        case'rawbeta2'
            regressor = rawbeta2_temporal_basis_function;
        case 'rawtheta'
            regressor = rawtheta_temporal_basis_function;
            %         figure; scatter(time_windows, rawtheta_temporal_basis_function);xlabel('time window'); ylabel('theta power')
        case 'rawalpha'
            regressor = rawalpha_temporal_basis_function;
            %         figure; scatter(time_windows, rawalpha_temporal_basis_function);xlabel('time window'); ylabel('alpha power')
        case 'cosine'
            which_cos = sum(cos_index(1:s));
            regressor = cosine_temporal_basis_function(:,which_cos);
        case 'monoexp'
            which_monoexp = sum(monoexp_index(1:s));
            regressor =  exp(-((1:n_bins)/tau(which_monoexp)))';
        case 'damped_cosine'
            which_damped_cosine = sum(damped_cosine_index(1:s));
            regressor =  exp(-(time_windows/tau(which_damped_cosine))).*cosine_temporal_basis_function(:,4);
        case 'task'
            regressor = task_regressor;
    end
    regressors = [regressors, regressor];
end
cd(home_dir)
end
