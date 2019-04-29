%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Find temporal function %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AIM: find how beta power changes over trials, for the average of the week
% or each day separately

%--------------------------------------------------------------------------
%% Set up analysis directory and initiate values
%--------------------------------------------------------------------------
home_dir = pwd;
subs = 1:22;
electrodes = 'right_frontal';
n_bins = 4;
n_runs = 6;
model_architecture = 'FWBW';
%--------------------------------------------------------------------------
%% Initiate
%--------------------------------------------------------------------------
csd_alpha_mean = NaN(n_bins, n_runs, length(subs), 2);
csd_beta1_mean = NaN(n_bins, n_runs, length(subs), 2);
csd_beta1_std = NaN(n_bins, n_runs, length(subs), 2);
csd_beta_mean = NaN(n_bins, n_runs, length(subs), 2);
csd_beta_std = NaN(n_bins, n_runs, length(subs), 2);
%--------------------------------------------------------------------------
%% Choose plots
%--------------------------------------------------------------------------

plot_beta1fit = 1;
plot_beta2fit = 0;
plot_beta2fit2 = 0;
plot_alphafit = 0;
plot_thetafit = 0;
%--------------------------------------------------------------------------
%% 2. spectral power matrix
%--------------------------------------------------------------------------
        analysis_dir = strcat(home_dir, '/analysis_DCM/', electrodes, '/', model_architecture);
        cd(analysis_dir)
for s = subs
    for el = 1:2

        for run = 1:n_runs
            
            address = strcat('DCM_pp', sprintf('%02d', s), '_bin4_run', num2str(run), '_full');
            for t = 1: n_bins %length(DCM.xY.csd)
                try
                    load([address,'_',num2str(t)])           
                    DCM_csd = DCM.xY.csd{1,1}(:,el,el);
                    csd_alpha_mean(t,run, s,el) = mean(DCM_csd(6:10));
                    csd_beta1_mean(t,run, s,el) = mean(DCM_csd(10:13)); % for beta1 = [12:15]
                    csd_beta1_std(t,run, s,el) = std(DCM_csd(10:13));
                    csd_beta_mean(t,run, s,el) = mean(DCM_csd(11:28));   % for beta = [13:30]
                    csd_beta_std(t,run, s,el) = std(DCM_csd(11:28));     
                catch
                    csd_alpha_mean(t,run, s,el) = NaN;
                    csd_beta1_mean(t,run, s,el) = NaN;
                    csd_beta1_std(t,run, s,el) = NaN;
                    csd_beta_mean(t,run, s,el) = NaN;
                    csd_beta_std(t,run, s,el) = NaN;
                end
            end
        end
    end
end
cd(strcat('../'))
%% define save names
switch model_architecture
    case '1el'
        add_model_architecture = '';
    otherwise
        add_model_architecture = '_2el';
end
save_name_beta1 = ['beta1_raw', add_model_architecture];
save_name_beta = ['beta_raw', add_model_architecture];

save_name_alpha = ['alpha_raw', add_model_architecture];
save_name_theta = ['theta_raw', add_model_architecture];


save(save_name_beta1, 'csd_beta1_mean')
save(save_name_beta, 'csd_beta_mean')
%save(save_name_theta, 'csd_theta_mean')
save(save_name_alpha, 'csd_alpha_mean')

cd(home_dir)
%end
