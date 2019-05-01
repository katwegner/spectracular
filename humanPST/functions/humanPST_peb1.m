function [PEB, BMA] = humanPST_peb1(s, run, electrodes, n_bins, spectra, model_architecture)

%% 1. Set parameters and path
% -------------------------------------------------------------------------
% Set paths and file names
% -------------------------------------------------------------------------
fprintf(['Subject: ', num2str(s) ,', Run: ', num2str(run),'\n\n'])
home_dir = pwd;
plot_parameters = 0;
analysis_dir =  strcat(home_dir, filesep, 'analysis_DCM', filesep, electrodes, filesep, model_architecture); 

%% 2.  Get GCM and create PEB regressors

GCM = humanPST_create_GCM(s, run, electrodes, n_bins, model_architecture);
clear PEB BMA parameters
M = struct();
M.X = humanPST_create_PEB_regressor(s, run,  electrodes, n_bins, spectra, model_architecture);
field = {'A' 'G' 'T'};

cd(analysis_dir)


% Estimate model
[PEB , GCM]    = spm_dcm_peb(GCM,M,field);

% Bayesian Model Averaging
BMA = spm_dcm_peb_bmc(PEB);

% plot parameters
if plot_parameters == 1
    n_parameters =  size(PEB.Ep,1);
    if size(M.X,2)>2
        figure;
        for regressor_index = 1:size(M.X,2)
            subplot(1,size(M.X,2), regressor_index)
            bar(BMA.Pp(1+ n_parameters*(regressor_index -1) : (n_parameters + n_parameters * (regressor_index -1)))); hold on
            plot(1:n_parameters, repmat(0.95, 1, n_parameters), '--r', 'MarkerSize', 14)
        end
        
        for regressor_index = 1:size(M.X,2)
            figure;
            var_for_ci = diag(BMA.Cp);
            spm_plot_ci(BMA.Ep(1+ n_parameters*(regressor_index -1) : (n_parameters + n_parameters * (regressor_index -1))),...
                var_for_ci(1+ n_parameters*(regressor_index -1) : (n_parameters + n_parameters * (regressor_index -1))))
            if regressor_index == 1
                %   title('constant', 'FontSize', 30)
                ylabel('Effect size','FontSize',30)
                %         else
                %             title(spectra{regressor_index - 1},'FontSize',30)
            end
            xlabel('Parameter', 'FontSize', 30)
            %         names = {'T(1)','T(2)','T(3)', 'T(4)', 'G(1)','G(2)','G(3)','G(4)','G(5)','G(6)','G(7)','G(8)','G(9)','G(10)'};
            
            %         if less_parameters == 1
            %             xticklabels(names(1,2,3,4,13, 6, 7, 8))
            %         else
            %             xticklabels({'T(1)','T(2)','T(3)', 'T(4)', 'G(1)','G(2)','G(3)','G(4)','G(5)','G(6)','G(7)','G(8)','G(9)','G(10)'})
            %         end
            axis([0 n_parameters+1 -2 3.5])
            axis square
        end
        
    end
end

cd(home_dir)
end