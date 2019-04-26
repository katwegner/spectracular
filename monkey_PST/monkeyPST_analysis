%% Time-varying DCM for ECoG data

% This script inverts DCMs per trial for the 1000 ms delay period before
% stimulus onset (step 1), compares competing models (step 2),
% extracts the log beta power per trial and uses PEB to
% identify consistency across days and monkeys. Data is recorded from two
% macaque monkeys (one left handed, one female) who performed a
% trial-and-error search task over days (monkey r for 19 days, monkey s for
% 21 days). More details about the data can be found in:

% Stoll, F. M., Wilson, C. R., Faraut, M. C., Vezoli, J., Knoblauch, K., &
% Procyk, E. (2015). The effects of cognitive control and time on frontal
% beta oscillations. Cerebral Cortex, 26(4), 1715-1732.


% Author: Katharina Wegner
% Last update: August, 15th  2018
%% add path
home_dir = pwd;
addpath([pwd, '/functions']);
path_to_spm = '/Users/katharinawegner/Documents/1.Ghent/Programs/spm12';               %%%%%%%%%%%% add spm path %%%%%%%%%%%%%%%%
addpath(path_to_spm);

%% choose parameters
clear all
task                = 'PEB2';            % 'run_DCM', 'model_comp', 'PEB3'
m                   = 1;                    % 1 for Risette, 2 for Satan, 1:2
el                  = 1:2;                  % 1(preSMA) or 2 (M1), 1:2
days                = 1:19;                 % span within 1 and 19
electrodes          = 'right_rfrontal_1114913';         %'left_lfrontal_47612'; %'right_rfrontal_1114913'    % 'right_motor' or 'left_motor' 'left_rfrontal_414613' 'left_rfrontal_21318', 'left_rfrontal_41268'
spatial_ref         = 'nobipolar';          % 'bipolar' or 'nobipolar'
phase               = 'search';             % 'search' or 'repeat'
binned_trials       = 1;                    % 0 or 1
n_bins              = 8;                    % 8 or 16
clip                = 0;                    % 0 (before), 1 (after 1st pause)
all_model_architecture = {'no_conn', 'FWBW', 'FW1'}; %'no_conn', 'FWBW', 'BWFW', 'FW1','FW2'
effect              = 'random';

% for plotting or 2nd level analysis
model_architecture  = 'FW1';
tbf                 = {'rawbeta1', 'rawbeta2', 'monoexp'};


%% choose analysis
% run_DCM             = 0;
% model_comp          = 1;
% trajectories        = 0;
% temporalbasis       = 0;
% beta_time_DCM       = 0;
% plot_beta           = 0;
% PEB1                = 0;
% PEB2                = 0;
% PEB3                = 0;


%% Step 1: single-trial DCMs per subject and day
switch task
    case 'run_DCM'
        for model = 1:length(all_model_architecture)
            model_architecture   = all_model_architecture{model};
            for c = clip
                for m_indx = m
                    [data_path, analysis_dir] = create_directories(m_indx, el, electrodes, phase, spatial_ref, binned_trials, n_bins, c, model_architecture);
                    for d = days
                       % try
                            switch phase
                                case {'search', 'repeat', 'complete'}
                                    [fname, DCM, analysis_dir,  Ns] = run_timeDCM_prep(m_indx, el, d, electrodes, phase, spatial_ref, binned_trials, n_bins, c, data_path, analysis_dir);
                                    run_timeDCM_inversion(el, analysis_dir, fname, DCM, Ns, model_architecture)
                                case {'search_repeat', 'repeat_search'}
                                    phase1 = {'search', 'repeat'};
                                    for p = 1:2
                                        [fname, DCM, analysis_dir,  Ns] = run_timeDCM_prep(m_indx, el, d, electrodes, phase1{p}, spatial_ref, binned_trials, n_bins, c, data_path, analysis_dir);
                                        run_timeDCM_inversion(el, analysis_dir, fname, DCM, Ns, model_architecture)
                                    end
                            end
%                         catch
%                            fprintf(['Inversion aborded!\n'], model, c, m_indx, d)
%                            %cd(home_dir)
%                         end
                    end
                end
            end
        end
        
        
        
        %% Step 2: Model comparison
    case 'model_comp'
        switch effect
            case 'fixed'
        for d = days
            switch phase
                case {'search', 'repeat', 'complete'}
                    [R2(d,:), F(d,:), index_best_model(d), winning(d), F_day(:,:,d)] = model_comparison(m, el, d, electrodes, phase, spatial_ref, binned_trials, n_bins, clip, all_model_architecture);
                case 'search_repeat'
                    phase1 = {'search', 'repeat'};
                    for p = 1:2
                        [R2(d+(p-1)*length(days),:), F(d+(p-1)*length(days),:), index_best_model(d+(p-1)*length(days)), winning(d+(p-1)*length(days))] = model_comparison(m, el, d, electrodes, phase1{p}, spatial_ref, binned_trials, n_bins, clip, all_model_architecture);
                    end
            end
        end
        index_best_model
        sum(F)
            case 'random'
                
                
                
        end
        
        %% Step 3:Plot DCM trajectories
    case 'trajectories'
        plot_DCM_trajectories_average(m, el, days, electrodes, spatial_ref, binned_trials, n_bins, clip, phase, model_architecture)
        
        
        %% Step 4: generate temporal basis functions
    case 'temporalbasis'
        find_temporal_function(days, electrodes, spatial_ref, binned_trials, n_bins, clip, phase, model_architecture)
        
    case 'PEB1'
        
        
        %% Step 7: PEB of PEBs
        
    case 'PEB2'
        plot_parameters = 0;
        [PEB, PEB_all, BMA] = run_peb_of_peb(m, el , days, electrodes, spatial_ref, binned_trials, n_bins, clip, phase, tbf, model_architecture, plot_parameters);
        %plot_peb_of_peb
        
        %% Step 7c: PEB of PEB of PEBs or multiple subjects
    case 'PEB3'
        plot_parameters = 0;
        m = 1:2;
        [PEB, PEB_all, BMA] = run_peb_of_peb_of_peb(m, el , days, electrodes, spatial_ref, binned_trials, n_bins, clip, phase, tbf, model_architecture, plot_parameters);
end
