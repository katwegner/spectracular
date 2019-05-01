%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HumanPST_analysis of time-varying spectra using Dynamic causal modeling (first level)
% and Parametric empirical Bayes (higher level analysis).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% About the Data:
% Participants were shown three stimuli (gray rectangle), of which one was
% the target. The task was to find the target by trial and error.
% The data were collected in 5 to 6 runs per subject.

% About the Analysis:
% A period of 1000 ms before stimulus onset is analyzed by investigating
% (1) whether spectral power changes within a run,
% (2) whether there are effective connections within the fronal cortex, and
% (3) whether effective connections change in accordance to the increase in spectral density.

% Format of Input Data:
% - task = 'preprocessing': bdf files to be preprocessed (task: 'preprocessing')
% - task = 'run_DCM': preprocessed EEG data in fieldtrip format and electrode location files
% - task = 'source_reconstruction': spm D structures
% - task = 'model_comp': DCM files
% - task = 'PEB3': DCM files

clear all
home_dir = pwd;
addpath ./functions
addpath ../functions
path_to_spm = which('spm');
addpath(path_to_spm(1:(end-5)));

%% choose parameters
% for first level analysis                                                  % choose among:
sub_index               = 1;                               % ALL: 1:22; R: [1:14 16:17 19:22], L: 15 AMB: 18
task                    = 'PEB1';                                        % 'preprocessing'run_DCM', 'model_comp', 'PEB3'
electrodes              = 'right_frontal';                                   % left_frontal, right_frontal
n_bins                  = 4;                                                % 4, 8 or 16
all_model_architecture  = {'no_conn', 'FWBW', 'FW'};                         % 'no_conn', 'FWBW','FW1'
effect                  = 'random';                                         % 'random' or 'fixed'

% for plotting or 2nd level analysis
model_architecture  = 'FWBW';                                                % 'no_conn', 'FWBW','FW1'
tbf                 = {'monoexp'};                  %'rawbeta1', 'rawbeta2', 'monoexp'
runs =               1:6;

%% Step 0: Create subject prefixes
subjects = [];
for subi = sub_index
    subjects{subi} = strcat('pp', sprintf('%02d', subi));
end

switch electrodes
    case 'left_frontal'
        load('coordinates_left_frontal_2el')
        positions = all_electrode_positions_left;
    case 'right_frontal'
        load('coordinates_right_frontal_2el')
        positions = all_electrode_positions_right;
    case 'left_frontal3'
       %positions = all_electrode_positions_left3;

end


switch task
    %% Step 1: preprocessing
    case 'preprocessing'
        for s = sub_index
            humanPST_preprocessing(subjects{s})
        end
        
        %% Step 2: single-trial DCMs per subject and day with chosen coordinates
    case 'run_DCM'
        for sub = sub_index
            for run = runs
                fname = humanPST_timeDCM_prep(subjects{sub}, n_bins, run);
                for model = 1:length(all_model_architecture)
                    model_architecture   = all_model_architecture{model};
                    humanPST_timeDCM_inv(subjects{sub}, fname, n_bins, run,electrodes, positions{sub}, model_architecture)
                end
            end
        end
        
        %% Step 3: Model comparison
    case 'model_comp'
        % random effect
        [F_matrix, exceedance_probabilities] = model_comparison_humanPST_random(sub_index, runs, electrodes, n_bins, all_model_architecture);
        if sum(exceedance_probabilities >= 0.90) == 1
            max_F_index = find(exceedance_probabilities == max(exceedance_probabilities));
            fprintf('% s is best model (random) \n.', all_model_architecture{max_F_index})
            exceedance_probabilities
        else
            fprintf('No best model found (random)\n')
            exceedance_probabilities
        end
        
        % fixed effect
        sum_F = sum(F_matrix);
        max_F_index = find(sum_F == max(sum_F));
        other_F_index = 1:3;
        other_F_index(max_F_index) = [];
        second_F_index = find(sum_F == max(sum_F(other_F_index)));
        if sum_F(max_F_index) - sum_F(second_F_index) > 3
            fprintf('% s is best model (fixed) \n.', all_model_architecture{max_F_index})
        else
            fprintf('No best model found (fixed)\n')
        end
        %% Step 4: Check increase in spectra over time on scalp level
    case 'time_effect_sensor'
        % see humanPST_freqanal
        
        %% Step 5: Source Reconstruction
    case 'source_reconstruction'
        % see gui: 'group inversion'
        
        
        %% Step 6: Check increase in spectra over time on source level
    case 'time_effect_source'
        effect = time_effect_source(4, length(sub_index), length(runs), electrodes, model_architecture);

        %% Step 7: Create temporal basis functions for PEB
    case 'tbf'
        get_temporal_function
        %% Step 8: PEB
    case 'PEB1'
        run = runs(1);
        [PEB, BMA] = humanPST_peb1(sub_index, run, electrodes, n_bins, tbf, model_architecture);
        %% Step 9: PEB of PEBs
    case 'PEB2'
        [PEB, PEB_all, BMA] = humanPST_peb2(sub_index, runs, electrodes, n_bins, tbf, model_architecture);
        %[PEB, BMA] = humanPST_peb_pooled(sub_index, runs, electrodes, n_bins, tbf, model_architecture)

    case 'PEB3'
        [PEB3, PEB_all3, BMA] = humanPST_peb3(sub_index, runs, electrodes, n_bins, tbf, model_architecture);
        
        
end

