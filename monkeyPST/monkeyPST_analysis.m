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
addpath ./functions
addpath ../functions
path_to_spm = which('spm');
addpath(path_to_spm(1:(end-5)));

%% choose parameters
clear all
task                = 'run_DCM';            % 'run_DCM', 'model_comp', 'PEB3'
m                   = 1;                    % 1 for Risette, 2 for Satan, 1:2
el                  = 1:2;                  % 1(preSMA) or 2 (M1), 1:2
days                = 1:19;                 % span within 1 and 19
electrodes          = 'right_frontal';      %'left_frontal'; %'right_frontal', 'right_motor', 'left_motor'
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

%% Step 1: single-trial DCMs per subject and day
switch task
    case 'run_DCM'
        for model = 1:length(all_model_architecture)
            model_architecture   = all_model_architecture{model};
            for c = clip
                for m_indx = m
                    [data_path, analysis_dir] = monkeyPST_create_directories(m_indx, el, electrodes, phase, spatial_ref, binned_trials, n_bins, c, model_architecture);
                    for d = days
                        % try
                        switch phase
                            case {'search', 'repeat', 'complete'}
                                [fname, DCM] = monkeyPST_timeDCM_prep(m_indx, el, d, electrodes, phase, spatial_ref, binned_trials, n_bins, c, data_path, analysis_dir);
                                monkeyPST_timeDCM_inversion(el, fname, DCM, analysis_dir, model_architecture)
                            case {'search_repeat', 'repeat_search'}
                                phase1 = {'search', 'repeat'};
                                for p = 1:2
                                    [fname, DCM] = monkeyPST_timeDCM_prep(m_indx, el, d, electrodes, phase{p}, spatial_ref, binned_trials, n_bins, c, data_path, analysis_dir);
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
%end(m, el, d, electrodes, phase, spatial_ref, binned_trials, n_bins, clip, data_path, analysis_dir)
home_dir = pwd;
monkey_name = {'Risette', 'Satan'};
monkey_init = {'r', 's'};

% -------------------------------------------------------------------------
% set unit index and n_trials
% -------------------------------------------------------------------------
n_trials = make_n_trials(binned_trials, n_bins, phase);

% -------------------------------------------------------------------------
% select electrodes
% -------------------------------------------------------------------------
switch spatial_ref
    case 'bipolar'
        switch electrodes
            case 'right_motor'
                if m == 1
                    electrode_index = [12 14; 19 20];
                elseif m == 2
                    electrode_index = [8 13;26 28];
                end
            case 'left_motor'
                if m == 1
                    electrode_index = [5 7; 15 16];
                elseif m == 2
                    electrode_index = [7 12; 21 27];
                end
            case 'right_left_motor'
                if m == 1
                    electrode_index = [5 7; 15 16; 12 14; 19 20];
                elseif m == 2
                    electrode_index = [7 12; 21 27;8 13;26 28];
                end
        end
    case 'nobipolar'
        switch electrodes
            case 'right_motor'
                if m == 1
                    electrode_index = [12 19];
                elseif m == 2
                    electrode_index = [8 26];
                end
            case 'left_motor'
                if m == 1
                    electrode_index = [5 15];
                elseif m == 2
                    electrode_index = [7 21];
                end
            case 'left_lfrontal'
                if m == 1
                    electrode_index = [2 6];
                elseif m == 2
                    electrode_index = [1 7];
                end
            case 'left_lfrontal_47612'
                if m == 1
                    electrode_index = [4 7];
                elseif m == 2
                    electrode_index = [6 12];
                end
            case 'right_rfrontal'
                if m == 1
                    electrode_index = [9 13];
                elseif m == 2
                    electrode_index = [4 13];
                end
            case 'right_rfrontal_1114913'
                if m == 1
                    electrode_index = [11 14];
                elseif m == 2
                    electrode_index = [9 13];
                end
            case 'right_lfrontal'
                if m == 1
                    electrode_index = [9 6];
                elseif m == 2
                    electrode_index = [4 7];
                end
            case 'right_lfrontal_96912'
                if m == 1
                    electrode_index = [9 6];
                elseif m == 2
                    electrode_index = [9 12];
                end
            case 'left_rfrontal'
                if m == 1
                    electrode_index = [2 13];
                elseif m == 2
                    electrode_index = [1 13];
                end
            case 'left_rfrontal_41268'
                if m == 1
                    electrode_index = [4 12];
                elseif m == 2
                    electrode_index = [6 8];
                end
            case 'left_rfrontal_414613'
                if m == 1
                    electrode_index = [4 14];
                elseif m == 2
                    electrode_index = [6 13];
                end
            case 'frontal'
                if m == 1
                    electrode_index = [4 11 7 14];
                elseif m == 2
                    electrode_index = [6 9 12 13];
                end
                
            otherwise
                fprintf('You have not specified the electrode location!')
        end
end

% -------------------------------------------------------------------------
% load data and name files
% -------------------------------------------------------------------------
cd(data_path)
file_name = [monkey_name{m}, '_c', num2str(clip), '_d_' num2str(d)];
load(file_name)
switch phase
    case 'search'
        data_phase = data_S;
    case 'repeat'
        data_phase = data_R ;
    case 'search_repeat'
        n_S = length(data_S)/781;
        n_R = length(data_R)/781;
        data_phase = [data_S, data_R];
    case 'complete'
        data_phase = data_C;
end
switch spatial_ref
    case 'bipolar' % bipolar montage
        for elec = 1:length(electrode_index)
            data_raw=data_phase(electrode_index(elec, :)', :); % S for search
            data_bip = data_raw(1,:) - data_raw(2,:);
            data1(elec, :) = data_bip;
        end
    case 'nobipolar'
        data1=data_phase(electrode_index, :); % S for search
end
fname = [analysis_dir 'spm_', monkey_init{m}, '_d_', num2str(d), '_c0_LFP_full'];                               % name of created spm file
DCM.name = make_address(m, d, phase);
cd(analysis_dir);
% -------------------------------------------------------------------------
% select electrode
% -------------------------------------------------------------------------
DATA2 = data1(el,:);
DATA1 = DATA2;
clear DATA2
clear data1

% -------------------------------------------------------------------------
% Set parameters
% -------------------------------------------------------------------------
[Nchannels,Nsamples,Ndays] = size(DATA1); % here trials is the number of days (or seizures etc)
Fsample=781;  %Hz                                                               % sampling frequency
Lw=1000; % in ms                                                                   % length of the window in ms
step=Lw;                                                                  % step of the sliding window in ms
Lw_p= Lw*Fsample/1000;                                                      % length of the window in samples
step_p=step*Fsample/1000;                                               % length of step in samples

% -------------------------------------------------------------------------
% Divide into [Lw] ms windows moving [step] ms at a time
% -------------------------------------------------------------------------
Ns = floor((Nsamples-Lw_p)/step_p) + 1;
for i=1:Ns
    l(i,:) = 1 + (i-1)*step_p:Lw_p + (i-1)*step_p ;
    DATA2(i,:,:) = DATA1(:,l(i,:));
end


% -------------------------------------------------------------------------
% Set more parameters
% -------------------------------------------------------------------------
Nsamples  = Lw_p;
TimeOnset = 0; % in sec
chlabels = cell(Nchannels,1);
for nel = 1:Nchannels; chlabels{nel} = strcat('Ch',num2str(el(nel))); end
Ic = [1:Nchannels];  % number of channels

% -------------------------------------------------------------------------
% Initialize SPM
% -------------------------------------------------------------------------
spm('defaults','EEG');

% -------------------------------------------------------------------------
% create the time axis (should be the same for all trials)
%--------------------------------------------------------------------------
timeaxis = [0:(Nsamples-1)]./Fsample + TimeOnset;
ftdata = [];

for i=1:Ns
    
    data2(Ic,:,:) = squeeze(DATA2(i,:,:,:));
    % -------------------------------------------------------------------------
    %% Create the Fieldtrip raw struct
    %--------------------------------------------------------------------------
    for j = 1:Ndays
        ftdata.trial{j + (i-1)*Ndays} = squeeze(data2(:, :, j));
        ftdata.time{j + (i-1)*Ndays} = timeaxis;
    end
    clear data2
end
ftdata.fsample = Fsample;
ftdata.label   = chlabels;
ftdata.label   = ftdata.label(:);

% -------------------------------------------------------------------------
% Convert the ftdata struct to SPM M\EEG dataset
% -------------------------------------------------------------------------
D = spm_eeg_ft2spm(ftdata, fname);
D = type(D, 'single');                                                  % Sets the dataset type
D = chantype(D, ':', 'LFP');                                            % Sets the channel type

%%%%%%%
if binned_trials == 1
    Nconditions = n_bins;
    condition_index = floor(linspace(1, Ns, Nconditions+1));
    k = 0;
    for c = 1:Nconditions
        if c < Nconditions
            upperboundary = condition_index(c+1)-1;
        elseif c == Nconditions
            upperboundary = condition_index(c+1);
        end
        for i = 1:length(condition_index(c):upperboundary)
            k = k + 1;
            D =  conditions(D, k, ['Condition ',num2str(c)]);
        end
    end
else
    switch phase
        case 'search_repeat'
            n_total = [n_S, n_R];
            for p = 1:2
                D =  conditions(D, [1:n_total(p)]+(p-1)*n_total(1), ['Condition ',num2str(p)]);
            end
        otherwise
            for i=1:Ns
                D =  conditions(D, [1:Ndays] + (i-1).*Ndays, ['Condition ',num2str(i)]);  % Sets the condition label
            end
    end
end


% -------------------------------------------------------------------------
% save (necessary step!)
save(D);
% -------------------------------------------------------------------------

%% step 2: add trial information
load(fname)

if binned_trials == 1
    Nconditions = n_bins;
    condition_index = round(linspace(1, Ns, Nconditions+1));
    k = 0;
    for c = 1:Nconditions
        if c < Nconditions
            upperboundary = condition_index(c+1)-1;
        elseif c == Nconditions
            upperboundary = condition_index(c+1);
        end
        for i = 1:length(condition_index(c):upperboundary)
            k = k + 1;
            D.condlist{1,k} = ['Condition ',num2str(c)];
            % [i, c, condition_index(c), condition_index(c+1)]
        end
    end
    Ns = Nconditions;
else
    switch phase
        case 'search_repeat'
            for i = 1:n_S
                D.condlist{1,i} = ['Condition ',num2str(1)];
            end
            for i = 1:n_R
                D.condlist{1,i+n_S} = ['Condition ',num2str(2)];
            end
        otherwise
            for i=1:Ns
                D.condlist{1,i} = ['Condition ',num2str(i)];
            end
    end
end
for i=1:size(D.trials,2)
    D.trials(i).onset = [];
end

save(fname,'D')
cd(home_dir)

