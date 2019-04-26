function [spm_fname, DCM]= monkeyPST_timeDCM_prep(m, el, d, electrodes, phase, spatial_ref, binned_trials, n_bins, clip, data_path, analysis_dir)

% -------------------------------------------------------------------------
% Set parameters
% -------------------------------------------------------------------------

home_dir = pwd;
monkey_name = {'Risette', 'Satan'};
monkey_init = {'r', 's'};
Fsample = 781;  %Hz                                                               % sampling frequency
Lw = 1000; % in ms
TimeOnset = 0; % in sec

% -------------------------------------------------------------------------
% set unit index and n_trials
% -------------------------------------------------------------------------
if binned_trials == 1
    n_trials = n_bins;
else
    n_trials = monkeyPST_make_n_trials(phase);
end

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
                
            case 'left_frontal'
                if m == 1
                    electrode_index = [4 7];
                elseif m == 2
                    electrode_index = [6 12];
                end
            case 'right_frontal'
                if m == 1
                    electrode_index = [11 14];
                elseif m == 2
                    electrode_index = [9 13];
                end                
            otherwise
                fprintf('You have not specified the electrode location!')
        end
end

% -------------------------------------------------------------------------
% load data and name files
% -------------------------------------------------------------------------
cd(data_path)
raw_fname = [monkey_name{m}, '_c', num2str(clip), '_d_' num2str(d)];
load(raw_fname)
switch phase
    case 'search'
        data_phase = data_S;
    case 'repeat'
        data_phase = data_R ;
    case 'search_repeat'
        n_S = length(data_S)/Fsample;
        n_R = length(data_R)/Fsample;
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
spm_fname = [analysis_dir 'spm_', monkey_init{m}, '_d_', num2str(d), '_c', num2str(clip),'_LFP_full'];                               % name of created spm file
cd(analysis_dir);
% -------------------------------------------------------------------------
% select electrode
% -------------------------------------------------------------------------
DATA2 = data1(el,:);
DATA1 = DATA2;
clear DATA2
clear data1

% -------------------------------------------------------------------------
% Time information
% -------------------------------------------------------------------------
[Nchannels,Nsamples,Ndays] = size(DATA1); % here trials is the number of days (or seizures etc)                                                                % length of the window in ms
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
D = spm_eeg_ft2spm(ftdata, spm_fname);
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
load(spm_fname)

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

address = ['DCM_', monkey_init{m}, '_', phase, '_d_', num2str(d),'_full'];
DCM.name = address;
save(spm_fname,'D')
cd(home_dir)

end