function humanPST_timeDCM_inv(subject, fname, n_bins, run,electrodes, electrode_positions, model_architecture)     


%--------------------------------------------------------------------------
% Data filename
%--------------------------------------------------------------------------
home_dir = pwd;
D = spm_eeg_load(fname);
analysis_dir = strcat(home_dir, '/analysis_DCM/', electrodes, '/',model_architecture, '/');
mkdir(analysis_dir)
cd(analysis_dir)
DCM.xY.Dfile = fname;

%--------------------------------------------------------------------------
% Parameters and options used for setting up model
%--------------------------------------------------------------------------
%T2 = floor(1000.*D.nsamples./D.fsample);                                    % number of time samples/window
Nconditions = length(unique(D.condlist));
clear D

DCM.options.analysis    = 'CSD'; % analyze evoked responses
DCM.options.model       = 'CMC';
DCM.options.spatial     = 'IMG'; % spatial model
DCM.options.Tdcm(1)     = -1000;     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)     = 0;   % end of peri-stimulus time to be modelled
DCM.options.Nmodes      = 2;     %previous:2 nr of modes for data selection
DCM.options.h           = 1;     % nr of DCT components
DCM.options.Fdcm(1)     = 3;
DCM.options.Fdcm(2)     = 40;
DCM.options.D           = 1;     % downsampling
DCM.options.hann        = 1;     % hanning
%--------------------------------------------------------------------------
% Location priors for dipoles
%--------------------------------------------------------------------------
DCM.Sname = {'dlPFC','preSMA'};   
DCM.Lpos = electrode_positions;

%--------------------------------------------------------------------------
% Spatial model
%--------------------------------------------------------------------------
DCM.xY.Ic = [1:64]';
DCM.xY.modality = 'EEG';
DCM = spm_dcm_erp_dipfit(DCM);

%--------------------------------------------------------------------------
% Specify connectivity model
%--------------------------------------------------------------------------

DCM.A = make_model_architecture(1:size(DCM.Sname,2), model_architecture);
DCM.B = {};
%DCM.C = [0];

%--------------------------------------------------------------------------
% Between trial effects
%--------------------------------------------------------------------------
DCM.xU.X = [];
DCM.xU.name = {};

%--------------------------------------------------------------------------
% Save
%--------------------------------------------------------------------------
DCM.name = ['DCM_', subject,'_bin', num2str(n_bins), '_run', num2str(run),'_full'];                                           % name of new DCM structure
DCM_name = DCM.name;
save(DCM.name,'DCM')
%cd(home_dir)

%--------------------------------------------------------------------------
% inversion
%--------------------------------------------------------------------------
tic
for tr = 1:Nconditions
    load(DCM.name)
    DCM.options.trials = tr;
    DCM = spm_dcm_csd(DCM);
    save([DCM.name '_' num2str(tr)], 'DCM', spm_get_defaults('mat.format'))
end
toc
cd(home_dir)
end