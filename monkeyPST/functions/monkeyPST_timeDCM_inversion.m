function monkeyPST_timeDCM_inversion(el, fname,  DCM, analysis_dir, model_architecture)
%% DCM analysis: inversion

% el = 1 or 2 (preSMA or M1)


%% step 3: run batch
home_dir = pwd;
cd(analysis_dir)
load(fname)
%--------------------------------------------------------------------------
% Parameters and options used for setting up model
%--------------------------------------------------------------------------
Nareas = size(el, 2);
NConditions = length(unique(D.condlist));
DCM.xY.Dfile = fname;
DCM.options.analysis = 'CSD'; % analyze c
DCM.options.model    = 'CMC';
DCM.options.spatial  = 'LFP'; % spatial model
DCM.options.Tdcm(1)  = 0;     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = 1000;   % end of peri-stimulus time to be modelled
DCM.options.Nmodes   = Nareas;     % nr of modes for data selection
DCM.options.h        = 1;     % nr of DCT components
DCM.options.D        = 1;     % downsampling
DCM.options.Fdcm(1) = 3;
DCM.options.Fdcm(2) = 40;
DCM.options.hann     = 1;     % hanning

% Location priors for dipoles
%--------------------------------------------------------------------------
DCM.Sname = cell(Nareas,1);
for nel = 1:Nareas; DCM.Sname{nel} = strcat('LFP',num2str(nel)); end


% Spatial model
        %--------------------------------------------------------------------------
        % Specify connectivity model
        %--------------------------------------------------------------------------
        DCM.A = make_model_architecture(el, model_architecture);
        DCM.B = {};
        
        %--------------------------------------------------------------------------
        % Between trial effects
        %--------------------------------------------------------------------------
        DCM.xU.X = [];
        DCM.xU.name = {};

        % Save
        save(DCM.name,'DCM')
        
        
        %% Step 3: inversion
        tic
        for tr = 1:NConditions
            load(DCM.name)
            DCM.options.trials = tr;
            DCM = spm_dcm_csd(DCM);
            save([DCM.name '_' num2str(tr)], 'DCM', spm_get_defaults('mat.format'))
        end
        toc
        cd(home_dir)
        clearvars -except d m el electrodes


end

