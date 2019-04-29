function [PEB3, PEB_all3, BMA] = humanPST_peb3(subs, runs, electrodes, n_bins, spectra, model_architecture)
    %% create structure
    PEB_all2 = cell(2,1);
    for s = subs
        [PEB, PEB_all, BMA] = run_peb_of_peb_humanPST(s, runs, electrodes, n_bins, spectra, model_architecture);

        PEB_all2{s} = PEB;
        clear PEB
    end
    [PEB3 , PEB_all3]    = spm_dcm_peb(PEB_all2);
    
    % Bayesian Model Averaging
    BMA = spm_dcm_peb_bmc(PEB3);

end