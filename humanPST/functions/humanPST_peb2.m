function [PEB, PEB_all, BMA] = humanPST_peb2(s, runs, electrodes, n_bins, spectra, model_architecture)
%% PEB per subject

PEB_all = cell(length(runs),1);
for r = runs
    clear PEB
    [PEB, BMA] = run_peb_humanPST(s, r, electrodes, n_bins, spectra, model_architecture);
    PEB_all{r} = PEB;
end

show_F = NaN(length(PEB_all),1); for i = 1:length(PEB_all); show_F(i) = PEB_all{i,1}.F; end

%% Estimate model
try
    PEB    = spm_dcm_peb(PEB_all);
    % Bayesian Model Averaging
    BMA = spm_dcm_peb_bmc(PEB);
catch
    PEB_all2 = cell(size(PEB_all));
    k = 0;
    for i = 1:length(PEB_all)
        if PEB_all{i,1}.F>0
            k = k +1;
            PEB_all2{k} = PEB_all{i};
        end
    end
    PEB_all2 = PEB_all2(~cellfun('isempty',PEB_all2));
    [PEB , PEB_all2]    = spm_dcm_peb(PEB_all2);
    % Bayesian Model Averaging
    BMA = spm_dcm_peb_bmc(PEB);
end

end

