function fname = humanPST_timeDCM_prep(subject, n_bins, run)

% -------------------------------------------------------------------------
% set path
% -------------------------------------------------------------------------
home_dir = pwd;
data_path = strcat(home_dir, '/data_preprocessed/');
cd(data_path)

data_name = strcat(subject, '_run', num2str(run),'_ft_preprocessed');
path_to_Dfiles = strcat(home_dir, '/data_preinversion/');      % change directory to folder with data

% -------------------------------------------------------------------------
% load data
% -------------------------------------------------------------------------
load(data_name)
ftdata = EEG_run_ft;
n_trials = length(ftdata.trial);

%--------------------------------------------------------------------------
% Initialize SPM
%--------------------------------------------------------------------------

    spm('defaults','EEG');

%--------------------------------------------------------------------------
% Convert the ftdata struct to SPM M\EEG dataset
%--------------------------------------------------------------------------
fname = [path_to_Dfiles, 'spm_', subject, '_bin', num2str(n_bins),'_run', num2str(run),'_indMRI_full'];                               % name of created spm file
D = spm_eeg_ft2spm(ftdata, fname);

%--------------------------------------------------------------------------
% Some details about the data
%--------------------------------------------------------------------------
D = units(D, D.indchantype('EEG'), 'uV');
Nconditions = n_bins;
condition_index = floor(linspace(1, n_trials, Nconditions+1));
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

%--------------------------------------------------------------------------
% Sensor position
%--------------------------------------------------------------------------
cd(strcat(home_dir, '/data_raw/'))
elec_name = strcat(subject,'Exported Electrodes');
elec = ft_read_sens(strcat(elec_name, '.sfp'));
cd(data_path)

% correct elec file
chanpos = elec.chanpos(4:end,:);
chantype = elec.chantype(4:end);
chanunit = elec.chanunit(4:end);
elecpos = elec.elecpos(4:end,:);
label = elec.label(4:end);
save(elec_name, 'chanpos', 'chantype', 'chanunit', 'elecpos', 'label')

% correct fiducials file
fid_name = strcat(subject,'_fiducials');
% % make fiducial file
chanpos = elec.chanpos(1:3,:);
chantype = elec.chantype(1:3);
chanunit = elec.chanunit(1:3);
elecpos = elec.chanpos;
label = elec.label(1:3);
type = elec.type;
unit = 'cm';
save(fid_name, 'chanpos', 'chantype', 'chanunit', 'elecpos', 'label')

S = [];
S.D = D;
S.task = 'loadeegsens';
S.source = 'mat';
S.sensfile = elec_name;
S.fidlabel = 'spmnas spmlpa spmrpa';
S.headshapefile = fid_name;
D = spm_eeg_prep(S);

% -------------------------------------------------------------------------
% mesh
% -------------------------------------------------------------------------
S = [];
S.D = D;
S.task = 'headshape';
S.source = 'mat';
S.headshapefile = fid_name;
S.fidlabel = 'nas lpa rpa';
D = spm_eeg_prep(S);
D.template = 1;
sMRI = strcat(home_dir, '/MRIs/', subject,'.img');
Msize = 2;
D.val = 1;
D.inv{D.val}.mesh = spm_eeg_inv_mesh(sMRI, Msize);
spm_eeg_inv_checkmeshes(D);
% -------------------------------------------------------------------------
% coreg: MRI
% -------------------------------------------------------------------------
D.template = 1;
D.useheadshape = 1;
D.sourcefid = D.fiducials;
D.targetfid  = D.inv{D.val}.mesh.fid;

M1= spm_eeg_inv_datareg(D);
D.inv{D.val}.datareg.sensors = ft_transform_sens(M1, D.sensors('EEG'));
D.inv{D.val}.datareg.fid_eeg = ft_transform_headshape(M1, D.sourcefid);
D.inv{D.val}.datareg.fid_mri = D.targetfid;
D.inv{D.val}.datareg.toMNI = D.inv{D.val}.mesh.Affine;
D.inv{D.val}.datareg.fromMNI = inv(D.inv{D.val}.datareg.toMNI);
D.inv{D.val}.datareg.modality = 'EEG';
%--------------------------------------------------------------------------
% forward model: headmodel
% %--------------------------------------------------------------------------
D.inv{D.val}.forward(1).voltype = 'EEG BEM';
D = spm_eeg_inv_forward(D,D.val);

% -------------------------------------------------------------------------
% save
% -------------------------------------------------------------------------
save(fname,'D')
cd(home_dir)


end