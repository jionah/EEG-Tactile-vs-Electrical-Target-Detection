%% Run Bayesian GLM estimation on electrode data

clear
close all

disp('Initialising...')
spm('defaults','eeg')
spm_jobman('initcfg');
spm_get_defaults('cmdline',true)

mydir    = '...\...\TacElec_EEG'; % Path to project folder

groups = {'Tac', 'Elec'};

%% Directories etc.
log_f = '_trial_log';
models = {'null', 'int', 'det', 'pf', 'unc', 'rep', 'cue'};

chantype = 'EEG';

src_folder  = 'preprocessed';
log_folder  = 'logs';
trg_folder  = 'bayesglm';

for g = 1:numel(groups)

    group = groups{g};
    disp(['Group: ' group])

    if strcmpi(group, 'Tac')
        SJs  = { 'VP02' 'VP03' 'VP04' 'VP05' 'VP06' 'VP07' 'VP08' 'VP09' 'VP10' 'VP11' 'VP12' 'VP15' 'VP16' 'VP17' 'VP18' 'VP19' 'VP20' 'VP21' 'VP22' 'VP24' 'VP25' 'VP27' 'VP28' 'VP29' 'VP30'};
    elseif strcmpi(group, 'Elec')
        SJs   = { 'S01' 'S02' 'S03' 'S06' 'S07' 'S08' 'S10' 'S12' 'S16' 'S19' 'S20' 'S21' 'S22' 'S25' 'S26' 'S27' 'S28' 'S29' 'S30' 'S31' 'S33' 'S35' 'S36' 'S37' 'S38'};
    end
    prefix = ['bfraeTMdff' group '_'];

    info_str = ['BayesGLM on all SJs\nModels: ' strjoin(models, ', ') '\n'];
    
    %% Loop through subjects and models
    for s = 1:numel(SJs)

        sj_eeg_file = fullfile(mydir, src_folder, group, [prefix,SJs{s},'.mat']);
        sj_trial_log = fullfile(mydir, log_folder, group, [SJs{s} log_f '.mat']);

        D = spm_eeg_load(sj_eeg_file);
        load(sj_trial_log);

        for m = 1:numel(models)

            fprintf('%s %s ', SJs{s}, models{m});

            % Get regressor
            if strcmp(models{m},'null')
                reg = [];
            else
                model_idx = strcmp(trial_log.labels, models{m});
                reg = zscore(trial_log.data(model_idx,:)');
            end

            S = [];
            S.D = D;
            S.chantype = chantype;
            S.reg = reg;
            S.alpha = 200;

            results = bayesglm_sensors(S);

            results.sj = SJs{s};
            results.model = models{m};
            results.sensors = D.chanlabels(D.indchantype('EEG'));
            results.time = D.time;

            sj_trg_dir = fullfile(mydir, trg_folder, group, SJs{s});
            if ~exist(sj_trg_dir,'dir')
                mkdir(sj_trg_dir)
            end
            save(fullfile(sj_trg_dir, ['bayesglm_results_' SJs{s} '_' models{m} '.mat']),'results')

        end

    end

end
disp('Done!')