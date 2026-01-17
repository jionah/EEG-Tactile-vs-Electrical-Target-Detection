%% Subsample ERPs to get intensity-matched detection ERP subsamples

clear
clc

mydir    = '...\...\TacElec_EEG'; % Path to project folder

groups   = {'Tac', 'Elec'}; % 'Tac' 'Elec'

rng(1234) % 1234 to reproduce the published subsample
n_subsamples = 40; % Number of subsamples to draw

trg_dir = fullfile(mydir,'2nd level','control analysis');

for g = 1:numel(groups)

    group = groups{g};
    disp('Drawing subsamples...')
    disp(['Group: ' group])
    
    if strcmpi(group,'Tac')
        SJs  = { 'VP02' 'VP03' 'VP04' 'VP05' 'VP06' 'VP07' 'VP08' 'VP09' 'VP10' 'VP11' 'VP12' 'VP15' 'VP16' 'VP17' 'VP18' 'VP19' 'VP20' 'VP21' 'VP22' 'VP24' 'VP25' 'VP27' 'VP28' 'VP29' 'VP30'};
        prefix = 'bfraeTMdffTac_';
    elseif strcmpi(group,'Elec')
        SJs   = { 'S01'  'S02' 'S03'  'S06' 'S07' 'S08' 'S10' 'S12' 'S16' 'S19' 'S20' 'S21' 'S22' 'S25' 'S26' 'S27' 'S28' 'S29' 'S30' 'S31' 'S33' 'S35' 'S36' 'S37' 'S38'};
        prefix = 'bfraeTMdffElec_';
    end

    log_f = '_trial_log';

    n_det = 2;
    n_int = 10;

    % Load subsamples file if existing
    filename = fullfile(mydir,'2nd level','control analysis','Subsamples.mat' );
    if exist(filename,'file')
        load(filename)
    end

    % Load EEG data dummy
    D = spm_eeg_load(fullfile(mydir,'preprocessed',group,[prefix SJs{1} '.mat']));

    channels  = D.indchantype('EEG');
    nChannels = length(channels);
    nSamples  = D.nsamples;

    %% Compute subsamples

    ntrials = nan(numel(SJs),1);
    for s = 1:numel(SJs)
        disp(SJs{s})

        % Get EEG data
        D = spm_eeg_load(fullfile(mydir,'preprocessed',group,[prefix SJs{s} '.mat']));

        % Get trial definitions
        trlog = fullfile(mydir, 'logs', group, [SJs{s} log_f, '.mat']);
        load(trlog)

        % make sure trial numbers match
        if D.ntrials ~= size(trial_log.data,2)
            error(['Trial numbers do not match!!! ' SJs{s}])
        end

        % Get detection and intensity conditions
        [~,int_idx] = ismember('int',trial_log.labels);
        [~,det_idx] = ismember('det',trial_log.labels);
        int_conds = unique(trial_log.data(int_idx,:));
        det_conds = unique(trial_log.data(det_idx,:));

        % Get int x det trial numbers
        n_intdet = nan(length(det_conds),length(int_conds));
        for d = 1:length(det_conds)
            for i = 1:length(int_conds)
                n_intdet(d,i) = sum(trial_log.data(det_idx,:) == det_conds(d) & trial_log.data(int_idx,:) == int_conds(i));
            end
        end
        min_trials = min(n_intdet);
        ntrials(s) = sum(min_trials);

        % Sample trials
        for it = 1:n_subsamples
            for i = 1:length(int_conds)
                for d = 1:n_det
                    filt = find(trial_log.data(det_idx,:) == det_conds(d) & trial_log.data(int_idx,:) == int_conds(i));
                    trial_sample = sort(randsample(n_intdet(d,i),min_trials(i)));
                    trial_idx = filt(trial_sample);
                    Subsamples.(group).(SJs{s}).trial_sample{d,i,it} = trial_idx;
                end
            end
        end
        Subsamples.(group).(SJs{s}).ntrials = ntrials(s);

    end

    % Save subsamples
    if ~exist(trg_dir,'dir')
        mkdir(trg_dir)
    end
    disp('Saving...')
    save(filename,'Subsamples');
    disp('Done!')
    
end % of groups loop

rng('shuffle') % Set back to random