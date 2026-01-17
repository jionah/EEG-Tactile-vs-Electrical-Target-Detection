%% Compute ERPs by intensity and detection

clear
close all


mydir    = '...\...\TacElec_EEG'; % Path to project folder

groups = {'Tac', 'Elec'};

log_f      = '_trial_log';
conditions = {'int', 'det'};

trg_dir = fullfile(mydir,'2nd level','ERP');
if ~exist(trg_dir,'dir')
    mkdir(trg_dir)
end

for g = 1:numel(groups)

    group = groups{g};
    disp(['Group: ' group])

    if strcmpi(group, 'Tac')
        SJs  = { 'VP02' 'VP03' 'VP04' 'VP05' 'VP06' 'VP07' 'VP08' 'VP09' 'VP10' 'VP11' 'VP12' 'VP15' 'VP16' 'VP17' 'VP18' 'VP19' 'VP20' 'VP21' 'VP22' 'VP24' 'VP25' 'VP27' 'VP28' 'VP29' 'VP30'};
    elseif strcmpi(group, 'Elec')
        SJs   = { 'S01' 'S02' 'S03' 'S06' 'S07' 'S08' 'S10' 'S12' 'S16' 'S19' 'S20' 'S21' 'S22' 'S25' 'S26' 'S27' 'S28' 'S29' 'S30' 'S31' 'S33' 'S35' 'S36' 'S37' 'S38'};
    end

    prefix     = ['bfraeTMdff' group '_'];

    % Load EEG data dummy
    D = spm_eeg_load(fullfile(mydir,'preprocessed',group,[prefix SJs{1} '.mat']));

    channels = D.indchantype('EEG');
    nChannels = length(channels);
    nSamples = D.nsamples;

    %% Get ERPs
    for c = 1:numel(conditions)
        condition = conditions{c};
        if strcmpi(condition,'int'), nConds = 10; elseif strcmpi(condition,'det'), nConds = 2; end
        ERP_file = fullfile(trg_dir,[group '_ERPs_' condition '.mat']);
        disp(['Computing grand mean ERPs: ' condition])

        ERP = nan(nChannels, nSamples, nConds, numel(SJs));

        for s = 1:numel(SJs)

            disp(SJs{s})

            % Get EEG data
            D = spm_eeg_load(fullfile(mydir,'preprocessed',group,[prefix SJs{s} '.mat']));

            % Get trial definitions
            trlog = fullfile(mydir,'logs',group, [SJs{s} log_f, '.mat']);
            load(trlog)

            % make sure trial numbers match
            if D.ntrials ~= size(trial_log.data,2)
                error(['Trial numbers do not match!!! ' SJs{s}])
            end

            % Get trial conditions
            [~,idx] = ismember(condition,trial_log.labels);
            conds = unique(trial_log.data(idx,:));

            for cc = conds
                filt = find(trial_log.data(idx,:) == cc);
                cond_data = D(channels,:,filt);
                if strcmp(condition,'int'), ERP(:,:,cc,s) = mean(cond_data,3);
                elseif strcmp(condition,'det'), ERP(:,:,cc+1,s) = mean(cond_data,3);
                end
            end
        end
        save(ERP_file,'ERP');
    end

end