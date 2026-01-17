%% Compute subsample ERPs and mean ERP over subsamples
clear
clc

mydir    = '...\...\TacElec_EEG'; % Path to project folder

ca_dir = fullfile(mydir, '2nd level', 'control analysis');

log_f    = '_trial_log';

for g = 1:numel(groups)

    group = groups{g};
    disp(['Group: ' group])

    if strcmpi(group,'Tac')
        SJs  = { 'VP02' 'VP03' 'VP04' 'VP05' 'VP06' 'VP07' 'VP08' 'VP09' 'VP10' 'VP11' 'VP12' 'VP15' 'VP16' 'VP17' 'VP18' 'VP19' 'VP20' 'VP21' 'VP22' 'VP24' 'VP25' 'VP27' 'VP28' 'VP29' 'VP30'};
        prefix    = 'bfraeTMdffTac_';
    elseif strcmpi(group,'Elec')
        SJs  = {  'S01'  'S02' 'S03'  'S06' 'S07' 'S08' 'S10' 'S12' 'S16' 'S19' 'S20' 'S21' 'S22' 'S25' 'S26' 'S27' 'S28' 'S29' 'S30' 'S31' 'S33' 'S35' 'S36' 'S37' 'S38'};
        prefix    = 'bfraeTMdffElec_';
    end

    trg_dir = fullfile(ca_dir,'Subsamples',group);
    if ~exist(trg_dir,'dir')
        mkdir(trg_dir)
    end

    % Load EEG data dummy
    D = spm_eeg_load(fullfile(mydir,'preprocessed',group, [prefix SJs{1},'.mat']));
    channels  = D.indchantype('EEG');
    nChannels = length(channels);
    nSamples  = D.nsamples;
    nConds    = 2; % hit/miss

    % Load subsamples file
    load(fullfile(ca_dir,'Subsamples.mat'))
    nSubsamples = size(Subsamples.(group).(SJs{1}).trial_sample,3);

    %% Get ERPs
    ERP = nan(nChannels, nSamples, nConds, numel(SJs));
    
    nbytes = fprintf('Computing ERPs for subsample 0/%d ',nSubsamples);
    for sss = 1:nSubsamples
        while nbytes > 0
            fprintf('\b');
            nbytes = nbytes - 1;
        end
        nbytes = fprintf('Computing grand mean ERPs for subsample %d/%d\n',sss,nSubsamples);

        for s = 1:numel(SJs)
            sj = SJs{s};

            % Get EEG data
            D = spm_eeg_load(fullfile(mydir,'preprocessed',group,[prefix sj '.mat']));

            % Get trial definitions
            trlog = fullfile(mydir, 'logs', group, [sj log_f  '.mat']);
            load(trlog)

            % make sure trial numbers match
            if D.ntrials ~= size(trial_log.data,2)
                disp(['Trial numbers do not match!!! ' sj])
                keyboard
            end

            % Get trial conditions
            [~, idx] = ismember('det', trial_log.labels);
            conds = unique(trial_log.data(idx,:));

            trial_samples = Subsamples.(group).(sj).trial_sample(:,:,sss);

            miss_idcs = sort([trial_samples{1,:}]); % miss
            hit_idcs = sort([trial_samples{2,:}]); % hit

            cond_data_miss = D(channels,:,miss_idcs);
            cond_data_hit = D(channels,:,hit_idcs);
            ERP(:,:,1,s) = mean(cond_data_miss,3,'omitnan');
            ERP(:,:,2,s) = mean(cond_data_hit,3,'omitnan');

        end % of subjects loop

        % Save
        
        ERP_file = fullfile(trg_dir,['ERPs_det_intmatched_' num2str(sss) '.mat']);
        save(ERP_file,'ERP');

    end % of subsamples loop

    %% Average over subsamples
    fprintf('Computing average over %d subsamples...\n',nSubsamples)
    ERP_ssavg = nan(nChannels, nSamples, nConds, numel(SJs), nSubsamples);
    for sss = 1:nSubsamples
        load(fullfile(ca_dir,'Subsamples',group,['ERPs_det_intmatched_' num2str(sss)]));
        ERP_ssavg(:,:,:,:,sss) = ERP;
    end
    ERP = mean(ERP_ssavg,5);

    % Save
    ERP_file = fullfile(ca_dir,sprintf('%s_ERPs_det_intmatched_avg%d.mat',group,nSubsamples));
    save(ERP_file,'ERP');

end

disp('Done!')