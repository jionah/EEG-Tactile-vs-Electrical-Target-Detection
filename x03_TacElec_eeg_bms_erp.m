%% Sensor level BMS
% On sensor level log evidence (from bayesglm_sensorlevel.m)

clear
close all

mydir    = '...\...\TacElec_EEG'; % Path to project folder

groups = {'Tac', 'Elec'};

models = {'null', 'int', 'det', 'pf', 'unc', 'rep', 'cue'};
partition = {1 [2 3 4] 5 6 7};
num_fams = numel(partition); % number of model families
mod_names = strjoin(models, '-');

glm_dir = 'bayesglm';
trg_dir = fullfile(mydir, '2nd level', 'BMS', mod_names);

if ~exist(trg_dir, 'dir')
    mkdir(trg_dir)
end

disp('BMS using VBA toolbox...')

for g = 1:numel(groups)

    group = groups{g};
    disp(['Group: ' group])
    
    if strcmpi(group, 'Tac')
        SJs  = { 'VP02' 'VP03' 'VP04' 'VP05' 'VP06' 'VP07' 'VP08' 'VP09' 'VP10' 'VP11' 'VP12' 'VP15' 'VP16' 'VP17' 'VP18' 'VP19' 'VP20' 'VP21' 'VP22' 'VP24' 'VP25' 'VP27' 'VP28' 'VP29' 'VP30'};
    elseif strcmpi(group, 'Elec')
        SJs   = { 'S01' 'S02' 'S03' 'S06' 'S07' 'S08' 'S10' 'S12' 'S16' 'S19' 'S20' 'S21' 'S22' 'S25' 'S26' 'S27' 'S28' 'S29' 'S30' 'S31' 'S33' 'S35' 'S36' 'S37' 'S38'};
    end

    %% Get dummy EEG file to retrieve channels and timing
    results_file = fullfile(mydir, glm_dir, group, SJs{1}, sprintf('bayesglm_results_%s_%s.mat', SJs{1}, models{1}) );
    load(results_file);
    channels = results.sensors;
    time = results.time;
    nTimepoints = length(time);
    nSensors = length(channels);

    %% Assemble data

    disp('----------------------- BMS -----------------------')
    disp(mod_names)

    disp('Retrieving data...')
    lme = nan(numel(SJs),numel(models),nTimepoints,nSensors);

    for s = 1:numel(SJs)
        for m = 1:numel(models)

            load(fullfile(mydir,glm_dir,group,SJs{s},sprintf('bayesglm_results_%s_%s.mat', SJs{s}, models{m})));

            lme(s,m,:,:) = results.LogEv;

            clear results
        end
    end
    fprintf('Done!\n')

    %% Run BMS per sensor and time point

    fprintf('Running BMS\n')

    exp_r       = nan(numel(models), nTimepoints, nSensors);
    xp          = nan(numel(models), nTimepoints, nSensors);
    exp_r_fam   = nan(num_fams, nTimepoints, nSensors);
    xp_fam      = nan(num_fams, nTimepoints, nSensors);

    % Set up BMS structure
    BMS.models      = models;
    BMS.SJs         = SJs;
    BMS.channels    = channels;
    BMS.time        = time;
    BMS.partition   = partition;
    BMS.exp_r       = exp_r;
    BMS.xp          = xp;
    BMS.exp_r_fam   = exp_r_fam;
    BMS.xp_fam      = xp_fam;

    for sensor = 1:nSensors

        fprintf('Sensor %d/%d ',sensor,nSensors)

        % Compute posterior model probabilities
        nbytes = fprintf('Time point 0/%d ',nTimepoints);
        for t = 1:nTimepoints
            while nbytes > 0
                fprintf('\b');
                nbytes = nbytes - 1;
            end
            nbytes = fprintf('Time point %d/%d',t,nTimepoints);

            % Prepare inputs to VBA_groupBMC() (L, options)
            L = permute(lme(:,:,t,sensor), [2 1]); % - L: Kxn array of log-model evidences (K models; n subjects)
            options.families = partition;
            options.modelNames = models;
            options.DisplayWin = 0;
            options.verbose = 0;

            % Extract results
            [~, out] = VBA_groupBMC(L,options);
            exp_r(:,t,sensor)     = out.Ef; % expectation of the posterior p(r|y)
            exp_r_fam(:,t,sensor) = out.families.Ef;
            xp(:,t,sensor)        = out.ep; % Exceedance probability
            xp_fam(:,t,sensor)    = out.families.ep;

        end

        % Update BMS structure and save
        BMS.exp_r       = exp_r;
        BMS.xp          = xp;
        BMS.exp_r_fam   = exp_r_fam;
        BMS.xp_fam      = xp_fam;

        
        fprintf('\n')
        fprintf('Sensor %d/%d done!\n',sensor,nSensors)
    end
    % Save
    save(fullfile(trg_dir, [group '_BMS_FamXPs.mat']),'BMS');
end

