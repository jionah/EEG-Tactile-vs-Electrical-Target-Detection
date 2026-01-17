%% VBA Group BMS for Tac/Elec
% Perform between-groups BMS (tactile vs. electrical), first on the family 
% level, then only within the +family (int, det, pf), then for (int, det) 
% and (det, pf).

clear
clc

mydir    = '...\...\TacElec_EEG'; % Path to project folder

groups   = {'Tac', 'Elec'};

SJs{1}  = { 'VP02' 'VP03' 'VP04' 'VP05' 'VP06' 'VP07' 'VP08' 'VP09' 'VP10' 'VP11' 'VP12' 'VP15' 'VP16' 'VP17' 'VP18' 'VP19' 'VP20' 'VP21' 'VP22' 'VP24' 'VP25' 'VP27' 'VP28' 'VP29' 'VP30'};
SJs{2}  = { 'S01'  'S02' 'S03' 'S06' 'S07' 'S08' 'S10' 'S12' 'S16' 'S19' 'S20' 'S21' 'S22' 'S25' 'S26' 'S27' 'S28' 'S29' 'S30' 'S31' 'S33' 'S35' 'S36' 'S37' 'S38'};
subs    =  { 1:numel(SJs{1}) 1:numel(SJs{2}) };

modspaces = { {'null', 'int', 'det', 'pf', 'unc', 'rep', 'cue'},   {'int', 'det', 'pf'},    {'int', 'det'},    {'det', 'pf'}};
mstring   = { 'all', '+fam', 'int-det', 'det-pf'};
partition = {1 [2 3 4] 5 6 7}; % is used only for 'all'

%% Get dummy EEG file to retrieve channels and timing
results_file = fullfile(mydir,'bayesglm', groups{1}, SJs{1}{1}, ['bayesglm_results_' SJs{1}{1} '_' modspaces{1}{1}] );
load(results_file);
channels    = results.sensors;
time        = results.time;
nTimepoints = length(time);
nSensors    = length(channels);
nGroups     = numel(groups);

%% Group BMS for 3 model subspaces: (int-det-pf), (int-det), (det-pf)
for ms = 1:numel(modspaces)

    models  = modspaces{ms};
    nModels = numel(models);

    %% Assemble subject-level LMEs
    lmes = cell(1,2);
    disp('Retrieve data...')
    for g = 1:nGroups
        lme = zeros(numel(SJs{g}),nModels,nTimepoints,nSensors);
        for s = subs{g}
            sj = SJs{g}{s};
            for m = 1:nModels
                load(fullfile(mydir,'bayesglm',groups{g},sj,['bayesglm_results_' sj '_' models{m}]));
                lme(s,m,:,:) = results.LogEv;
            end
        end
        lmes{g} = lme;
    end
    disp('Done.')


    %% Perform BMS with and without pooling groups
    disp('Running Between-Groups BMS-RFX...')
    disp(['Model space: ' mstring{ms}])
    trg_dir = fullfile(mydir,'2nd level','BMS');
    if ~exist(trg_dir,'dir')
        mkdir(trg_dir)
    end

    L1      = permute(lmes{1}, [2 1 4 3]);
    L2      = permute(lmes{2}, [2 1 4 3]);
    [p_bg, bf10_bg] = deal(nan(nSensors,nTimepoints));

    % VBA parameters
    options = [];
    options.MinIter = 100;
    options.MaxIter = 1000;
    options.TolFun  = 1e-5;
    if strcmpi(mstring{ms},'all'), options.families = partition; end
    options.modelNames = models;
    options.DisplayWin = 0;
    options.verbose    = 0;

    for sensor = 1:nSensors

        fprintf('Sensor %d/%d ',sensor,nSensors)

        % Compute posterior model probabilities
        % nbytes = fprintf('Time point 0/%d ',nTimepoints);
        for t = 1:nTimepoints
        %     while nbytes > 0
        %         fprintf('\b');
        %         nbytes = nbytes - 1;
        %     end
            % nbytes = fprintf('Time point %d/%d',t,nTimepoints);

            [~, p] = VBA_groupBMC_btwGroups({L1(:,:,sensor,t) L2(:,:,sensor,t)}, options); % VBA internally computes the BF10 (exp(Fd-Fe)), but converts to posterior probability for output
            p_bg(sensor,t) = p;
            bf10_bg(sensor,t) = (1-p)/p; % convert back to BF10 (same as exp(Fd-Fe) in VBA_groupBMC_btwGroups(), line 32)        
        end

        % fprintf('\n')

        fprintf('Sensor %d/%d done!\n',sensor,nSensors)
    
    end
    % Save
    save(fullfile(trg_dir, sprintf('BMS_BetweenGroups_%s.mat',mstring{ms})), 'p_bg', 'bf10_bg');

end