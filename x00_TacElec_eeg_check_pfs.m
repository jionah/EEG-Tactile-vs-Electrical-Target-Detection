%% Plot PFs
% Plot psychometric functions and apply data exclusion criteria (min P(detected) > 10% || max P(not detected) < 90%)

clear
close all

mydir    = '...\...\TacElec_EEG'; % Path to project folder

groups   = {'Tactile', 'Electrical'};
ana_dir  = fullfile(mydir,'analysis');
trg_dir = fullfile(mydir,'2nd level','Behaviour');

nRuns   =  8;
nInts   = 10;
nTrials = 200;

%% 1. Load data
% =========================================================================
for g = 1:numel(groups)
    group = groups{g};

    if strcmp(group,'Tactile')
        % 29 SJs (4 excluded; logs of excluded VP26 [1:3 5:7] lost in HDD crash
        % after publication, can mostly be recovered from raw data
        SJs  = { 'VP01' 'VP02' 'VP03' 'VP04' 'VP05' 'VP06'      'VP07' 'VP08' 'VP09' 'VP10' 'VP11' 'VP12' 'VP13' 'VP14'   'VP15'  'VP16' 'VP17' 'VP18'    'VP19'     'VP20'    'VP21'  'VP22'   'VP24' 'VP25'   'VP27'  'VP28'  'VP29'  'VP30'};
        runs = {[1 3:7]  2:7   2:8  [1:5 7] [1:5 7] [1:2 4:8]   1:7    1:7   [1:5 8]  1:7    1:7    2:7    1:6 [1:4 6:7] [1:3 5:7] 1:7    1:6 [1:4 6:7] [1:2 4 6:8] [1:4 7:8] [1 3:8] [1:2 4:8]  1:6    3:8   [1 3:8]  [1 3:8]  2:7  [1 3:5 7:8]};
        meancol = [0.49, 0.73, 0.91]; % Maya blue
    elseif strcmp(group,'Electrical')
        % 28 SJs (3 excluded; excluded S05 not plotted because intensities were changed after runs 1:2)
        SJs  = {    'S01'  'S02' 'S03' 'S06'   'S07' 'S08' 'S10' 'S12'   'S13'  'S16' 'S18'    'S19'   'S20' 'S21' 'S22' 'S25'  'S26'  'S27' 'S28' 'S29' 'S30' 'S31' 'S33' 'S35' 'S36' 'S37' 'S38'};
        runs = { [1:4 6:7]  1:6   1:7 [1:2 4:8] 1:7   1:7   1:7  [1:6 8] 1:7    2:8  [1:5 7] [1 3 5:8]  1:7   2:8   1:7   1:7 [1:4 6:7] 1:7   1:7   1:7   1:7  1:7    1:7   1:7   1:7   1:7   1:6};
        meancol = [0.89, 0.15, 0.21]; % Alizarin crimson
    end

    nSub    = numel(SJs);

    cd(ana_dir)
    [Data, PFs] = load_logs(SJs,nRuns,runs,mydir,ana_dir,group);

    %% 2. Get psychometric functions
    % =========================================================================

    logistic = @(c,x) (1./(1+exp(-c(2)*(x-c(1)))));

    % Normalise PFs to intensity range 1-10
    iniT = 5.5;
    ints = cell(nSub,nRuns);
    resp = cell(nSub,nRuns);
    norm_PFs = cell(nSub,nRuns);
    norm_slopes = nan(nSub,nRuns);
    norm_threshs = nan(nSub,nRuns);

    for s = 1:nSub
        disp(SJs{s})
        for r = runs{s}
            if strcmp(group,'Tactile')
                [~,ints{s,r}] = ismember(Data{s,r}.behaviour.PF.Intensities,Data{s,r}.Exp.Intensities); % Get Intensity levels
            elseif strcmp(group,'Electrical')
                [~,ints{s,r}] = ismember(Data{s,r}.behaviour.PF.ISIs,Data{s,r}.Exp.Intensities); % Get Intensity levels (field was mistakenly called "ISIs" in Electrical group)
            end
            resp{s,r} = Data{s,r}.behaviour.PF.Responses;                                           % Get responses
            norm_PFs{s,r} = fit_logistic(iniT,ints{s,r},resp{s,r},SJs{s},0);                        % Fit logistic function --> normalised PF
            norm_slopes(s,r) = norm_PFs{s,r}.fit_logistic(1,2);                                     % Get normalised slope
            norm_threshs(s,r) = norm_PFs{s,r}.fit_logistic(1,1);                                    % Get normalised T50
        end
    end

    mean_norm_slopes = mean(norm_slopes,2,'omitnan');
    mean_norm_threshs = mean(norm_threshs,2,'omitnan');

    normPFs.norm_PFs = norm_PFs;
    normPFs.norm_slopes = norm_slopes;
    normPFs.norm_threshs = norm_threshs;
    normPFs.mean_norm_slopes = mean_norm_slopes;
    normPFs.mean_norm_threshs = mean_norm_threshs;

    % Fit mean PF for every participant
    plot_range = 1:nInts;
    mean_fit = [];
    for s = 1:nSub
        mean_fit(s,:) = logistic([mean_norm_threshs(s), mean_norm_slopes(s)],plot_range);
    end

    normPFs.mean_fit = mean_fit;

    % Apply exclusion criterion: mean fitted det probability min > 10% or max < 90%
    excl = [];
    excl = round(mean_fit(:,1)*100) > 10 | round(mean_fit(:,end)*100) < 90;

    %% 3. Plot PFs
    fig_pos = [1,2,14,14];
    plot_range = 0:0.01:nInts;
    f = figure('Units','centimeters','position',fig_pos);
    hold on
    set(gca,'FontName','calibri','FontSize',12)

    area([0 10],[1 1],'FaceColor',[1 1 1], 'LineStyle', 'none')
    area([0 10],[.9 .9],'FaceColor',[.9 .9 .9], 'LineStyle', 'none')
    area([0 10],[.1 .1],'FaceColor',[1 1 1], 'LineStyle', 'none')
    
    for s = 1:nSub
        sub_pf = logistic([mean_norm_threshs(s), mean_norm_slopes(s)],plot_range);
        if excl(s)
            fprintf('%s min: %f max: %f\n',SJs{s}, mean_fit(s,1), mean_fit(s,10));
            line(plot_range,sub_pf,'Color',[.5 .5 .5],'LineWidth',1,'LineStyle','--');
        else
            line(plot_range,sub_pf,'Color',[0 0 0],'LineWidth',1);
        end
    end

    % Plot mean across subjects
    avg_pf = logistic([mean(mean_norm_threshs), mean(mean_norm_slopes)],plot_range);
    line(plot_range,avg_pf,'Color',meancol,'LineWidth',3);

    axis([1 10 0 1]);
    line([1 1],[0 1],'color',[0 0 0])
    line([1 10],[0 0],'color',[0 0 0])
    ylabel('P(detected)');
    xlabel('Intensity level')
    title(group,'fontsize',18)
    set(gca,...
        'xtick',1:10,...
        'ytick',[0,.5,1],...
        'ticklength',[.015 .025],...
        'tickdir','both',...
        'yminortick','on',...
        'fontsize',18)

    a = gca;
    a.XRuler.TickLabelGapOffset = -4;
    % a.YRuler.TickLabelGapOffset = -8;

    saveas(f,fullfile(trg_dir,[group '_normPFs.tiff']))

    % Store mean slopes of included participants for between-groups comparison
    mean_slopes{g} = mean_norm_slopes(~excl);
   

end
% Compare mean slopes between the two groups
bf10 = bf.ttest2( mean_slopes{1}, mean_slopes{2} );
fprintf('BF10 for a difference in slopes between the two groups: %.2f',bf10);
