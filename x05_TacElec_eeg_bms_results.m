%% Plot ERPs along with BMS results and beta estimates

clear
close all

mydir = '...\...\TacElec_EEG';  % Directory containing project folder

groups = {'Tac', 'Elec'};

xp_threshold = 0.99;
bf_threshold = 150;

%% Directories etc.

xticks = 0:100:600;
xuselims = [-50 600];

fig_pos = [1,2,5.6,7];
rect_width = 1.5;
ax_width = 1;
erp_width = 1;

cond = 'int';
ncond = 10;
chansel = { 'CP4' 'C6' 'CPz' };

% For supplement figure:
% chansel = { 'C2' };

bms_dir = fullfile(mydir,'2nd level','BMS','null-int-det-pf-unc-rep-cue');
erp_dir = fullfile(mydir,'2nd level','ERP');
glm_dir = fullfile(mydir,'bayesglm');

recompute_betas = 0;
recompute_betatest = 0;

for g = 1:numel(groups)

    group = groups{g};
    disp(['Group: ' group])

    %% Load data
    ERP_file = fullfile(erp_dir, [group '_ERPs_' cond '.mat']);
    load(ERP_file)

    mean_ERP    = mean(ERP,4);
    sd_ERP      = std(ERP,[],4);
    se_ERP      = std(ERP,[],4)/sqrt(size(ERP,4));

    BMS_file = fullfile(bms_dir, [group '_BMS_FamXPs.mat']);
    load(BMS_file)

    SJs         = BMS.SJs;
    models      = BMS.models;
    channels    = BMS.channels;
    time        = BMS.time*1000;
    plot_time   = 1:length(time);
    time        = time(plot_time);
    nChannels   = length(channels);
    nSamples    = length(time);
    nModels     = numel(models);
    nSubjects   = numel(SJs);
    nConds      = size(ERP,3);

    %% Get betas

    beta_file = fullfile(bms_dir, [group '_betas.mat']);

    if exist(beta_file,'file') && ~recompute_betas
        fprintf('\nLoading betas...')
        load(beta_file);
        fprintf('Done!\n')
    else
        betas = nan(nSubjects, nSamples, nModels, nChannels);

        for s = 1:nSubjects

            fprintf('Retrieving betas: %s ', SJs{s})

            sj_glm_dir = fullfile(glm_dir,group,SJs{s});

            for m = 1:nModels

                fprintf('%s ', models{m})

                load(fullfile(sj_glm_dir, ['bayesglm_results_' SJs{s} '_' models{m} '.mat']));

                if ~strcmp(models{m},'null')
                    betas(s,:,m,:) = permute(squeeze(results.beta(2,:,plot_time)), [3,2,4,1]);
                else
                    betas(s,:,m,:) = inf(1,nSamples,1,nChannels);
                end

            end

            fprintf('Done!\n')
        end
        save(beta_file,'betas');
    end

    beta_mean = permute(mean(betas,1),[3,2,4,1]);
    beta_mean(end,:,:) = -1*beta_mean(end,:,:);
    beta_std = permute(std(betas,[]),[3,2,4,1]);
    beta_se = beta_std./sqrt(nSubjects);

    %% Compute Bayesian t-tests on betas

    beta_test_file = fullfile(bms_dir, [group '_beta_test.mat']);

    if exist(beta_test_file,'file') && ~recompute_betatest
        fprintf('\nLoading beta test...')
        load(beta_test_file);
        bf10 = beta_test.bf10;
        pValue = beta_test.pValue;
        fprintf('Done!\n')
    else
        fprintf('\nComputing t-tests...\n')
        bf10 = nan(size(beta_mean));
        pValue = nan(size(beta_mean));

        for c = 1:nChannels
            fprintf('Channel %d: %s',c,channels{c})
            for m = 1:numel(models)
                fprintf(' %s',models{m})
                for s = 1:nSamples
                    if betas(:,s,m,c) == Inf
                        bf10(m,s,c) = Inf;
                        pValue(m,s,c) = Inf;
                    else
                        [bf10(m,s,c),pValue(m,s,c)] = bf_ttest(betas(:,s,m,c));
                    end
                end
            end
            fprintf('\n')
        end
        beta_test.bf10 = bf10;
        beta_test.pValue = pValue;
        save(beta_test_file,'beta_test');
    end

    %% Switch det and pf for plotting

    if ismember('pf',models)

        % models
        tmp = models;
        tmp(3) = models(4);
        tmp(4) = models(3);
        models = tmp;

        % betas
        tmp = beta_mean;
        tmp(3,:,:) = beta_mean(4,:,:);
        tmp(4,:,:) = beta_mean(3,:,:);
        beta_mean = tmp;

        % beta test
        tmp = bf10;
        tmp(3,:,:) = bf10(4,:,:);
        tmp(4,:,:) = bf10(3,:,:);
        bf10 = tmp;

        % XPs
        tmp = BMS.xp;
        tmp(3,:,:) = BMS.xp(4,:,:);
        tmp(4,:,:) = BMS.xp(3,:,:);
        BMS.xp = tmp;

    end

    %% Prepare plots

    erp_cols = parula(nConds+1);
    erp_cols = erp_cols(1:nConds,:);

    if strcmp(cond,'det')
        erp_cols = [0.89, 0.15, 0.21
            0.89, 0.15, 0.21];   % Alizarin crimson
        linestyle = {'--','-'};
    end

    model_cols  = [ 0  1  0   % int - green
                    0  0  1   % pf - blue
                    1  0  0   % det - red
                  ];

    fam_cols    = [.5 .5 .5   % null - grey
                   .8 .8 .8   % +fam - light grey
                    0  1  1   % unc - cyan
                    1  0  1   % rep - magenta
                    1  1  0   % cue - yellow
                  ];

    dot_cols = [ fam_cols(1,:)
        model_cols
        fam_cols(3:end,:)
        ];

    % Family info
    partition   = [1 2 2 2 3 4 5]; % BMS.partition
    nFam        = length(unique(partition));
    fam_idx     = cell(1,nFam);
    fam_size    = nan(1,nFam);
    for i = 1:nFam
        fam_idx{i} = find(partition == i);
        fam_size(i) = length(fam_idx{i});
    end

    %% Set figure positions

    offset = .1;
    width = .8;
    height1 = .44;
    height2 = .18;
    gap = .02;

    %% Plot -50ms to 600ms

    if isempty(chansel)
        cidx = 1:64;
    else
        [~,cidx] = ismember(chansel,channels);
    end

    for c = cidx

        f = figure('units','centimeters','position',fig_pos);

        chan = channels{c};

        %-- 1. Plot ERPs and BMS ------------------------------------------------
        subplot('position',[offset,offset+2*gap+2*height2,width,height1])
        hold on
        set(gca,'FontName','Calibri','FontSize',10)

        line([0 0],[-10 10],'color',[0 0 0],'linewidth',ax_width);

        h = nan(ncond,1);
        max_y = nan;
        min_y = nan;
        for l = 1:ncond
            x = time';
            y = squeeze(mean_ERP(c,plot_time,l))';
            se = squeeze(se_ERP(c,plot_time,l))';
            fill([x;flipud(x)],[y-se;flipud(y+se)],erp_cols(l,:),'linestyle','none','facealpha',.2);
            if strcmp(cond,'int')
                h(l) = line(x,y,'Color',erp_cols(l,:));
            else
                h(l) = line(x,y,'Color',erp_cols(l,:),'linestyle',linestyle{l});
            end
            max_y = max([max_y max(y+se)]);
            min_y = min([min_y min(y-se)]);
        end
        range_y = max_y - min_y;

        % Plot XPs of each family
        min_y   = min_y - (range_y/10);
        xlim(xuselims)
        ylim([min_y-range_y/7 max_y + range_y/10])
        for fam = 1:nFam
            % Plot winning models constrained by fam exceedance and beta test
            famXPs  = BMS.xp_fam(fam,plot_time,c)';
            exc     = famXPs >= xp_threshold;
            XPs     = BMS.xp(fam_idx{fam},plot_time,c)';
            for s = 1:length(time)
                if exc(s)
                    sbf10 = bf10(fam_idx{fam},plot_time(s),c);
                    if fam_size(fam) > 1
                        sXP = XPs(s,:);
                        [~,max_sXP] = max(sXP);
                        beta_test = sbf10(max_sXP) >= bf_threshold;
                        % plot model mix
                        if fam_size(fam) == 3
                            col = sXP([3 1 2]);
                        else
                            col = [sXP(2) sXP(1) 0];
                        end
                    else
                        col = fam_cols(fam,:);
                        beta_test = sbf10 >= bf_threshold;
                    end
                    if beta_test
                        try
                            plot(time([s s+1]),[min_y min_y],'Color',col,'LineWidth',7)
                        catch
                        end
                    end
                end
            end
        end

        set(gca,...
            'xtick',xticks,...
            'xticklabel',{},...
            'ticklength',[.015 .025],...
            'tickdir','both',...
            'xminortick','on')

        ylims = ylim;
        xlims = xlim;
        line([xlims(1) xlims(1)],[ylims(1) ylims(2)],'color',[0 0 0],'linewidth',ax_width);
        line([xlims(1) xlims(2)],[ylims(1) ylims(1)],'color',[0 0 0],'linewidth',ax_width);

        % Plot channel name
        text(10,ylims(2)-.5,chan,'fontweight','bold')

        %-- 2. Plot XP time courses ---------------------------------------------
        subplot('position',[offset,offset+gap+height2,width,height2])
        hold on
        set(gca,'FontName','Calibri','FontSize',10)

        % Get fam/model XPs
        famXPs = BMS.xp_fam(:,plot_time,c);
        XPs = BMS.xp(fam_idx{fam_size>1},plot_time,c);
        XPs = XPs.*repmat(famXPs(fam_size>1,:),fam_size(fam_size>1),1);  % Scale by family XP

        % Get corresponding colors
        data = zeros(nModels,nSamples,3);
        idx = 1;
        for ff = 1:nFam
            if fam_size(ff) > 1
                for m = 1:fam_size(ff)
                    col_mat = repmat(permute(model_cols(m,:),[1,3,2]),1,nSamples,1);
                    col_mat = col_mat.*repmat(XPs(m,:),1,1,3);
                    data(idx,:,:) = col_mat;
                    idx = idx+1;
                end
            else
                col_mat = repmat(permute(fam_cols(ff,:),[1,3,2]),1,nSamples,1);
                col_mat = col_mat.*repmat(famXPs(ff,:),1,1,3);
                data(idx,:,:) = col_mat;
                idx = idx+1;
            end
        end

        % plot as image
        image(time,1:nModels,data)
        set(gca,'YDir','reverse','ytick',1:nModels,'YTickLabel',{},'xtick',xticks,'XTickLabel',{},'tickdir','out');
        set(gca,...
            'YDir','reverse',...
            'ytick',1:nModels,...
            'YTickLabel',{},...
            'xtick',xticks,...
            'XTickLabel',{},...
            'ticklength',[.015 .025],...
            'tickdir','both',...
            'xminortick','on');
        ylims = ylim;
        line([0 0],ylims,'color',[1 1 1],'linewidth',ax_width)
        xlim(xuselims)
        ylim([.5 nModels+.5])

        % Mark threshold exceedance
        exc = famXPs >= xp_threshold;
        exc_diff = diff(exc,1,2);
        idx = .5;
        for ff = 1:nFam
            onsets = find(exc_diff(ff,:) == 1);
            offsets = find(exc_diff(ff,:) == -1) - 1;
            if ~isempty(onsets) && ~isempty(offsets)
                if offsets(1) == 0
                    offsets = offsets(2:end);
                end
                if offsets(1) < onsets(1)
                    onsets = [1 onsets];
                end
                if offsets(end) < onsets(end)
                    offsets = [offsets nSamples];
                end
            elseif ~isempty(onsets)
                offsets = nSamples;
            elseif ~isempty(offsets) && ~(offsets(1) == 0)
                onsets = 1;
            else
                idx = idx + fam_size(ff);
                continue
            end

            onsets = time(onsets);
            offsets = time(offsets);
            durations = offsets - onsets;
            durations(durations == 0) = 1;

            for r = 1:length(onsets)
                rectangle('Position',[onsets(r) idx durations(r) fam_size(ff)],'edgecolor',[1 1 1],'linewidth',rect_width)
            end

            idx = idx + fam_size(ff);
        end

        ylims = ylim;
        xlims = xlim;
        line([xlims(1) xlims(1)],[ylims(1) ylims(2)],'color',[0 0 0],'linewidth',ax_width);
        line([xlims(1) xlims(2)],[ylims(2) ylims(2)],'color',[0 0 0],'linewidth',ax_width);

        % Plot colored dots in y axis
        ylocs = get(gca, 'YTick');
        xloc = xlim;
        xText = xloc(1) - 30;  % Adjust as needed

        hold on;
        for i = 1:length(ylocs)
            text(xText, ylocs(i), '●', ...
                'Color', dot_cols(i,:), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 10);
        end
        hold off;


        %-- 3. plot betas -------------------------------------------------------
        subplot('position',[offset,offset,width,height2])
        hold on
        set(gca,'FontName','Calibri','FontSize',10)

        plot_betas = beta_mean(:,plot_time,c);
        plot_betas(1,:) = 0;

        imagesc(time,1:nModels,plot_betas);
        set(gca,...
            'YDir','reverse',...
            'ytick',1:nModels,...
            'YTickLabel',{},...
            'xtick',xticks,...
            'ticklength',[.015 .025],...
            'tickdir','both',...
            'xminortick','on');
        set(gca,'CLim',[-.1 .1])
        colorbar('position',[offset+width+.01,offset,.03,height2],'box','off','ticks',[],'color',[1 1 1])
        xlim(xuselims)
        ylim([.5 nModels+.5])
        line([0 0],ylim,'color',[0 0 0],'linewidth',ax_width)

        set(gca,...
            'xticklabel',{0, [], 200, [], 400, [], [] },...
            'xticklabelrotation',0);

        % Mark threshold exceedance
        exc = bf10(:,plot_time,c) >= bf_threshold;
        exc(1,:) = 0;
        exc_diff = diff(exc,1,2);
        idx = 0.5;
        for m = 1:size(exc,1)
            onsets = find(exc_diff(m,:) == 1);
            offsets = find(exc_diff(m,:) == -1) - 1;
            if ~isempty(onsets) && ~isempty(offsets)
                if offsets(1) == 0
                    offsets = offsets(2:end);
                end
                if offsets(1) < onsets(1)
                    onsets = [1 onsets];
                end
                if offsets(end) < onsets(end)
                    offsets = [offsets nSamples];
                end
            elseif ~isempty(onsets)
                offsets = nSamples;
            elseif ~isempty(offsets) && ~(offsets(1) == 0)
                onsets = 1;
            else
                idx = idx + 1;
                continue
            end

            onsets = time(onsets);
            offsets = time(offsets);
            durations = offsets - onsets;
            durations(durations == 0) = 1;

            for r = 1:length(onsets)
                rectangle('Position',[onsets(r) idx durations(r) 1],'edgecolor',[1 1 1],'linewidth',rect_width)
            end

            idx = idx + 1;
        end

        ylims = ylim;
        xlims = xlim;
        line([xlims(1) xlims(1)],[ylims(1) ylims(2)],'color',[0 0 0],'linewidth',ax_width);
        line([xlims(1) xlims(2)],[ylims(2) ylims(2)],'color',[0 0 0],'linewidth',ax_width);

        % Plot colored dots in y axis
        ylocs = get(gca, 'YTick');
        xloc = xlim;
        xText = xloc(1) - 30;  % Adjust as needed

        hold on;
        for i = 1:length(ylocs)
            text(xText, ylocs(i), '●', ...
                'Color', dot_cols(i,:), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 10);
        end
        hold off;

        % 4. Save
        set(f,'PaperPositionMode','auto')
        img_name = [group '_' chan '_BMS_results.tiff'];
        saveas(f, fullfile(bms_dir,img_name));
        clear beta_test
    end

end
