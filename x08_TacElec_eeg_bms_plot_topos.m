%% Plot BMS results as topoplot

clear
close all

mydir = '...\...\TacElec_EEG'; % Path to project folder

groups = {'Tac', 'Elec'};

xp_threshold    = .99;
beta_threshold = 150;
THRESH         = 1;

topo_template = fullfile(mydir,'templates','topo.set');  % Use template EEGLAB data set to get channel topography
EEG = pop_loadset(topo_template);

bms_dir = fullfile(mydir,'2nd level','BMS','null-int-det-pf-unc-rep-cue');

models      = {'null' 'int', 'det', 'pf', 'unc', 'rep', 'cue'};
isfam       = [0 1 1 1 0 0 0];
plus_fam    = 2;

for g = 1:numel(groups)
    group = groups{g};

    prefix = ['bfraeTMdff' group '_'];

    % BMS file and beta_test
    bms_file = fullfile(bms_dir, [group, '_BMS_FamXPs.mat']);
    load(bms_file)
    beta_test_file = fullfile(bms_dir, [group, '_beta_test.mat']);
    load(beta_test_file)

    % Dummy data sets
    D = spm_eeg_load(fullfile(mydir,'preprocessed',group,[prefix BMS.SJs{1} '.mat']));
    nChannels = numel(BMS.channels);


    %% Prepare plots

    fam_cols = [.5 .5 .5   % null - grey
        .8 .8 .8   % +fam - light grey
        0  1  1   % unc - cyan
        1  0  1   % rep - magenta
        1  1  0   % cue - yellow
        ];

    % Family info
    partition   = [1 2 2 2 3 4 5];
    nFam        = numel(BMS.partition);
    fam_idx     = cell(1,nFam);
    fam_size    = nan(1,nFam);
    for i = 1:nFam
        fam_idx{i} = find(partition == i);
        fam_size(i) = length(fam_idx{i});
    end

    fam_xps = BMS.xp_fam;
    xps = BMS.xp(isfam==1,:,:);
    xps = xps.*repmat(fam_xps(plus_fam,:,:),fam_size(plus_fam),1,1);

    %% Plot
    toi = [0 0.05 0.075 0.1 0.125 0.15 0.25 0.35 0.5];
    plotsamples = D.indsample(toi);

    fig = figure('color','white');
    for s = 1:length(plotsamples)

        subplot(1,9,s);
        ps = plotsamples(s);

        bf10 = nan(nChannels,1);
        tp_chan_cols = cell(nChannels,1);
        for c = 1:nChannels
            [~,max_fam] = max(fam_xps(:,ps,c),[],1);  % Determine winning family

            if max_fam == plus_fam
                [~,max_mod] = max(xps(:,ps,c),[],1);
                bf10(c) = beta_test.bf10(fam_idx{max_fam}(max_mod),ps,c);
                if fam_size(fam_size == 3)
                    tp_chan_cols{c} = xps([2,1,3],ps,c)';
                else
                    tp_chan_cols{c} = [xps(2,ps,c) xps(1,ps,c) 0];
                end
            else
                bf10(c) = beta_test.bf10(fam_idx{max_fam},ps,c);
                tp_chan_cols{c} = fam_cols(max_fam,:)*fam_xps(max_fam,ps,c);
            end
        end

        if THRESH
            data = any(squeeze(fam_xps(:,ps,:))' >= xp_threshold,2) & bf10 >= beta_threshold;
        else
            data = ones(nChannels,1);
        end

        topoplot(data,EEG.chanlocs,'style','blank','efontsize',20,'emarkercolors',tp_chan_cols);
        title(sprintf('%.0f ms',toi(s)*1000),'position',[-.05 .6],'FontSize',14,'FontWeight','b')
        
    end
    
    set(gcf,'Units','centimeters','position',[0,1,30,6])
    set(gcf,'color','white')
    saveas(fig,fullfile(bms_dir, [group '_BMS_Topo.tiff']));

end