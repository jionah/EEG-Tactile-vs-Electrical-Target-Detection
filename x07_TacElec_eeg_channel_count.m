%% Plot temporal evolution of model fits across the scalp
% For each model: significant electrodes/all electrodes

clear
close all

mydir = '...\...\TacElec_EEG'; % Path to project folder

groups = {'Tac', 'Elec'};

xp_threshold   = .99;
beta_threshold = 150;

%% Directories etc.

bms_dir = fullfile(mydir,'2nd level','BMS','null-int-det-pf-unc-rep-cue');

uselims = [-50 600];
xticks = 0:100:600;
yticks = 0:10:100;
plot_mods = [1 5 6 7 2 4 3]; % plot models in this order

for g = 1:numel(groups)

    group = groups{g};

    load(fullfile(bms_dir,[group, '_BMS_FamXPs.mat']))
    load(fullfile(bms_dir,[group, '_beta_test.mat']))

    nChannels = length(BMS.channels);
    nSamples = length(BMS.time);

    %% Family info
    partition   = [1 2 2 2 3 4 5];
    nMods       = numel(BMS.models);
    nFam        = numel(BMS.partition);
    fam_idx     = cell(1,nFam);
    fam_size    = nan(1,nFam);
    for i = 1:nFam
        fam_idx{i} = find(partition == i);
        fam_size(i) = length(fam_idx{i});
    end

    %% Count channels
    XPs = BMS.xp(fam_idx{fam_size > 1},:,:);
    FamXPs = BMS.xp_fam;
    model_count = nan(nMods,nSamples);
    [~,max_fams] = max(FamXPs,[],1);

    idx = 1;
    for f = 1:nFam

        if fam_size(f) > 1
            [~,max_mods] = max(XPs,[],1);
            beta = nan(size(max_mods));
            for s = 1:nSamples
                for c = 1:nChannels
                    beta(1,s,c) = beta_test.bf10(fam_idx{f}(max_mods(1,s,c)),s,c);
                end
            end
            %     winning family & XP threshold exceedance      & beta threshold exceedance
            filt = max_fams == f & FamXPs(f,:,:) >= xp_threshold & beta >= beta_threshold;
            max_mods(~filt) = 0;
            for ff = 1:fam_size(f)
                model_count(idx,:) = sum(max_mods == ff,3)./nChannels;
                idx = idx+1;
            end
        else
            beta = beta_test.bf10(fam_idx{f},:,:);
            filt = max_fams == f & FamXPs(f,:,:) >= xp_threshold & beta >= beta_threshold;
            model_count(idx,:) = sum(filt,3)./nChannels;
            idx = idx+1;
        end
    end

    %% Plot

    cols = [.7 .7 .7
        0  1  0
        1  0  0
        0  0  1
        0  1  1
        1  0  1
        1  1  0];

    fig_pos = [1 2 16 4]*2;
    subplot_pos = [.05 .1 .9 .8];
    fig = figure('units','centimeters','position',fig_pos);
    subplot('Position',subplot_pos)
    hold on
    set(gca,'FontName','Calibri','FontSize',18,'ytick',yticks,'xtick',xticks,'xminortick','on','yminortick','on')

    for f = plot_mods
        plot(BMS.time*1000,model_count(f,:)*100,'color',cols(f,:),'linewidth',2.5)
    end

    % xlabel('Time (ms)')
    ylabel('% of electrodes')
    xlim(uselims)
    ylim([0 45])

    % Plot vertical lines at time points of interest
    toi = [0 0.05 0.075 0.1 0.125 0.15 0.25 0.35 0.5];
    xline(toi*1e3,'color',[.8 .8 .8],'layer','bottom','linew',1.5)

    % Plot time courses
    line([0 0],ylim,'color',[0 0 0],'linewidth',1)
    line([-50 -50],ylim,'color',[0 0 0],'linewidth',1)
    line(xlim,[0 0],'color',[0 0 0],'linewidth',1)

    saveas(fig,fullfile(bms_dir, [group '_BMS_channel_count.tiff']))

end