%% Plot ERPs averaged over subsamples

clear
close all
clc

mydir    = '...\...\TacElec_EEG'; % Path to project folder

ca_dir = fullfile(mydir, '2nd level', 'control analysis');
erp_dir = fullfile(mydir,'2nd level','ERP');
trg_dir = fullfile(mydir,'2nd level','control analysis');

N = 25; % 25 participants in each group
nConds = 2; % hit/miss

leg = {'Miss (Tac)' 'Hit (Tac)' 'Miss (Elec)    test' 'Hit (Elec)'};

fig_pos = [1,1,15,10];
channels = {'CP4' 'C6' 'CPz'};

xlims = [-50 600];


% Load EEG data dummy
D = spm_eeg_load(fullfile(mydir,'preprocessed','Tac','bfraeTMdffTac_VP02.mat'));

%% Get ERPs
load(fullfile(ca_dir,'Tac_ERPs_det_intmatched_avg40.mat'))
mean_ERP_tac = squeeze(mean(ERP,4));
sd_ERP_tac = squeeze(std(ERP,[],4));
se_ERP_tac = sd_ERP_tac./sqrt(N);
clear ERP;

load(fullfile(ca_dir,'Elec_ERPs_det_intmatched_avg40.mat'))
mean_ERP_elec = squeeze(mean(ERP,4));
sd_ERP_elec = squeeze(std(ERP,[],4));
se_ERP_elec = sd_ERP_elec./sqrt(N);
clear ERP;

%% Plot

erp_cols = [ 0.49, 0.73, 0.91;  % Maya blue
             0.89, 0.15, 0.21]; % Alizarin crimson
erp_lines = {':','-'};

offset = .1;
width  = .8;
height = .88;

cidx = D.indchannel(channels);
for c = cidx'

    f = figure('units','centimeters','position',fig_pos);
    subplot('position',[offset,offset,width,height])
    hold on
    set(gca,'FontName','Calibri','FontSize',10)

    min_y = 0;
    max_y = 0;

    x = (D.time*1000)';

    h = gobjects(4,1); % Preallocate for 4 handles

    % Tactile
    for d = 1:nConds
        y = mean_ERP_tac(c,:,d)';
        se = se_ERP_tac(c,:,d)';
        fill([x;flipud(x)],[y-se;flipud(y+se)],erp_cols(1,:),'linestyle','none','facealpha',.2);
        h(d) = line(x,y,'Color',erp_cols(1,:),'Linestyle',erp_lines{d},'LineWidth',2);
        max_y = max([max_y max(y(x<=xlims(2))+se(x<=xlims(2)))]);
        min_y = min([min_y min(y(x<=xlims(2))-se(x<=xlims(2)))]);
    end

    % Electrical
    for d = 1:nConds
        y = mean_ERP_elec(c,:,d)';
        se = se_ERP_elec(c,:,d)';
        fill([x;flipud(x)],[y-se;flipud(y+se)],erp_cols(2,:),'linestyle','none','facealpha',.2);
        h(d+nConds) = line(x,y,'Color',erp_cols(2,:),'Linestyle',erp_lines{d},'LineWidth',2);
        max_y = max([max_y max(y(x<=xlims(2))+se(x<=xlims(2)))]);
        min_y = min([min_y min(y(x<=xlims(2))-se(x<=xlims(2)))]);
    end

    ylim([min_y max_y])
    ylims = ylim;
    xlim(xlims)

    set(gca,...
        'xtick',0:200:600,...
        'fontsize',20,...
        'ticklength',[.015 .025],...
        'tickdir','both',...
        'xminortick','on')

    name = [D.chanlabels{c} '_det_intmatched_avg40.tiff'];

    line([0 0],ylims,'color',[0 0 0],'linewidth',1)
    line(xlims,[0 0],'color',[0 0 0],'linewidth',1)
    line(xlims,[ylims(1) ylims(1)],'color',[0 0 0],'linewidth',1)
    line([xlims(1) xlims(1)],[ylims(1) ylims(2)],'color',[0 0 0],'linewidth',1)
    ylim(ylims);
    % legend(h, leg);

    % Plot channel label
    text(3, ylims(1)+diff(ylims)*0.96, D.chanlabels{c},'fontsize',20) 

    % Save
    tightfig;
    saveas(f,fullfile(trg_dir,[D.chanlabels{c} '_det_intmatched_avg40.tiff']))
end
