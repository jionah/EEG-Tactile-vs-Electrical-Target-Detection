%% Plot ERPs with between-groups comparison BF10 time courses

clear
clc

mydir    = '...\...\TacElec_EEG'; % Path to project folder

bms_dir  = fullfile(mydir,'2nd Level','BMS');
erp_dir  = fullfile(mydir,'2nd level','ERP');
loc_dir  = 'templates';

bf_threshold  = 3; % Bayes factor threshold

leg = {'Miss (Tac)' 'Hit (Tac)' 'Miss (Elec)' 'Hit (Elec)'};
channels = {'CP4' 'C6' 'CPz'};

% Load MEEG dummy object for time & channel labels
D = spm_eeg_load(fullfile(mydir,'preprocessed','Tac','bfraeTMdffTac_VP02.mat'));

% Time stuff
x     = (D.time*1000)';
xlims = [-50 600];

% Get indices of channel labels
cidx = D.indchannel(channels)';

% Model colors
cols = [1 0.4 0.7,  % pink
    1.0, 0.85, 0.4, % yellow
    0.6, 1.0, 0.8,  % mint green
    0.75 0.7 1.0   % lavender
    ];

% Load data to plot
load(fullfile(bms_dir,'BMS_BetweenGroups_all.mat')); % Family level
bf10_allFams = bf10_bg;
load(fullfile(bms_dir,'BMS_BetweenGroups_+Fam.mat')); % only +family
bf10_plusFam = bf10_bg;
load(fullfile(bms_dir,'BMS_BetweenGroups_int-det.mat')); % only int vs. det
bf10_int_det = bf10_bg;
load(fullfile(bms_dir,'BMS_BetweenGroups_det-pf.mat')); % only det vs. pf
bf10_det_pf = bf10_bg;

%--------------------------------------------------------------------------
%% Plot ERPs
% Get ERPs
load(fullfile(erp_dir,'Tac_ERPs_det.mat'),'ERP');
mean_ERP_tac = squeeze(mean(ERP,4));
sd_ERP_tac = squeeze(std(ERP,[],4));
se_ERP_tac = sd_ERP_tac./sqrt(size(ERP,4));
clear ERP;

load(fullfile(erp_dir,'Elec_ERPs_det.mat'),'ERP');
mean_ERP_elec = squeeze(mean(ERP,4));
sd_ERP_elec = squeeze(std(ERP,[],4));
se_ERP_elec = sd_ERP_elec./sqrt(size(ERP,4));
clear ERP;

nConds = size(mean_ERP_tac,3); % miss/hit
fig_pos = [1,1,15,10];

erp_cols = [ 0.49, 0.73, 0.91;  % Maya blue
             0.89, 0.15, 0.21]; % Alizarin crimson

erp_lines = {':','-'};

offset = .1;
width = .8;
height = .88;

for c = cidx

    f = figure('units','centimeters','position',fig_pos);
    subplot('position',[offset,offset,width,height])
    hold on
    set(gca,'FontName','Calibri','FontSize',10)

    min_y = 0;
    max_y = 0;

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
    
    chanlabel = D.chanlabels{c};
    name = [chanlabel '_det.tiff'];

    line([0 0],ylims,'color',[0 0 0],'linewidth',1)
    line(xlims,[0 0],'color',[0 0 0],'linewidth',1)
    line(xlims,[ylims(1) ylims(1)],'color',[0 0 0],'linewidth',1)
    line([xlims(1) xlims(1)],[ylims(1) ylims(2)],'color',[0 0 0],'linewidth',1)
    ylim(ylims);
    % legend(h, leg);

    % Plot channel label
    text(520, ylims(1)+diff(ylims)*0.08, chanlabel,'fontsize',20')

    %% Plot evidence time courses of between-groups differences
    yl = ylim;
    ymin = yl(1);
    ymax = yl(2);
    range_y = ymax - ymin;
    stack_rel = [0.01, 0.04, 0.07];  % controls the vertical height of the time courses relative to ymin

    exc_fams    = bf10_allFams(c,:) <= bf_threshold; % family level: no evidence for a difference
    exc_plusFam = bf10_plusFam(c,:) >= bf_threshold; % within +fam: evidence for a difference
    exc_int_det = bf10_int_det(c,:) >= bf_threshold; % only int-det: evidence for a difference
    exc_det_pf  = bf10_det_pf(c,:)  >= bf_threshold; % only det-pf: evidence for a difference

    for s = 1:length(x)-1
        if exc_fams(s)
            yline_pos = ymin + stack_rel(1)*range_y;
            plot(x([s s+1]),[yline_pos yline_pos],'Color',cols(1,:),'LineWidth',7)
        end
        if strcmpi(chanlabel,'CP4')
            if exc_int_det(s)
                yline_pos = ymin + stack_rel(2)*range_y;
                plot(x([s s+1]),[ yline_pos  yline_pos],'Color',cols(3,:),'LineWidth',7)
            end
        else
            if exc_det_pf(s)
                yline_pos = ymin + stack_rel(2)*range_y;
                plot(x([s s+1]),[ yline_pos  yline_pos],'Color',cols(4,:),'LineWidth',7)
            end
        end
        if exc_plusFam(s)
            yline_pos = ymin + stack_rel(3)*range_y;
            plot(x([s s+1]),[ yline_pos  yline_pos],'Color',cols(2,:),'LineWidth',7)

        end
    end
    
    % Plot vertical lines at time points of interest (closest possible to
    % conventional nominal [50 ms, 100 ms, 300 ms] if all criteria are met;
    % conventional nominal [80 ms, 140 ms] if criteria are never met around
    % component)
    if strcmpi(chanlabel,'cp4'), conv_times = 50;
        time_points = find(bf10_plusFam(c,:) > bf_threshold & bf10_int_det(c,:) > bf_threshold & bf10_allFams(c,:) < bf_threshold); % all criteria met
    elseif strcmpi(chanlabel,'c6'), conv_times = [80; 140];
        time_points = dsearchn(x,conv_times)'; % criteria not met in relevant time range; use conventional nominal time points
    elseif strcmpi(chanlabel,'cpz'), conv_times = [100 ;300];
        time_points = find(bf10_plusFam(c,:) > bf_threshold & bf10_det_pf(c,:) > bf_threshold & bf10_allFams(c,:) < bf_threshold); % all criteria met
    end
    time_points = time_points(dsearchn(time_points', dsearchn(x,conv_times))); % closest to conventional nominal time
    
    % Plot lines
    xline(x(time_points),'color',[.7 .7 .7],'layer','bottom','linew',1.75)
    
    tightfig;
    fname = [D.chanlabels{c} '_det_BMS_BetweenGroups.tiff'];
    saveas(f,fullfile(fullfile(bms_dir,fname)))
   
end % of cidx loop
