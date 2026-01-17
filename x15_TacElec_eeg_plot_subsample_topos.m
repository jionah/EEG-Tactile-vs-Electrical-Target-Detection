%% Plot topographies for control analysis (voltage comparisons)

clear
close all
clc

mydir    = '...\...\TacElec_EEG'; % Path to project folder

ca_dir  = fullfile(mydir,'2nd level','control analysis');
loc_dir = fullfile(mydir,'templates');

EEG     = pop_loadset(fullfile(loc_dir,'topo.set'));

components  = { 'P50', 'N80', 'P100', 'N140', 'P300', 'P50'};
times       = [.05 .08 .1 .14 .3 .05]; % conventional nominal ERP component peak times (in seconds)
threshold   = 3;  % BF10 threshold for the black asterisks
high_thresh = 10; % BF10 threshold for the green asterisks

% Load subsample average ERPs
load(fullfile(ca_dir, 'Tac_ERPs_det_intmatched_avg40.mat'));
ERP_tac = ERP; clear ERP
load(fullfile(ca_dir, 'Elec_ERPs_det_intmatched_avg40.mat'));
ERP_elec = ERP; clear ERP

% Define time vector
D = spm_eeg_load(fullfile(mydir,'preprocessed','Tac','bfraeTMdffTac_VP02.mat')); % Load dummy MEEG
x = D.time;
t_idcs = dsearchn(x', times'); % find indices ERP component peak times

% % Convert electrode labels to indices
chlabels = D.chanlabels;

% Topo & marker parameters
maplims  = [-1 1]; % Color map boundaries
splims   = [-.6 .6]; % subplot limits (to get whole nose & ears)
emarker = '*'; % marker: asterisk
dotsize = 7; % marker size
emk     = 'k'; % black for BF10 > 3
emk2    = 'g'; % green for BF10 > 10

%--------------------------------------------------------------------------
% Loop over components
for c = 1:numel(components)
    component = components{c};

    if c < 6 % Plot topographies at nominal component peak times
        t = t_idcs(c);

        % Data (hit minus miss)
        tac_data  = squeeze(ERP_tac(:,t,2,:) - ERP_tac(:,t,1,:));
        elec_data = squeeze(ERP_elec(:,t,2,:) - ERP_elec(:,t,1,:));

    else % Plot P50 topography again with actual peak times
        
        % Find actual peaks as local maxima nearest to 50 ms
        ch_idx  = find(ismember(chlabels,'CP4')); % Get channel index

        % Time courses of condition-wise grand average ERPs for CP4
        y_th  = mean(ERP_tac(ch_idx,:,2,:),4); % tactile, hit
        y_tm  = mean(ERP_tac(ch_idx,:,1,:),4); % tactile, miss
        y_eh = mean(ERP_elec(ch_idx,:,2,:),4); % electrical, hit
        y_em = mean(ERP_elec(ch_idx,:,1,:),4); % electrical, miss

        % Local maxima vectors
        max_th = find(islocalmax(y_th));
        max_tm = find(islocalmax(y_tm));
        max_eh  = find(islocalmax(y_eh));
        max_em = find(islocalmax(y_em));

        % Closest to nominal P50 peak time
        nom_time = .05;
        nom_samp = dsearchn(x', nom_time);

        peaks = [ max_th(dsearchn(max_th', nom_samp)),
            max_tm(dsearchn(max_tm', nom_samp)),
            max_eh(dsearchn(max_eh', nom_samp)),
            max_em(dsearchn(max_em', nom_samp))
            ];

        % Print actual component peak times to command window
        fprintf('Actual P50 peak times: %.1f (tac hit), %.1f (tac miss), %.1f (elec hit), %.1f (elec miss), ', x(peaks)*1000);

        % Data (hit minus miss)
        tac_data  = squeeze(ERP_tac(:,peaks(1),2,:) - ERP_tac(:,peaks(2),1,:)); % hit peak minus miss peak
        elec_data = squeeze(ERP_elec(:,peaks(3),2,:) - ERP_elec(:,peaks(4),1,:)); % hit peak minus miss peak

    end

    % Compute BF10 between group-specific difference topographies
    for ch = 1:size(tac_data,1)
        bf10_tac(ch)  = bf.ttest( tac_data(ch,:) );
        bf10_elec(ch) = bf.ttest( elec_data(ch,:) );
        bf10_diff(ch) = bf.ttest2( tac_data(ch,:), elec_data(ch,:) );
    end

    %% Plot voltage topographies with BF10 markers:
    % tac hit-miss (paired), elec hit-miss (paired), diffs-diff (two-sample)
    fig = figure;
    fig.Position = [1 1 500 250];

    %--------------------------------
    % Plot tac voltage difference map
    subplot(131)
    plotchans = find(ismember(chlabels,chlabels(bf10_tac > threshold))); % find marker positions for lower BF10 threshold
    plotchans2 = find(ismember(chlabels,chlabels(bf10_tac > high_thresh))); % find marker positions for higher BF10 threshold

    topoplot(mean(tac_data,2), EEG.chanlocs, 'maplimits',maplims,'style','map','electrodes','on',...
        'emarker2',{plotchans emarker emk dotsize .8});

    % Overlay BF10 > 10 markers in green
    hold on
    topoplot(mean(tac_data,2), EEG.chanlocs, 'maplimits',maplims,'style','map','electrodes','on',...
        'emarker2',{plotchans2 emarker emk2 dotsize .8});

    %---------------------------------
    % Plot elec voltage difference map
    subplot(132)
    plotchans = find(ismember(chlabels,chlabels(bf10_elec > threshold)));
    plotchans2 = find(ismember(chlabels,chlabels(bf10_elec > high_thresh)));

    topoplot(mean(elec_data,2), EEG.chanlocs, 'maplimits',maplims,'style','map','electrodes','on',...
        'emarker2',{plotchans emarker emk dotsize .8});

    % Overlay BF10 > 10 markers in green
    hold on
    topoplot(mean(elec_data,2), EEG.chanlocs, 'maplimits',maplims,'style','map','electrodes','on',...
        'emarker2',{plotchans2 emarker emk2 dotsize .8});

    %-----------------------------------
    % Plot difference of differences map
    subplot(133)
    plotchans  = find(ismember(chlabels,chlabels(bf10_diff > threshold)));
    plotchans2 = find(ismember(chlabels,chlabels(bf10_diff > high_thresh)));

    topoplot(mean(tac_data,2)-mean(elec_data,2), EEG.chanlocs, 'maplimits',maplims,'style','map','electrodes','on',...
        'emarker2',{plotchans emarker emk dotsize .8});

    % Overlay BF10 > 10 markers in green
    hold on
    topoplot(mean(tac_data,2)-mean(elec_data,2), EEG.chanlocs, 'maplimits',maplims,'style','map','electrodes','on',...
        'emarker2',{plotchans2 emarker emk2 dotsize .8});

    set(gcf,'color','w');

    % Save
    tightfig;
    filename = sprintf('%s_voltage_topos_BF10',component);
    if c == 6, filename = [filename '_peaktimes']; end
    saveas(fig,fullfile(ca_dir, [filename '.tiff']));

end % of components loop





