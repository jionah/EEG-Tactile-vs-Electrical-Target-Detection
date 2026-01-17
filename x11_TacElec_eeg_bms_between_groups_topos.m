%% Plot BF10 topographies for two-model comparisons at selected time points
clear
close all
clc

mydir    = '...\...\TacElec_EEG'; % Path to project folder

bms_dir = fullfile(mydir,'2nd Level','BMS');
erp_dir = fullfile('2nd level','ERP');
loc_dir = fullfile(mydir,'templates');

bf_threshold  = 3;   % Bayes factor threshold (positive evidence)
bf_highthresh = 10;  % High sthreshold (stronger evidence)

% Components and their time points to plot topographies for
components  = {'P50' 'N80' 'N140' 'P100' 'P300'};
time_points = [54 68 99 78 190]; % see x10_TacElec_eeg_bms_between_groups_results.m

% Model colors
cols = [0.6, 1.0, 0.8,  % mint green (int-det)
        0.75 0.7 1.0   % lavender (det-pf)
        0 0 0]; % black (high BF threshold)

% Load data to plot
load(fullfile(bms_dir,'BMS_BetweenGroups_int-det.mat')); % int vs. det
bf10_int_det = bf10_bg;
load(fullfile(bms_dir,'BMS_BetweenGroups_det-pf.mat')); % det vs. pf
bf10_det_pf = bf10_bg;

%% Plot BF10 topographies for two-model comparisons
EEG        = pop_loadset(fullfile(loc_dir,'topo.set'));
nChannels  = EEG.nbchan;

fig_pos    = [1 1 40 8];
f = figure('units','centimeters','position',fig_pos);
for c = 1:numel(components)
    comp = components{c};
    subplot(1,5,c)
    if strcmp(comp,'P50')
        markchans  = repmat({[.5 .5 .5]}, nChannels, 1);
        markchans(bf10_int_det(:,time_points(c)) > bf_threshold)  = {cols(1,:)}; % [BF10 > threshold] marker positions
        markchans(bf10_int_det(:,time_points(c)) > bf_highthresh) = {cols(3,:)}; % [BF10 > higher threshold] marker positions
        topoplot(bf10_int_det(:,time_points(c)) > bf_threshold,EEG.chanlocs,'style','blank','efontsize',20,'emarkercolors',markchans);
    else
        markchans  = repmat({[.5 .5 .5]}, nChannels, 1);
        markchans(bf10_det_pf(:,time_points(c)) > bf_threshold)  = {cols(2,:)}; 
        markchans(bf10_det_pf(:,time_points(c)) > bf_highthresh) = {cols(3,:)};
        topoplot(bf10_det_pf(:,time_points(c)) > bf_threshold,EEG.chanlocs,'style','blank','efontsize',20,'emarkercolors',markchans);
    end

    % Plot ERP component name as subplot title
    title(comp,'fontsize',20)
    
end
set(gcf,'color','w');
tightfig;

 % Save
saveas(f,fullfile(bms_dir,'BMS_BetweenGroups_topos.tiff'));
