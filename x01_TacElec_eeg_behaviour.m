%% TacElec EEG behaviour
% Detection rates
% Reaction times
% Response associations: Detection - Match reports

clear 
clc

mydir = '...\...\TacElec_EEG'; % Path to project folder

group = 'Tac'; % 'Tac' 'Elec' (run separately for each group)

%% Directories etc
if strcmpi(group, 'Tac')
    SJs  = { 'VP02' 'VP03' 'VP04' 'VP05' 'VP06' 'VP07' 'VP08' 'VP09' 'VP10' 'VP11' 'VP12' 'VP15' 'VP16' 'VP17' 'VP18' 'VP19' 'VP20' 'VP21' 'VP22' 'VP24' 'VP25' 'VP27' 'VP28' 'VP29' 'VP30'};
    gend = {  'w'    'w'    'm'    'm'    'w'    'w'    'm'    'w'    'w'    'w'    'w'    'w'    'm'    'm'    'w'    'w'    'm'    'm'    'm'    'w'    'w'    'm'    'm'    'w'    'm'  };
    age  = [   19     24     28    21     25      27    24      21     28     25     22    29      28    26      25     32    28     23      27     22     26    24     23     23     32  ];
elseif strcmpi(group, 'Elec')
    SJs   = { 'S01' 'S02' 'S03' 'S06' 'S07' 'S08' 'S10' 'S12' 'S16' 'S19' 'S20' 'S21' 'S22' 'S25' 'S26' 'S27' 'S28' 'S29' 'S30' 'S31' 'S33' 'S35' 'S36' 'S37' 'S38'};
    gend  = {  'm'   'w'   'm'   'w'   'w'   'w'   'w'   'w'   'w'   'm'   'w'   'm'   'w'   'm'   'w'   'm'   'w'   'm'   'w'   'm'   'w'   'w'   'm'   'w'   'w'};
    age   = [   28    23    35    25    29    25    23    30    21    33    24    28    25    27    24    33    33    25    31    24    26    26    28    23    21];
end

log_dir = 'logs';
log_f   = '_trial_log';
trg_dir = fullfile(mydir,'2nd level','Behaviour');
if ~exist(trg_dir,'dir')
    mkdir(trg_dir)
end

%% Detection rates
det_rates = nan(numel(SJs),1);

for s = 1:numel(SJs)
    
    sj_trl_log = fullfile(mydir, log_dir, group, [SJs{s} log_f '.mat']);
    load(sj_trl_log)
    
    [~,det_idx] = ismember('det',trial_log.labels);
    
    det_rates(s) = sum(trial_log.data(det_idx,:))/size(trial_log.data,2)*100;
  
end

mean_det_rates = mean(det_rates);
sd_det_rates = std(det_rates);

DetRates.all = det_rates;
DetRates.mean = mean_det_rates;
DetRates.sd = sd_det_rates;

save(fullfile(trg_dir, [group '_DetRates.mat']),'DetRates')

%% Reaction times
RTs.all = [];
RTs.sj_mean = nan(1,numel(SJs));

for s = 1:numel(SJs)
    
    sj_trl_log = fullfile(mydir, log_dir, group, [SJs{s} log_f '.mat']);
    load(sj_trl_log)
    
    [~,rt_idx] = ismember('rts',trial_log.labels);
    rt_data = trial_log.data(rt_idx,:)*1000;
    
    RTs.all = [RTs.all rt_data];
    RTs.sj_mean(s) = mean(rt_data);
    
    [~,cidx] = ismember('det',trial_log.labels);
    c_data = trial_log.data(cidx,:);
    cs = unique(c_data);

    [~,ridx] = ismember('rep',trial_log.labels);
    r_data = trial_log.data(ridx,:);
    rs = unique(r_data);
    
    for cc = 1:length(cs)
        RTs.det.sj_mean(s,cc) = mean(rt_data(c_data==cs(cc)));
        RTs.rep.sj_mean(s,cc) = mean(rt_data(r_data==rs(cc)));
        for rr = 1:length(rs)
            RTs.detrep.sj_mean(s,cc,rr) = mean(rt_data(c_data==cs(cc) & r_data==rs(rr)));
        end
    end
end

RTs.group_mean = mean(RTs.sj_mean);
RTs.group_sd = std(RTs.sj_mean);

RTs.det.labels = {'miss' 'hit'};
RTs.det.group_mean = mean(RTs.det.sj_mean);
RTs.det.group_sd = std(RTs.det.sj_mean);
RTs.det.BF10 = bf_ttest(RTs.det.sj_mean(:,1),RTs.det.sj_mean(:,2)); % Uses the ttest.m function from the bayesFactor toolbox (Bart Krekelberg (2021).
                                                                    % bayesFactor (https://github.com/klabhub/bayesFactor)), renamed to bf_ttest.m
RTs.rep.labels = {'mismatch' 'match'};                              % to avoid naming conflicts  
RTs.rep.group_mean = mean(RTs.rep.sj_mean);
RTs.rep.group_sd = std(RTs.rep.sj_mean);
RTs.rep.BF10 = bf_ttest(RTs.rep.sj_mean(:,1),RTs.rep.sj_mean(:,2));

RTs.detrep.labels = {'miss-mismatch' 'miss-match'; 'hit-mismatch' 'hit-match'};
RTs.detrep.group_mean = squeeze(mean(RTs.detrep.sj_mean));
RTs.detrep.group_sd = squeeze(std(RTs.detrep.sj_mean));
RTs.detrep.BF10 = bf_ttest(RTs.detrep.sj_mean(:,1),RTs.detrep.sj_mean(:,2));

save(fullfile(trg_dir,[group '_RTs.mat']),'RTs')

%% Response associations: Detection (det) - Match reports (rep)
BF10_det_rep = nan(numel(SJs),1);

for s = 1:numel(SJs)
    
    sj_trl_log = fullfile(mydir, log_dir, group, [SJs{s} log_f '.mat']);
    load(sj_trl_log)
    
    [~,det_idx] = ismember('det',trial_log.labels);
    [~,rep_idx] = ismember('rep',trial_log.labels);
    
    dets = trial_log.data(det_idx,:);
    reps = trial_log.data(rep_idx,:);

    det_rep(1,1) = sum(dets==0 & reps==0);
    det_rep(1,2) = sum(dets==0 & reps==1);
    det_rep(2,1) = sum(dets==1 & reps==0);
    det_rep(2,2) = sum(dets==1 & reps==1);
    BF10_det_rep(s) = c_table(det_rep);     % Uses the c_table.m function from the companion software to Johnson, V. E., & Albert, J. H. 
                                            % (2006). Ordinal data modeling. Springer Science & Business Media. Available at: 
                                            % https://de.mathworks.com/matlabcentral/fileexchange/2264-ordinal-data-modeling 
    
end

BF01_det_rep = 1./BF10_det_rep;

Det_Rep.BF10 = BF10_det_rep;
Det_Rep.BF01 = BF01_det_rep;

save(fullfile(trg_dir, [group '_ResponseAssociations.mat']),'Det_Rep')

%% Display stats
disp('---------------------------------------------------------------------------------')

disp('Descriptive stats')
disp('---------------------------------------------------------------------------------')
fprintf('N=%d (%d female)\n',numel(SJs), sum(ismember(gend,'w')))
fprintf('Mean age: %.2f +- %.2f years\n', mean(age), std(age))
fprintf('Mean detection rate (%%): %.2f +- %.2f\n',DetRates.mean,DetRates.sd)
disp('Reaction times:')
fprintf('    Hits: %.2f +- %.2f\n    Misses: %.2f +- %.2f\n    BF10: %.2f\n',RTs.det.group_mean(2), ...
    RTs.det.group_sd(2), RTs.det.group_mean(1), RTs.det.group_sd(1), RTs.det.BF10)
fprintf('Dissociation (target detection from overt report): ')
thresh_range = [3 10];
if all(thresh_range(1) < Det_Rep.BF01) & all(Det_Rep.BF01 < thresh_range(2))
    fprintf('%d < BF01 < %d for all participants\n', thresh_range)
end


