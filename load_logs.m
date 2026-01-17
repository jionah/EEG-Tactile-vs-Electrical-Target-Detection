function [Data, PF] = load_logs(SJs,nRuns,runs,data_dir,ana_dir,group)
%% loads all runs of all subjects into one cell array *Data*

Data = cell(length(SJs), nRuns);
PF = cell(length(SJs),1);
if strcmpi(group,'Tactile'), sfx = 'SomaVibro'; 
elseif strcmpi(group,'Electrical'), sfx = 'periT_FTip'; end

for s = 1:length(SJs)
    cd(fullfile(data_dir, SJs{s}, 'logs'));
    for r = runs{s}   
        log = dir(['log_' sfx '_eeg_' SJs{s} '_' num2str(r) '*.mat']);
        log = log.name;
        Log = load(log);
        evalc(['Data{s,r} = Log.log_' sfx]);
    end
    pf = dir(['psychFunc_' SJs{s} '*.mat']);
    if ~isempty(pf)
        pf = pf.name;
        pf = load(pf);
        PF{s} = pf.PF;
    end
end

cd(ana_dir);
