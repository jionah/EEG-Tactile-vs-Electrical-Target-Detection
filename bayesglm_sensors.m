function results = bayesglm_sensors(S)
% Fit GLM to sensor data using Bayesian estimation scheme
%
% Input:
%   S.D         - MEEG object containing trialwise time courses
%   S.chantype  - Channel type that should be included (EEG|LFP)
%   S.reg       - nxm model regressor (n trials, m regressors)
%   S.alpha     - prior precision of regression coefficients
%
% Adapted from spm_spm_vb.m to analyse sensor-level EEG data
%

try
    D = S.D;
    chantype = S.chantype;
    reg = S.reg;
    alpha = S.alpha;
catch
    error('Missing inputs!')
end

%% get data

channels = D.selectchannels(chantype);
Y = D(channels,:,:);
Y = zscore(Y,[],3);

nSamples = D.nsamples;
nSensors = length(channels);
nTrials = size(Y,3);

%% Set up design matrix

X = [ones(nTrials,1), reg];

nBeta = size(X,2);
nPsd = nBeta;
AR_P = 0;

%% Define block template
block_template               = spm_vb_init_volume(X,AR_P);
block_template.maxits        = 4;
block_template.tol           = 0.0001;
block_template.compute_det_D = 0;
block_template.verbose       = 0;
block_template.update_w      = 1;
block_template.update_lambda = 1;
block_template.update_F      = 1;

%% Specify priors
priors.W = 'Voxel - Uninformative'; % prior for regression coefficients
priors.A = 'Voxel - Uninformative'; % prior for AR coefficients

%% Run VB

% Prepare data structures
dims  = 'betas x channels x samples'; 
beta  = zeros(nBeta,nSensors,nSamples);
Psd   = zeros(nPsd,nSensors,nSamples);
LogEv = zeros(nSamples,nSensors);
Hp = zeros(nSamples,nSensors); % Noise variance
AR = zeros(AR_P,nSensors,nSamples);

fprintf('Sample ');
nbytes = fprintf('0/%d',nSamples);
for t = 1:nSamples
    
    % Display block number
    while nbytes > 0
        fprintf('\b');
        nbytes = nbytes - 1;
    end
    nbytes = fprintf('%d/%d',t,nSamples);
    
    % Get data
    Y_t = permute(Y(:,t,:),[3,1,2]);
    
    % Estimate model
    
    % Set priors
    block = spm_vb_set_priors(block_template,priors,(1:nSensors)');
    block.mean_alpha = alpha*ones(block.k,1);
    block.b_alpha_prior  = block.mean_alpha ./ block.c_alpha_prior;

    % Fit model
    block = spm_vb_glmar(Y_t,block);
    
    % Report AR values
    if AR_P > 0
        AR(1:AR_P,:,t) = block.ap_mean;
    end
    
    % Update LogEv
    progress(t).F = block.F;
    Fn = spm_vb_Fn(Y_t,block);
    LogEv(t,:) = LogEv(t,:)+Fn;
    
    % Update regression coefficients
    beta(:,:,t) = block.wk_mean;
    
    % Report noise variance
    Hp(t,:) = sqrt(1./block.mean_lambda');
    
    % Store regression coefficient posterior standard deviations
    Psd(:,:,t) = block.w_dev;
    
    % Prior precision
    progress(t).mean_alpha = block.mean_alpha;
    
    % Get block-wise Taylor approximation to posterior correlation
    block = spm_vb_taylor_R(Y_t,block);
    progress(t).mean = block.mean;
    progress(t).elapsed_seconds = block.elapsed_seconds;
    
    % Save coefficient RESELS and number of voxels
    progress(t).gamma_tot = block.gamma_tot;
    progress(t).N = block.N;
    
    clear block Y_t
    
end

%% Assemble results
results.betadims        = dims;         % beta dimensions
results.beta            = beta;         % Estimated regression coefficients
results.Psd             = Psd;          % Posterior standard deviations of regression coefficients
results.LogEv           = LogEv;        % Log evidence
results.Hp              = Hp;           % Noise variance
results.AR_P            = AR_P;         % AR model order
results.AR              = AR;           % AR values
results.progress        = progress;     % Estimation progress

fprintf(' Done!\n')

