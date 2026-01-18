# EEG-Tactile-vs-Electrical-Target-Detection
Matlab scripts to analyze data from an EEG study on somatosensory target detection:

Förster, J., Vardiero, G., Nierhaus, T., & Blankenburg, F. (2025). ERP responses reveal different neural mechanisms for perception of electrical and tactile stimuli. Consciousness and Cognition 135:103935. https://doi.org/10.1016/j.concog.2025.103935

Preprocessed data to use with this code are available on Figshare: https://dx.doi.org/10.6084/m9.figshare.31081684

## Note:
Scripts x00-x08 rely heavily on work by Pia Schröder, to which I owe a great deal: https://github.com/PiaSchroeder/SomatosensoryTargetDetection_EEG/tree/main

## Requires:
- SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
- EEGLAB: https://eeglab.org/download/
- bayesFactor Toolbox: https://github.com/klabhub/bayesFactor/
- Variational Bayesian Analysis Toolbox: https://mbb-team.github.io/VBA-toolbox/

## Scripts

### Behavioral analysis:

  - x00_TacElec_eeg_check_pfs.m
  - x01_TacElec_eeg_behaviour.m

### Bayesian 1st level GLMs:

  - x02_TacElec_eeg_bayesglm_erp.m
  - bayesglm_sensors.m

### Bayesian Model Selection within groups:

  - x03_TacElec_eeg_bms_erp.m
  - x04_TacElec_eeg_compute_int_det_erps.m
  - x05_TacElec_eeg_bms_results.m
  - x06_TacElec_eeg_compute_intmatched_det_erps.m
  - x07_TacElec_eeg_channel_count.m
  - x08_TacElec_eeg_bms_plot_topos.m

### Bayesian Model Selection between groups:

  - x09_TacElec_eeg_bms_between_groups.m
  - x10_TacElec_eeg_bms_between_groups_results.m
  - x11_TacElec_eeg_bms_between_groups_topos.m

### Control analysis:

   - x13_TacElec_eeg_compute_subsample_average_erps.m
   - x14_TacElec_eeg_plot_subsample_erps.m
   - x15_TacElec_eeg_plot_subsample_topos.m
