function [boldpat_vet, N_trial, subj] = arbMBMF_boldpat(Exp, ROI, ID, varargin)
% return EPI volume as masked time series

% boldpat_vet: roi_size x length(event) x N_trial
% multi-dimensional array with Voxel, Event, Trial dimension

% =========================================================================
% Name-value pairs
% 'TrialCondition' - 
% 'Event' - within-trial event index (1:f1, 2:S1, 3:A1, 4:f2, 5:S2, 6:A2, 7:f3, 8:R)
%           possible range: [-7:0, 1:8, 9:16] (prev, curr, next trial)
% 'roi_type' - 'sphere10', 'sphere5', 'AAL3'

% =========================================================================

% root = '/home/ydsung'; is_cluster = 1;
root = 'data'; is_cluster = 0;

% default input parameters
options = struct('SignalType', 'zscore', ...
                'MaskNum', [], ...
                'Event', 1:8, ...
                'TR', 2.78, ...
                'Delay', 6.1, ...
                'TrialCondition', 'none', ...
                'strict', [], ...
                'interpolation', [], ...
                'roi_type', 'AAL3', ...
                'fextension', '.nii', ...
                'prefix', 'wra');
% read input parameters
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(arbMBMF_boldpat) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(arbMBMF_boldpat) %s is not a recognized parameter name', pair{1})
    end
end

%% Experiment type

switch Exp
    
    case {'Lee2014', 'Lee2014pp2'}
        % Lee2014 subject
        % sbj ids
        IDs = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,23,24};
        % IDs = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24};  

        % SBJstructure
%         sbj_struct = load([root '/A_Research/Dim_control_PFC_metaRL/M3_2014subs_ori/SBJ_structure_each_exp_BESTcollectionJan29']);
        sbj_struct = ... 
            load(fullfile(root, 'human_behavior', 'SBJ_structure_each_exp_BESTcollectionJan29'));
        SBJ = sbj_struct.SBJ{ID};
        
        % Session : 1~2, 1~3, 1~4 or 1~5
        N_sessions = 5*ones(1,24); N_sessions(11)=4; N_sessions(16)=3; N_sessions(20)=2; N_sessions(21)=4;
        N_session = N_sessions(IDs{ID});
        % Session marks
        tmpload = load(fullfile(root, 'fmri_supp', 'sess_mark_Lee2014.mat'));
        sess_mark_Lee2014 = tmpload.sess_mark_Lee2014;
        sess_marks = sess_mark_Lee2014(IDs{ID},:);
        
        
end


%% BOLDpattern

boldpat_dir = [root '/subj_masked_EPI/' ...
    Exp '/' options.roi_type '/' erase(options.prefix, '*') '/' ROI];

fname_format = ['subj_EPI_' Exp '_%s_%s'];

% load BOLD pattern
boldpat_name = [boldpat_dir '/' sprintf(fname_format, ROI, num2str(IDs{ID})) '.mat'];

switch options.SignalType

    case 'raw'
        
        % load BOLD pattern
        subj_struct = load(boldpat_name);     % try to load 'subj' struct
        subj = subj_struct.subj;
        boldpat = get_mat(subj, 'pattern', 'epi');
        
    case 'dt'

        % load BOLD pattern
        subj_struct = load(boldpat_name);     % try to load 'subj' struct
        subj = subj_struct.subj;
        boldpat = get_mat(subj,'pattern','epi_dt');

    case 'z2'

        % load BOLD pattern
        subj_struct = load(boldpat_name);     % try to load 'subj' struct
        subj = subj_struct.subj;
        boldpat = get_mat(subj,'pattern','epi_dt');

    case 'zscore'

        try
            subj_struct = load(boldpat_name);     % try to load 'subj' struct
            subj = subj_struct.subj;
            boldpat = get_mat(subj,'pattern','epi_z');
        catch Err
            disp(Err.message)
            % Activated only for the first run & use saved data in later runs
            disp(['---------- EPI pattern loading for ' ROI ' ----------'])

            boldpat = [];
            for is = 1:N_session
                
                % For the new mask
                % subj = roi_nifti_load_2(exp, sbj_id, session, roi_type, roi_name, varargin)
                subj = roi_nifti_load_2(Exp, ID, is, options.roi_type, ROI, ...
                    'is_cluster', is_cluster, ...
                    'preproc_prefix', options.prefix, ...
                    'mask_num', options.MaskNum, ...
                    'save_dir', []);                

                sel_name = sprintf('sel%d', is);
                which_session = ones(1, size(get_mat(subj,'pattern','epi'), 2));
                subj = initset_object(subj, 'selector', ...
                    sel_name, which_session);
                
                % Detrending
                subj = detrend_pattern(subj, 'epi', sel_name);
                
                % Z scoring
                subj = apply_to_runs(subj, 'epi_dt', sel_name, 'apply_zscore', ... 
                    'new_patname', sprintf('epi_z_%d', is));
                
                ts_z = get_mat(subj, 'pattern', sprintf('epi_z_%d', is));
                boldpat = cat(2, boldpat, ts_z);

            end % session
            subj = init_object(subj,'pattern', 'epi_z');
            subj = set_mat(subj, 'pattern', 'epi_z', boldpat);

            try
                save(boldpat_name, 'subj');
            catch
                mkdir(boldpat_dir)
                save(boldpat_name, 'subj');
            end
            
        end
        % boldpat = get_mat(subj,'pattern','epi');
        % boldpat = get_mat(subj,'pattern','epi_dt');
        % boldpat = get_mat(subj,'pattern','epi_z');
        
    case 'percent'
        boldpat = [];
        for is = 1:N_session
            subj_save_path = ['Z:\JR\dimensionality/subj/' options.roi_type '/' ROI '/subj_' num2str(ID) '_' num2str(is)];
            try
                subj_struct = load(subj_save_path);
                subj = subj_struct.subj;
            catch
                % For the new mask
                % subj = roi_nifti_load_2(exp, sbj_id, session, roi_type, roi_name, varargin)
                subj = roi_nifti_load_2(Exp, ID, is, options.roi_type, ROI, ...
                    'is_cluster', is_cluster, ...
                    'mask_num', options.MaskNum, ...
                    'save_dir', subj_save_path);
            end
            
            % get percent signal
            ts = get_mat(subj,'pattern','epi');
            ts_p = (ts - mean(ts,2)) ./ (mean(ts,2) * ones(1,size(ts,2)));
            
            boldpat = [boldpat, ts_p];
            
        end % session
        
        % voxel screening
        var_ts_p = var(boldpat, [], 2);
        ts_p_thr = prctile(var_ts_p, 95);
        boldpat = boldpat(var_ts_p < ts_p_thr, :) * 100;
        
end


%% Finding fMRI vol index of cue1/cue2/reward/reward_end in each trial (based on HIST_event_info)
TR = options.TR;
delay = options.Delay;

MRI_ids = [];
blk_con = [];
for sess = 1:N_session
%     mri_ids = sess_marks(sess) + SARtime_to_MRIidx(SBJ.HIST_event_info{sess}, TR, delay);
    mri_ids = sess_marks(sess) + ... 
        SARtime_to_MRIidx(SBJ.HIST_event_info{sess}, TR, delay, ... 
        'Exp', Exp, ...
        'strict', options.strict, ...
        'forInterpol', options.interpolation);
    mri_ids(mri_ids>sess_marks(sess+1)) = sess_marks(sess+1);
    MRI_ids = [MRI_ids, mri_ids]; % 8 x N_trials
    blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
end
% 'MRI_ids' expansion (add previous & next trial)
MRI_prev = MRI_ids; MRI_next = MRI_ids;
MRI_prev(:,2:end) = MRI_ids(:,1:(end-1)); MRI_prev(:,1) = nan;
MRI_next(:,1:(end-1)) = MRI_ids(:,2:end); MRI_next(:,end) = nan;
MRI_ids = [MRI_prev; MRI_ids; MRI_next]; % 24 x N_trials

switch options.TrialCondition
    case 'none'
        N_trial = size(MRI_ids,2);
        tri = 1:N_trial;
    case 1 % 1:specific goal, 0:flexible goal
        tri = ismember(blk_con, [1,2]);
    case 2 % 1:flexible goal, 0:specific goal
        tri = ismember(blk_con, [3,4]);
    case 3 % 1:low uncertainty, 0:high uncertainty
        tri = ismember(blk_con, [1,4]);
    case 4 % 1:high uncertainty, 0:low uncertainty
        tri = ismember(blk_con, [2,3]);
end

% possible input events: -7~16, corresponding idx: 1~24
event = options.Event;

if iscell(event)
    % example: event = {[2 3 4], [5 6 7],[8 9]};  
    boldpat_vet = nan * ones(size(boldpat,1), length(event), ...
        size(MRI_ids(:,tri), 2));
    for ei = 1:length(event)
        full_idx = reshape(MRI_ids(event{ei}+8, tri), 1, []);
        null_idx = isnan(full_idx);
        new_idx = full_idx; new_idx(null_idx) = 1;
        boldpat_temp = boldpat(:, new_idx);
        boldpat_temp(:, null_idx) = nan;
        boldpat_temp = reshape(boldpat_temp, ...
        size(boldpat,1), length(event{ei}), []);
%         boldpat_vet(:, ei, :) = mean(boldpat_temp, 2);
        boldpat_vet(:, ei, :) = nanmean(boldpat_temp, 2);
    end
else
    % example: event = 1:8; event = [-1 3 15];  
    % multi-dimensional BOLD pattern (N_voxel x N_event x N_trial)
    full_idx = reshape(MRI_ids(event+8, tri), 1, []);
    null_idx = isnan(full_idx);
    new_idx = full_idx; new_idx(null_idx) = 1; % arbitrary positive integer
    boldpat_vet = boldpat(:, new_idx);
    boldpat_vet(:, null_idx) = nan;
    boldpat_vet = squeeze(reshape(boldpat_vet, ...
        size(boldpat,1), length(event), []));
    N_trial = size(boldpat_vet, 3);
    % disp(size(boldpat_event)); % check
    % boldpat_event = permute(boldpat_event, [1 3 2]); % 200602
end


