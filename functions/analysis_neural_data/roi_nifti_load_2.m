function subj = roi_nifti_load_2(Exp, sbj_id, session, roi_type, roi_name, varargin)
% makes subj structure and load nifti files to creat multi-voxel pattern
% -------------------------------------------------------------------------
% input:
%     exp       (string)    - experiment name (Lee2014)
%     sbj_id    (num)       - subject ID
%     session   (num)       - imaging session
%     roi_type  (string)    - type of ROI (AAL, AICHAmc, sphere10)
%     roi_name  (string)    - name of ROI 
%                 (lilPFC, rilPFC, F3TL, F3TR, V1L, V1R, HIPPOL, HIPPOR,
%                 VentStrL, VentStrR, )
% 
% output: subj (princeton-mvpa-toolbox-master)
% subj fields
%     mask      - ROI mask
%     pattern   - masked multi-voxel pattern
%     ...
% -------------------------------------------------------------------------
% <Name-value pairs>
% is_cluster
% preproc_prefix
% mask_num
% hemisphere
% save_dir


%% Setting & initialization

% default
options = struct('is_cluster', 1, ...
                'preproc_prefix', 'wra', ...
                'mask_num', [], ...
                'hemisphere', [], ...
                'save_dir', []);
option_names = fieldnames(options);
if mod(length(varargin), 2) == 1
    error('(roi_nifti_load) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(roi_nifti_load) %s is not a recognized parameter name', pair{1})
    end
end

% root
if options.is_cluster
    root = '/home/ydsung';
else
    root = '//143.248.73.45/bmlsamba/ydsung';
end

% experiment setting

atlas_dir = fullfile(root, 'fmri_supp', 'atlas');
% atlas_dir = fullfile(root, 'ydsung', 'A_Research', 'princeton-mvpa-toolbox-master', 'ROI_mask', 'atlas');

switch Exp
    
    case 'Lee2014'
        % original subject IDs
        IDs = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,23,24];
        % IDs = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24];
                
        % neuroimaging file (functional EPI volume) directory
        if options.is_cluster
            nifti_dir = [root sprintf('/2014fmri/fmri_arbitration/od-arbitration-%03d', IDs(sbj_id))];
        else
            nifti_dir = ['D:/fmri_arbitration/' sprintf('od-arbitration-%03d', IDs(sbj_id))];
        end
        session_dir = sprintf('/func/run_%04d', session);

    case 'Lee2014pp2'

        % original subject IDs
        IDs = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,23,24];
        % IDs = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24];

        % neuroimaging file (functional EPI volume) directory
        nifti_dir = fullfile(root, 'fmri_Lee2014', sprintf('od-arbitration-%03d', IDs(sbj_id)));
        session_dir = sprintf('/func/run_%04d', session);

        
end

% subj structure initialization
subj = init_subj(Exp, ['subject', num2str(sbj_id)]);


%% load ROI mask (.nii file)

switch roi_type
    
    case 'AAL'
        subj = load_spm_mask_2(subj, [roi_type '_' roi_name], ...
            [atlas_dir '/aal_' Exp '_coord.nii'], ...
            'mask_num', options.mask_num, ...
            'hemisphere', options.hemisphere);
        
    case 'AAL3'
        
        subj = load_spm_mask_2(subj, [roi_type '_' roi_name], ...
            [atlas_dir '/ROI_MNI_V7_' Exp '_coord.nii'], ...
            'mask_num', options.mask_num, ...
            'hemisphere', options.hemisphere);
    
    case 'AICHAmc'
        
        subj = load_spm_mask_2(subj, [roi_type '_' roi_name], ...
            [atlas_dir '/AICHAmc_' Exp '_coord.nii'], ...
            'mask_num', options.mask_num, ...
            'hemisphere', options.hemisphere);

    case 'JuBrain'

        % only for preSMA
        subj = load_spm_mask_2(subj, [roi_type '_' roi_name], ...
            fullfile(atlas_dir, [roi_type '_' roi_name '_' Exp '_coord.nii']));
        
end



%% load fMRI data (.nii file)

% file names
% fileNames = {};
% listing_s = dir(fullfile(nifti_dir, session_dir, [options.preproc_prefix, '*.nii']));
% for it = 1:length(listing_s)
%     tmp_file =listing_s(it).name;
%     if ~isempty(tmp_file)
%         fileNames = cat(1, fileNames, fullfile(nifti_dir, session_dir, tmp_file));
%     end
% end

fileList = dir(fullfile(nifti_dir, session_dir, ... 
    [options.preproc_prefix, '*.nii'])); % (nFiles, 6)
nFiles = size(fileList, 1);
filePaths = cell(nFiles, 1);
for i = 1:nFiles
    filePaths{i} = fullfile(fileList(i).folder, fileList(i).name);
end

% read spm volume using file names
subj = load_spm_pattern(subj, 'epi', [roi_type '_' roi_name], filePaths);


%% save subj structure

if ~isempty(options.save_dir)
    save(options.save_dir, 'subj')
end



end


