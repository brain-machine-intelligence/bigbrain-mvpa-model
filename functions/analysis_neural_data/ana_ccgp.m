
function [] = ana_ccgp(id)

% comment out or set your own path
% addpath(genpath('/home/ydsung/A_Research/princeton-mvpa-toolbox-master')); 

% Setting =================================================================

% ===== Experiment =====
Exp = 'Lee2014'; Exp_sfx = []; prefix = 'wra';


% ===== ROI mask for multi-voxel pattern (MVP) =====

% main ROI set
% ROI = {'F3TL', 'F3TR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR', 'DLPFCL', 'DLPFCR', 'OFC'};

% additional ROIs
% ROI = {'ACCL', 'ACCR', 'SMAL', 'SMAR'};
% ROI = {'preSMAL', 'preSMAR'};

% full ROI set
ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'ACCL', 'ACCR', 'preSMAL', 'preSMAR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};

% AAL only
% ROI = {'DLPFCL', 'DLPFCR', 'OFC'};
% ROI = {'OFC'};


% ====================== labelling by task variable =======================
% ----- task-relevant context mixing -----
LABEL = {'Goal(6,7,8)'}; N_allclass = 3 * ones(1, length(LABEL)); trialSelect = {[]}; 
% LABEL = {'Goal(6,7,8) perm'}; N_allclass = 3 * ones(1, length(LABEL)); trialSelect = {[]};


CClabel_name = 'UncCond'; ccSelect = [];
% CClabel_name = []; ccSelect = [];
categorize = [];
label_pairing = false;

% to categorize continuous variables into several classes
if isempty(categorize)
    categorize = cell(1, length(LABEL));
end
% categorize = cell(1, length(LABEL));  % empty cell for no categorization
% categorize = num2cell(2 * ones(1, length(LABEL)));
% categorize = num2cell(3 * ones(1, length(LABEL)));
% categorize = num2cell(N_allclass);

% lab_shuffle = 1;
lab_shuffle = 0;



%% ========================= MVP setting =========================
% prefix = 'wra*';  % preprocessing prefix
% prefix = 'swar*';
% prefix = 'swa*';
% prefix = 'wa*';

% SignalType = 'percent';
SignalType = 'zscore'; % default setting

% ROI_types = {};
% ROI_types = {'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL'}; % +DLPFC & R/L integrated OFC
ROI_types = {'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL', 'AAL3', 'AAL3', 'JuBrain', 'JuBrain', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'};
% roi_type = 'AAL3';
% roi_type = 'AAL';
% roi_type = 'JuBrain';
roi_type = [];

MaskNum_map = containers.Map({'F3TL', 'F3TR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR', ...
    'DLPFCL', 'DLPFCR', 'OFC', 'ACCL', 'ACCR', 'SMAL', 'SMAR', 'preSMAL', 'preSMAR'}, ...
    {9, ... % F3TL
    10, ... % F3TR
    47, ... % V1L
    48, ... % V1R
    41, ... % HIPPOL
    42, ... % HIPPOR
    157, ... % VentStrL
    158, ... % VentStrR
    7, ... % DLPFCL
    8, ... % DLPFCR
    [5 9 15 27 6 10 16 28], ... % OFC
    [153 155], ... % ACCL
    [154 156], ... % ACCR
    15, ... % SMAL
    16, ... % SMAR
    nan, ... % preSMAL
    nan}); % preSMAR

% EVENT = 5;
% EVENT = 8;
% EVENT = [-4 13];        % for test 201020
% EVENT = {[1 2], [4 5], [7 8]};
% EVENT = {[1 2 4 5 7 8]};
% EVENT = [-1 0];                 % make-up run for F3TR S3
% EVENT = [9 10];               % make-up run for F3TL & F3TR S3
% EVENT = [1 2 4 5 7 8];        % mainly for 'strict' MVP
% EVENT = [1 2 4 5 7 9 10];        
EVENT = [1 2 4 5 7 8 9 10];        % main setting
% EVENT = [-1 0 1 2 4 5 7 8 9 10];
% EVENT = [7 8 9 10];     % make-up run for F3TL blkCond (slave1 too slow)
% EVENT = -3:8;
% EVENT = -3:12;
% EVENT = 1:8;
% EVENT = {[2 3], [5 6], 8};

strict = 0;
% strict = 1;
interpolation = 0;



% ========================= fitclinear setting =========================
% Repeat = 20;    % # of undersample
% Repeat = 50;    % # of undersample
Repeat = 100;    % # of undersample

rRepeat = 0;
% rRepeat = 10;
% rRepeat = 100;

% Coding = 'onevsone';
% Learners = 'svm';
% KFold = 10;
% KFold = 5;
% KFold = 3;
KFold = [];

% CVtype = 'LOOV';
CVtype = 'LOROV';
% CVtype = '10foldCV';

WeightSave = 0;
% WeightSave = 1;

% LOOV = 1;
% LOOV = 0;

alpha = 0.0001;


% ========================= save setting =========================
save_path = 'ccgp_result';

if WeightSave
    basic_name = 'fitclinear_weights';
else
    basic_name = 'fitclinear_accbox';
end
if ~isempty(CClabel_name)
    save_name = [basic_name '_' CClabel_name ccSelect 'CCGP' Exp_sfx];
else
    save_name = [basic_name Exp_sfx];
end



% main ====================================================================

fprintf('%d ================================================================================================ \n', id)


% #########################################################################
% =========================== loop1: label load ===========================
% #########################################################################

labels = [];
trial_select_sfx_cell = cell(1, length(LABEL));
min_class_sizes = zeros(1, length(LABEL));

for labi = 1:length(LABEL)
    
    var_name = erase(LABEL{labi}, ' shuff');
    var_name = erase(var_name, ' LOROV');
    var_name = erase(var_name, ' 5fold');
    var_name = erase(var_name, ' paired');
    var_name = erase(var_name, ' perm');
    
    if contains(var_name, ' run')
        runs_temp = var_name((strfind(var_name, ' run')+3):end);
        runs = nan(1,length(runs_temp));
        for i=1:length(runs); runs(i) = str2double(runs_temp(i)); end
        session = arbMBMF_load_var(Exp, 'Session', id, []);
        selector = 1 * ismember(session, runs);
        var_name = erase(var_name, var_name(strfind(var_name, ' run'):end));
        trial_select_sfx = [];
    elseif ~isempty(trialSelect{labi})
        selector = trialSelect{labi};
        trial_select_sfx = [' ' trialSelect{labi}];
        trial_select_sfx_cell{labi} = trial_select_sfx;
    else
        selector = [];
        trial_select_sfx = [];
    end
    
    % labelling by task variable
    label = arbMBMF_load_var(Exp, var_name, id, selector);
        
    if lab_shuffle || contains(LABEL{labi}, 'shuff')
        valid_lab_idx = ~isnan(label);
        valid_lab_shuffled = label(valid_lab_idx);
        valid_lab_shuffled = valid_lab_shuffled(randperm(length(valid_lab_shuffled)));
        label(valid_lab_idx) = valid_lab_shuffled;
    end
    
    if ~isempty(categorize{labi})
        label = regr2class(label, categorize{labi});
    end
    
    labels = [labels; label];
    [~, class_sizes] = unq_elms(label);
    min_class_sizes(labi) = min(class_sizes);
    
end % labi

% Determining the smallest class size for undersampling
if label_pairing
    
    if contains(CVtype, 'LOROV')
        N_session = unq_elms(arbMBMF_load_var(Exp, 'Session', id, []));
        min_class_size_ratio = (N_session - 1) / N_session;
    else
        min_class_size_ratio = 1;
    end
    
    min_class_size = ceil(min(min_class_sizes) * min_class_size_ratio);
else
    min_class_size = [];
end

% For CCGP: cross-condition labeling for each trial
if ~isempty(CClabel_name)
    CClabel = arbMBMF_load_var(Exp, CClabel_name, id, ccSelect);
else 
    CClabel = [];
end

% #########################################################################
% ============== loop2: decoding label for each ROI & event ===============
% #########################################################################

for roii = 1:length(ROI)
    roi = ROI{roii};
    fprintf([roi ' ============================================= \n'])
    
    if ~isempty(ROI_types)
        roi_type = ROI_types{roii};
    end
    
    if contains(roi, 'us')  % for voxel undersampling analysis
        k = strfind(roi, 'us');
        % temp_roi_name = roi(1:k-numel(num2str(seed_us))-1);
        temp_roi_name = roi(1:k-numel(anaAxis)-1);
        % roi(1:k-numel(num2str(seed_us))-1): roi_name
        % roi(k+2): underset id
        % sub_vox: voxel index subsets
        MVPfull = arbMBMF_boldpat(Exp, temp_roi_name, id, ... % 3D (N_voxel x N_event x N_trial)
            'SignalType', SignalType, ...
            'roi_type', roi_type, ...
            'MaskNum', MaskNum_map(temp_roi_name), ...
            'Event', EVENT, ...
            'strict', strict, ...
            'interpolation', interpolation, ...
            'prefix', prefix);
        MVPfull = MVPfull(sub_vox{str2double(roi(k+2:end))},:,:);
%         disp(sub_vox{str2double(roi(k+2))})
        fprintf('%d ', sub_vox{str2double(roi(k+2:end))})
        disp(' ')
    else
        switch roi
            case 'ilPFC'
                % multi-voxel pattern MVP: [N_voxel x N_event x N_trial]
                lilPFC = arbMBMF_boldpat(Exp, 'lilPFC', id, ...
                    'SignalType', SignalType, ...
                    'roi_type', roi_type, ...
                    'Event', EVENT, ...
                    'strict', strict, ...
                    'interpolation', interpolation);
                rilPFC = arbMBMF_boldpat(Exp, 'rilPFC', id, ...
                    'SignalType', SignalType, ...
                    'roi_type', roi_type, ...
                    'Event', EVENT, ...
                    'strict', strict, ...
                    'interpolation', interpolation);
                MVPfull = cat(1, lilPFC, rilPFC);
                
            case 'V1'
                lV1 = arbMBMF_boldpat(Exp, 'lV1', id, ...
                    'SignalType', SignalType, ...
                    'roi_type', roi_type, ...
                    'Event', EVENT, ...
                    'strict', strict, ...
                    'interpolation', interpolation);
                rV1 = arbMBMF_boldpat(Exp, 'rV1', id, ...
                    'SignalType', SignalType, ...
                    'roi_type', roi_type, ...
                    'Event', EVENT, ...
                    'strict', strict, ...
                    'interpolation', interpolation);
                MVPfull = cat(1, lV1, rV1);
                
            otherwise
                % multi-voxel pattern MVP: [N_voxel x N_event x N_trial]
                MVPfull = arbMBMF_boldpat(Exp, roi, id, ... 
                    'SignalType', SignalType, ...
                    'roi_type', roi_type, ...
                    'MaskNum', MaskNum_map(roi), ...
                    'Event', EVENT, ...
                    'strict', strict, ...
                    'interpolation', interpolation, ...
                    'prefix', prefix);
        end
    end
        
    for labi = 1:length(LABEL)
        
        fprintf([LABEL{labi} ' ================ \n'])
        
        % loaded label
        label = labels(labi, :);
        
        for evi = 1:length(EVENT)
            
            % event specific MVP
            if ndims(MVPfull) > 2
                MVP = squeeze(MVPfull(:,evi,:));
            else
                MVP = MVPfull;
            end
            % event name for saving
            if iscell(EVENT)
                ev_name = EVENT{evi};
            else
                ev_name = EVENT(evi);
            end
            
            % sanity check
            trial_len = size(label, 2);
            if size(MVP,2)~=trial_len
                disp('label-pattern size unmatch'); 
                continue
            end
            unqSet = unq_elms(label);
            if length(unqSet)~=N_allclass(labi) 
                disp('(arbMBMF_shattering) label element error'); 
                continue
            end
            
            % linear shattering
            fprintf('%d ', id)
            if isempty(KFold)
                acc_box = linear_shattering(MVP, label, ...
                    'Exp', Exp, ...
                    'ID', id, ...
                    'CVtype', CVtype, ...
                    'CClabel', CClabel, ...
                    'Repeat', Repeat, ...
                    'rRepeat', rRepeat, ...
                    'MinimalClassSize', min_class_size, ...
                    'WeightSave', WeightSave, ...
                    'Alpha', alpha, ...
                    'AllClassOnly', 1, ...
                    'verbose', 1);
            else
                acc_box = linear_shattering(MVP, label, ...
                    'Kfold', KFold, ...
                    'CClabel', CClabel, ...
                    'Repeat', Repeat, ...
                    'rRepeat', rRepeat, ...
                    'MinimalClassSize', min_class_size, ...
                    'WeightSave', WeightSave, ...
                    'Alpha', alpha, ...
                    'AllClassOnly', 1, ...
                    'verbose', 1);
            end 
            
            % running information
            run_info.Exp = Exp;
            run_info.ROI = ROI;
            run_info.LABEL = LABEL;
            run_info.N_allclass = N_allclass;
            run_info.label_pairing = label_pairing;
            run_info.categorize = categorize;
            run_info.lab_shuffle = lab_shuffle;
            run_info.MVP_setting.SignalType = SignalType;
            run_info.MVP_setting.roi_type = roi_type;
            run_info.MVP_setting.EVENT = EVENT;
            run_info.MVP_setting.strict = strict;
            run_info.MVP_setting.interpolation = interpolation;
            run_info.MVP_setting.prefix = prefix;
            run_info.stat_setting.Repeat = Repeat;
            run_info.stat_setting.rRepeat = rRepeat;
            run_info.stat_setting.CVtype = CVtype;
            run_info.stat_setting.KFold = KFold;
            run_info.stat_setting.WeightSave = WeightSave;
            run_info.stat_setting.Alpha = alpha;
%             run_info.fitcecoc_setting.Coding = Coding;
%             run_info.fitcecoc_setting.Learners = Learners;
%             run_info.fitcecoc_setting.LOOV = LOOV;
                        
            % save information
            save_info.id = id;
            save_info.roi = ROI{roii};
            save_info.roi_type = roi_type;
            save_info.label_name = LABEL{labi};
            save_info.label_true = label;
            save_info.trialSelect = trialSelect{labi};
            save_info.event = ev_name;
            save_info.save_path = save_path;
            save_info.save_name = save_name;
            save_info.acc_box = acc_box;             
            
            % save
            save_dir = [save_path '/' roi_type '/' ROI{roii} '/' , ...
                LABEL{labi} trial_select_sfx_cell{labi} '/' num2str(ev_name) '/'];
            
%             % Rename saved file names (2023-03-30) - run one arbitrary sbj, then change all sbjs
%             file_list = dir([save_dir, '*.mat']);
%             for i = 1:length(file_list)
%                 old_name = file_list(i).name;
%                 % Find the numbers in the filename using regular Expressions
%                 numbers = regExp(old_name, '\d+', 'match');
%                 % Check if there are any numbers in the filename
%                 if ~isempty(numbers)
%                     % Get the first number in the filename
%                     number = numbers{1};
%                     % Change the filename using the number
%                     new_name = [save_name, number, '.mat'];
%                     % Rename the file
%                     movefile([save_dir, old_name], [save_dir, new_name]);
%                 end
%             end
            
            try
                save([save_dir, ...
                    save_name num2str(id) '.mat'], ...
                    'run_info', 'save_info')
                
            catch
                mkdir(save_dir)
                save([save_dir, ...
                    save_name num2str(id) '.mat'], ...
                    'run_info', 'save_info')
            end
            
            disp([save_dir, save_name num2str(id) ' saved'])
            
        end % evi
        
    end % labi
    
end % roii

