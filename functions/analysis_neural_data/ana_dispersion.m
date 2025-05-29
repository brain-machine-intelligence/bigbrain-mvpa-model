
function [] = ana_dispersion(id, varargin)

% comment out or set your own path
% addpath(genpath('/home/ydsung/A_Research/princeton-mvpa-toolbox-master'));

%%
% ========================= Experiment =========================
Exp = 'Lee2014'; Exp_sfx = []; prefix = 'wra';
% Exp = 'Lee2014pp2'; Exp_sfx = '_Lee2014pp2'; prefix = 'wra';
% Exp = 'Kim2019'; Exp_sfx = '_Kim2019_wa'; prefix = 'wa';
% Exp = 'Heo2021'; Exp_sfx = '_Heo2021_wra'; prefix = 'wra';
% Exp = 'Heo2021'; Exp_sfx = '_Heo2021_swra'; prefix = 'swra';
% Exp = 'Shin2025'; Exp_sfx = '_Shin2025_swar'; prefix = 'swar';


% ================ ROI mask for multi-voxel pattern (MVP) =================
% main ROI set
% ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};

% additional ROIs
% ROI = {'ACCL', 'ACCR', 'SMAL', 'SMAR'};
% ROI = {'ACCL', 'ACCR', 'preSMAL', 'preSMAR'};

% full ROI set
ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'ACCL', 'ACCR', 'preSMAL', 'preSMAR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};

% AAL only
% ROI = {'DLPFCL', 'DLPFCR', 'OFC'};
% ROI = {'OFC'};

% test
% ROI = {'HIPPOL', 'HIPPOR'};



% ====================== labelling by task variable =======================

% LABEL_NULL = {'S3'};
LABEL_NULL = {};


% LABEL = {'UncCond'}; N_allclass = 2; trialSelect = {[]}; lab_shuffle = 0; % for fig.2c & d
LABEL = {'UncCond', 'UncCond'}; N_allclass = [2 2]; trialSelect = {'G', 'H'}; lab_shuffle = 0; % for fig.2c & d
% LABEL = {'UncCond', 'UncCond'}; N_allclass = [2 2]; trialSelect = {'cL', 'cH'}; lab_shuffle = 0;
% LABEL = {'UncCond', 'UncCond'}; N_allclass = [2 2]; trialSelect = {'binPMB1', 'binPMB0'}; lab_shuffle = 0;
% LABEL = {'UncCondd1', 'UncCondd2', 'UncCondd1', 'UncCondd2', ...
%     'UncConda1', 'UncConda2', 'UncConda1', 'UncConda2'}; N_allclass = 2 * ones(size(LABEL)); trialSelect = {'G', 'G', 'H', 'H', 'G', 'G', 'H', 'H'}; lab_shuffle = 0;
% LABEL = {'UncCond delay1', 'UncCond delay2', 'UncCond delay1', 'UncCond delay2', ...
%     'UncCond advance1', 'UncCond advance2', 'UncCond advance1', 'UncCond advance2'}; N_allclass = 2 * ones(size(LABEL)); trialSelect = {'G', 'G', 'H', 'H', 'G', 'G', 'H', 'H'}; lab_shuffle = 0;
% LABEL = {'UncCond delay1', 'UncCond delay2', 'UncCond delay1', 'UncCond delay2'}; N_allclass = [2 2 2 2]; trialSelect = {'G', 'G', 'H', 'H'}; lab_shuffle = 0;
% LABEL = {'UncCond advance1', 'UncCond advance2', 'UncCond advance1', 'UncCond advance2'}; N_allclass = [2 2 2 2]; trialSelect = {'G', 'G', 'H', 'H'}; lab_shuffle = 0;
% LABEL = {'UncCond delay1', 'UncCond delay2', 'UncCond delay1', 'UncCond delay2'}; N_allclass = [2 2 2 2]; trialSelect = {'cL', 'cL', 'cH', 'cH'}; lab_shuffle = 0;
% LABEL = {'UncCond shuff', 'UncCond shuff'}; N_allclass = [2 2]; trialSelect = {'G', 'H'}; lab_shuffle = 1;
% LABEL = {'UncCond perm', 'UncCond perm'}; N_allclass = [2 2]; trialSelect = {'G', 'H'}; lab_shuffle = 0;
% LABEL = {'Goal(6,7,8)', 'S2', 'S3'}; N_allclass = [3 4 4]; trialSelect = {[], 'G', 'G'}; lab_shuffle = 0; % for fig.2a & b


% to categorize continuous variables into several classes
categorize = cell(1, length(LABEL));  % empty cell for no categorization
% categorize = num2cell(2 * ones(1, length(LABEL)));
% categorize = num2cell(4 * ones(1, length(LABEL)));
% categorize = num2cell(5 * ones(1, length(LABEL)));



% ========================= MVP setting =========================
% prefix = 'wra*';  % preprocessing prefix
% prefix = 'swa*';

% SignalType = 'percent';
SignalType = 'zscore'; % default setting

% ROI_types = {};
% (main ROI set)
% ROI_types = {'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'}; % +DLPFC & R/L integrated OFC
% (full ROI set)
ROI_types = {'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL', 'AAL3', 'AAL3', 'JuBrain', 'JuBrain', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'};
% (additional ROIs)
% ROI_types = {'AAL3', 'AAL3', 'JuBrain', 'JuBrain'};

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

% EVENT = {[1 2], [4 5], [7 8]};
EVENT = [1 2 4 5 7 8 9 10];      % mainly for 'strict' MVP
% EVENT = [9 10];      % mainly for 'strict' MVP
% EVENT = 4;      % 2021-07-07 test for weight map
% EVENT = [1 2 5 7 8 9 10];      % 2021-07-07 test for weight map
% EVENT = -3:8;
% EVENT = -3:12;
% EVENT = 1:8;
% EVENT = {[2 3], [5 6], 8};

strict = 0;
% strict = 1;
interpolation = 0;



% ========================= Undersampling setting =========================
% Repeat = 1;    % # of undersample
% Repeat = 5;    % # of undersample
% Repeat = 10;    % # of undersample
Repeat = 20;    % # of undersample
% Repeat = 100;    % # of undersample

rRepeat = 0;
% rRepeat = 5;
% rRepeat = 10;
% rRepeat = 100;


% ========================= save setting =========================
save_path = 'dispersion_result';
save_name = ['dispersion_' Exp_sfx];  trained_mdl_info = [];

% --- for CCGP ---
trained_label = [];
% trained_label = 'Goal(6,7,8) uH';
% trained_label = 'Goal(6,7,8) uH LOROV';
% trained_label = 'Goal(6,7,8) uL';
% trained_label = 'Goal(6,7,8) uL LOROV';
% trained_label = 'Goal uH LOROV';
% trained_label = 'Goal uL LOROV';
% trained_label = 'Goal cH LOROV';
% trained_label = 'Goal cL LOROV';
% trained_label = 'Unc ChoOpt 1';

% % 2021-03-24
% trained_roi = [];
% trained_label = 'S3';
% trained_event = 8;

% 2021-04-30
trained_roi = [];
% trained_label = 'S3';
% trained_label = 'S3 run24';
% trained_label = 'S3 run13';
% trained_label = 'S3 run234';
% trained_label = 'S3 run134';
% trained_label = 'S3 run124';
% trained_label = 'S3 run123';
% trained_label = 'S3(6,9)_H';
% trained_label = 'S3(6,7,9)_H';
% trained_label = 'S3(6,7,9)_H LOROV';
% trained_label = 'S3(6,9)_G';
% trained_label = 'S3(6,7,9)_G';
% trained_label = 'S3(6,7,9)_G LOROV';
% trained_label = 'UncCond';
% trained_label = 'GoalCond uL';
% trained_label = 'GoalCond uH';
% trained_label = 'UncCond H';
% trained_label = 'UncCond G';
% trained_label = 'GoalCond P3';
% trained_label = 'GoalCond P1';
% trained_label = 'GoalCond P2';
% trained_label = 'UncCond P1';
% trained_label = 'UncCond P2';
% trained_label = 'UncCond P3';
% trained_label = 'CmplxCond P1';
% trained_label = 'CmplxCond P2';
% trained_label = 'CmplxCond P3';
% trained_label = 'prevUC P3';
% trained_label = 'Goal(6,7,8)';
trained_event = [];
% trained_event = 8;  % 2021-05-31

% trained_mdl_info = [trained_roi trained_label num2str(trained_event)];
% save_name = ['fitcecoc_CCGP_' Exp_sfx trained_mdl_info];


% main ====================================================================

fprintf('%d ================================================ \n', id)

% ===== loop1: label load =====
labels = [];

for labi = 1:length(LABEL)
    
    if contains(LABEL{labi}, 'delay')
        var_name = split(LABEL{labi}, ' delay');
        var_name = var_name{1};
    elseif contains(LABEL{labi}, 'advance')
        var_name = split(LABEL{labi}, ' advance');
        var_name = var_name{1};
    else
        var_name = LABEL{labi};
    end    
    var_name = erase(var_name, ' shuff');
    var_name = erase(var_name, ' LOROV');
    var_name = erase(var_name, ' 5fold');
    var_name = erase(var_name, ' perm');
    
    if contains(var_name, ' run')
        runs_temp = var_name((strfind(var_name, ' run')+3):end);    % 'run keyword must be the last
        runs = nan(1,length(runs_temp));
        for i=1:length(runs); runs(i) = str2double(runs_temp(i)); end
        session = arbMBMF_load_var(Exp, 'Session', id, []);
        selector = 1 * ismember(session, runs);
        var_name = erase(var_name, var_name(strfind(var_name, ' run'):end));
    elseif ~isempty(trialSelect{labi})
        selector = trialSelect{labi};
        trial_select_sfx = [' ' trialSelect{labi}];
        trial_select_sfx_cell{labi} = trial_select_sfx;
    else
        selector = [];
        trial_select_sfx = [];
        trial_select_sfx_cell{labi} = trial_select_sfx;
    end
   
    label = arbMBMF_load_var(Exp, var_name, id, selector);

    if contains(LABEL{labi}, ' delay')
        tokens = regexp(LABEL{labi}, 'delay\s*(\d+)', 'tokens');
        nDelay = str2double(tokens{1}{1}); % # of trials to delay the label
        % labelOri = label; % only used for test
        label = [nan(1, nDelay), label( 1:(end-nDelay) )];
        % imagesc([labelOri; label]); % only used for test
    end

    if contains(LABEL{labi}, ' advance')
        tokens = regexp(LABEL{labi}, 'advance\s*(\d+)', 'tokens');
        nAdvance = str2double(tokens{1}{1}); % # of trials to advance the label
        % labelOri = label; % only used for test
        label = [label( (1+nAdvance):end ), nan(1, nAdvance)];
        % imagesc([labelOri; label]); % only used for test
    end
    
    if lab_shuffle || contains(LABEL{labi}, 'shuff')
        valid_lab_idx = ~isnan(label);
        valid_lab_shuffled = label(valid_lab_idx);
        valid_lab_shuffled = valid_lab_shuffled(randperm(length(valid_lab_shuffled)));
        label(valid_lab_idx) = valid_lab_shuffled;
    end
    
    if ~isempty(categorize{labi})
        label = regr2class(label, categorize{labi});
    end

    labels = cat(1, labels, label);

end % labi

% load label to be nullified
if ~isempty(LABEL_NULL)
    label_nulls = [];
    for ni = 1:length(LABEL_NULL)
        %     label_null = arbMBMF_load_var(Exp, var_name, id, selector);
        label_null = arbMBMF_load_var(Exp, LABEL_NULL{ni}, id, []);
        label_nulls = [label_nulls; label_null];
    end
    save_name = [save_name, 'null', LABEL_NULL{:} '_'];
end


% ===== loop2: decoding label for each ROI & event =====

for roii = 1:length(ROI)
    roi = ROI{roii};
    fprintf([roi ' ================================ \n'])
    
    if ~isempty(ROI_types)
        roi_type = ROI_types{roii};
    end


    % multi-voxel pattern MVP: [N_voxel x N_event x N_trial]
    MVPfull = arbMBMF_boldpat(Exp, roi, id, ...
        'SignalType', SignalType, ...
        'roi_type', roi_type, ...
        'MaskNum', MaskNum_map(roi), ...
        'Event', EVENT, ...
        'strict', strict, ...
        'interpolation', interpolation, ...
        'prefix', prefix);


    % removing mean neural acvitity evoked by null-label
    if ~isempty(LABEL_NULL)
        for ni = 1:length(LABEL_NULL)
            MVPfull = lab_pat_nullify(label_nulls(ni, :), MVPfull, ...
                'NullDim', 3);
        end
    end
    
    for labi = 1:length(LABEL)
        
        fprintf([LABEL{labi} ' ================ \n'])
        
        for evi = 1:length(EVENT)
            
            % loaded label 2021-03-21 - under the event loop for some reason
            label = labels(labi, :);
            
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
            
            if ~contains(save_name, 'fitrlinear')
                unqSet = unq_elms(label);
                if length(unqSet)~=N_allclass(labi)
                    disp('(arbMBMF_decoding) label element error');
                    continue
                end
            end
            
            % running information
            run_info.Exp = Exp;
            run_info.ROI = ROI;
            run_info.LABEL = LABEL;
            run_info.N_Allclass = N_allclass;
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
            
            % save information
            save_info.id = id;
            save_info.roi = ROI{roii};
            save_info.roi_type = roi_type;
            save_info.label_name = LABEL{labi};
            save_info.event = ev_name;
            save_info.save_path = save_path;
            save_info.save_name = save_name;


            % ================ within-class dispersion analysis ================

            MVP_ori   = MVP;          % original MVP timecourse
            label_ori = label;        % original label timecourse

            % Identify all runs & their lengths
            runs = arbMBMF_load_var(Exp, 'Session', id, []);
            [runSet, lenRuns] = unq_elms(runs);
            nRuns = length(runSet);

            % Identify the full set of classes present in 'label_ori'
            allClasses = unique(label_ori(~isnan(label_ori)));
            nClasses   = length(allClasses);

            % We will store run-level dispersions in a 3D matrix:
            % disp_run_mat(r, runIdx, c)
            %   = dispersion for class c in run "runIdx" for the r-th repetition
            disp_run_mat = nan(Repeat, nRuns, nClasses);
            var_run_mat = nan(Repeat, nRuns, nClasses);

            for rpi = 1:Repeat
                rng(rpi);

                for runi = 1:nRuns
                    runNum = runSet(runi);

                    % --- Extract & undersample data for this run ---
                    maskRun = (runs == runNum);
                    [label_run, MVP_run] = lab_pat_undersample(label_ori(maskRun), MVP_ori(:, maskRun));

                    % For each class c, compute average distance from class centroid
                    % -------------------------------------------------------------
                    [uniqueRunClasses, ~, idxClass] = unique(label_run);  % classes in this run
                    nRunClasses = length(uniqueRunClasses);

                    % Pre-allocate for this run (some classes might not appear in every run)
                    runDisp = nan(1, nClasses);  % eventually fill only the relevant classes
                    runVar = nan(1, nClasses);

                    for cIdx = 1:nRunClasses
                        thisClass = uniqueRunClasses(cIdx);

                        % Find columns for this class
                        maskC = (idxClass == cIdx);
                        X_c   = MVP_run(:, maskC);  % data points for this class
                        varVox = var(X_c, 0, 2); % (nVox, 1)

                        % Mean pattern of class c
                        mu_c = mean(X_c, 2);

                        % Distances of each sample from mu_c
                        distVec = sqrt( sum( (X_c - mu_c).^2, 1 ) );
                        % distVec is [1 x #samples_in_class]

                        % Average distance
                        avgDist = mean(distVec);
                        avgVar = mean(varVox);

                        % Identify the global index of this class among allClasses
                        globalClassIdx = find(allClasses == thisClass);

                        % Store in runDisp
                        runDisp(globalClassIdx) = avgDist;
                        runVar(globalClassIdx) = avgVar;
                    end

                    % Store run-level dispersion for this repetition
                    disp_run_mat(rpi, runi, :) = runDisp;
                    var_run_mat(rpi, runi, :) = runVar;
                end % for runIdx
            end % for r

            % --- Now compute the weighted average of run-level dispersions for each repeat ---
            % disp_repeat(r, c) = weighted average (by run length) of run-level dispersion for class c
            disp_repeat = nan(Repeat, nClasses);
            var_repeat = nan(Repeat, nClasses);

            for rpi = 1:Repeat
                for c = 1:nClasses
                    classRunValues1 = squeeze(disp_run_mat(rpi, :, c)); % [1 x nRuns]
                    classRunValues2 = squeeze(var_run_mat(rpi, :, c)); % [1 x nRuns]
                    % If a class was missing in a particular run, you'll have NaN there.
                    % We only weight over runs where it's non-NaN, for example.

                    validMask = ~isnan(classRunValues1);
                    if ~any(validMask)
                        % If no run had data for this class, keep NaN
                        continue
                    end
                    if ~all(validMask == ~isnan(classRunValues2))
                        error('(ana_dispersion) invalid measures')
                    end
                    
                    % Weighted average across valid runs:
                    w = lenRuns(validMask);             % run lengths for valid runs (col)
                    v1 = classRunValues1(validMask);      % dispersion values for valid runs (row)
                    v2 = classRunValues2(validMask);      % dispersion values for valid runs (row)
                    disp_repeat(rpi, c) = (v1 * w) / sum(w);
                    var_repeat(rpi, c) = (v2 * w) / sum(w);
                end
            end

            % --- Final: produce a single measure per class across repeats (optional) ---
            % For example, you can average across repeats:
            Dispersion = mean(disp_repeat, 1, 'omitnan')';  % (nClasses x 1) column vector
            Variance = mean(var_repeat, 1, 'omitnan')';  % (nClasses x 1) column vector

            % === Save the results in your "save_info" structure, just like you did before ===
            save_info.disp_repeat      = disp_repeat;       % (Repeat x nClasses)
            save_info.Dispersion       = Dispersion;        % (nClasses x 1) final average
            save_info.var_repeat       = var_repeat;        % (Repeat x nClasses)
            save_info.Variance         = Variance;          % (nClasses x 1) final average
            save_info.allClasses       = allClasses;        % for reference

            % save
            save_dir = [save_path '/' roi_type '/' ROI{roii} '/' , ...
                LABEL{labi} trial_select_sfx_cell{labi} '/' num2str(ev_name) '/'];


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


