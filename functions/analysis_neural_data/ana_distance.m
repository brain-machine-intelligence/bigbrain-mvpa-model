
function [] = ana_distance(id, varargin)

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
save_path = 'distance_result';
save_name = ['distance_' Exp_sfx];  trained_mdl_info = [];

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

            
            % ================ main coding distance analysis ================
            MVP_ori   = MVP;      % original MVP timecourse
            label_ori = label;    % original label timecourse

            runs = arbMBMF_load_var(Exp, 'Session', id, []);
            [runSet, lenRuns] = unq_elms(runs);
            nRuns = length(runSet);

            % Prepare to store coding-distance results
            dists              = nan(1, Repeat);            % scalar distance per "Repeat"
            dists_run_mat      = nan(Repeat, nRuns);        % scalar distance per run & repeat
            pairwiseDist_run   = cell(Repeat, nRuns);       % store the full distance matrix in each run & repeat

            for r = 1:Repeat

                rng(r); % Keep the same randomization approach if needed

                dists_run = nan(1, nRuns);  % distance in each run for this “Repeat” iteration

                for runIdx = 1:nRuns
                    runNum = runSet(runIdx);  % e.g., runSet might be [1, 2, 3, ...]

                    % === Extract data for this run ===
                    [label_run, MVP_run] = lab_pat_undersample( ...
                        label_ori(runs == runNum), MVP_ori(:, runs == runNum) ...
                        );

                    % === Compute mean pattern per class ===
                    [uniqueLabels, ~, idxClass] = unique(label_run);
                    nClass = length(uniqueLabels);

                    meanPatterns = zeros(size(MVP_run, 1), nClass);
                    for c = 1:nClass
                        meanPatterns(:, c) = mean(MVP_run(:, idxClass == c), 2);  % (nVoxel, nClass)
                    end

                    % === Compute full pairwise distance matrix ===
                    % Option 1: "pdist" + "squareform" for Euclidean distances
                    % 'pdist': The distances are arranged in column order in the lower left triangular part of the m×m distance matrix: (2,1), (3,1), ..., (m,1), (3,2), ..., (m,2), ..., (m,m–1).
                    distMat = squareform(pdist(meanPatterns', 'euclidean'));
                    % distMat will be [nClass x nClass]

                    % Store this run's full pairwise distance matrix
                    pairwiseDist_run{r, runIdx} = distMat;

                    % === Compute coding distance ===
                    if nClass == 2
                        % Distance is just that between the 2 class means
                        dists_run(runIdx) = distMat(1,2);
                    else
                        % Distance is the average of pairwise distances among classes
                        upperTriIdx = find(triu(ones(nClass), 1));
                        dists_run(runIdx) = mean(distMat(upperTriIdx));
                    end
                end % for runIdx (loop over runs)

                % Weighted average of distances across runs
                dists(r) = dists_run * lenRuns / sum(lenRuns);
                dists_run_mat(r, :) = dists_run;

            end % for r = 1:Repeat

            % === Save info in the same structure you used before ===
            save_info.dists            = dists;          % distance across "Repeat" runs
            save_info.dists_run_mat    = dists_run_mat;  % distance per run for each "Repeat"
            save_info.pairwiseDist_run = pairwiseDist_run;
            % cell array {Repeat x #Runs}, each cell = full distance matrix among classes


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


