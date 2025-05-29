
function [] = ana_decoding(id, varargin)

% comment out or set your own path
% addpath(genpath('/home/ydsung/A_Research/princeton-mvpa-toolbox-master'));

% Setting =================================================================


% default input parameters
options = struct('ROI', [], ...
                'LABEL', [], ...
                'EVENT', [], ...
                'save_name', []);
% read input parameters
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(ana_decoding) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(ana_decoding) %s is not a recognized parameter name', pair{1})
    end
end

%%
% ========================= Experiment =========================
Exp = 'Lee2014'; Exp_sfx = []; prefix = 'wra';


% ================ ROI mask for multi-voxel pattern (MVP) =================
if ~isempty(options.ROI)
    ROI = options.ROI;
else
    
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

end


% ====================== labelling by task variable =======================
if ~isempty(options.LABEL)
    LABEL = options.LABEL;
    LABEL_NULL = {};    % 2021-11-11 temporary
else
    
%
% LABEL_NULL = {'S3'};
LABEL_NULL = {};


LABEL = {'UncCond', 'UncCond'}; N_allclass = [2 2]; trialSelect = {'G', 'H'}; lab_shuffle = 0; % for fig.2c & d
% LABEL = {'UncCond perm', 'UncCond perm'}; N_allclass = [2 2]; trialSelect = {'G', 'H'}; lab_shuffle = 0;
% LABEL = {'Goal(6,7,8)', 'S2', 'S3'}; N_allclass = [3 4 4]; trialSelect = {[], 'G', 'G'}; lab_shuffle = 0; % for fig.2a & b


end


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

if ~isempty(options.EVENT)
    EVENT = options.EVENT;
else
    
% EVENT = {[1 2], [4 5], [7 8]};
EVENT = [1 2 4 5 7 8 9 10];      % mainly for 'strict' MVP
% EVENT = [9 10];      % mainly for 'strict' MVP
% EVENT = 4;      % 2021-07-07 test for weight map
% EVENT = [1 2 5 7 8 9 10];      % 2021-07-07 test for weight map
% EVENT = -3:8;
% EVENT = -3:12;
% EVENT = 1:8;
% EVENT = {[2 3], [5 6], 8};

end

strict = 0;
% strict = 1;
interpolation = 0;



% ========================= fitcecoc setting =========================
% Repeat = 20;    % # of undersample
Repeat = 100;    % # of undersample

rRepeat = 0;
% rRepeat = 5;
% rRepeat = 10;
% rRepeat = 100;

% Coding = 'onevsone';
Coding = 'onevsall';
% ---------------------------------
% for fitrlinear
% Learners = 'svm';
% for fitcecoc
Learners = templateLinear('Learner', 'svm', 'Regularization', 'ridge', 'Solver', 'dual');
% ---------------------------------
% KFold = 10;
% KFold = 5;
KFold = [];

% CVtype = 'LOOV';
CVtype = 'LOROV';

% % LOOV = 1;
% LOOV = 0;



% ========================= save setting =========================
save_path = 'decoding_result';

if ~isempty(options.save_name)
    save_name = options.save_name;
%     trained_mdl_path = [];
    trained_mdl_info = [];
else
    
% save_name = 'CVMdl_';     trained_mdl_info = [];
% save_name = ['fitrlinear_corrcoef_' Exp_sfx];  trained_mdl_info = [];
save_name = ['fitcecoc_CVacc_' Exp_sfx];  trained_mdl_info = [];
% save_name = ['fitcecoc_CVacc_PredLab_ova_' Exp_sfx];  trained_mdl_info = [];
% save_name = 'fitcecoc_CVMdl_';    trained_mdl_info = [];
% save_name = ['fitcecoc_weights_' Exp_sfx];  trained_mdl_info = [];

% --- for CCGP ---
trained_label = [];

% 2021-04-30
trained_roi = [];
trained_event = [];

% trained_mdl_info = [trained_roi trained_label num2str(trained_event)];
% save_name = ['fitcecoc_CCGP_' Exp_sfx trained_mdl_info];

end

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
                    disp('(ana_decoding) label element error');
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
            run_info.stat_setting.Coding = Coding;
            run_info.stat_setting.Learners = Learners;
            % run_info.stat_setting.LOOV = LOOV;
            run_info.stat_setting.CVtype = CVtype;
            run_info.stat_setting.KFold = KFold;
            
            % save information
            save_info.id = id;
            save_info.roi = ROI{roii};
            save_info.roi_type = roi_type;
            save_info.label_name = LABEL{labi};
            save_info.event = ev_name;
            save_info.save_path = save_path;
            save_info.save_name = save_name;
                        
            % ================ main decoding ================
            MVP_ori = MVP;
            label_ori = label;
            
            if isempty(trained_mdl_info)
                temp_name = save_name;
            else
                temp_name = erase(save_name, trained_mdl_info);
            end
            
            if contains(temp_name, 'fitcecoc_CVacc_') || contains(temp_name, 'fitcecoc_weights_')
                
                if contains(temp_name, 'fitcecoc_weights_')
%                     betas_r = cell(1, Repeat);
                    beta_foldmeans_r = nan * ones(size(MVP, 1), Repeat);
%                     betas = nan * ones(size(MVP, 1), Repeat, KFold);
                end
                
                accs = nan * ones(1, Repeat);
                accs_r = nan * ones(Repeat, rRepeat);
                accs_run_mat = nan * ones(Repeat, length(unq_elms(arbMBMF_load_var(Exp, 'Session', id, []))));
                pred_labs = nan(Repeat, length(label_ori));
                for r = 1:Repeat
                    
                    % rng(r);
                    rng(r*423);
                    % rng(r+421);
                    % rng(2024 + r);
                    % rng('shuffle');
                    % undersampling
                    [label, MVP] = lab_pat_undersample(label_ori, MVP_ori);
                    
                    if length(label) < 2 * KFold
                        disp('insufficient samples for CV')
                        continue
                        
                    else
                        
                        % ===== linear multiclass classification =====
                        if isempty(KFold) || (KFold==0)
                            switch CVtype
                                case 'LOOV'     % leave-one-out validation
                                    [CVMdl, Param] = fitcecoc(MVP, label, ...
                                        'Coding', Coding, ...
                                        'Learners', Learners, ...
                                        'ObservationsIn', 'columns', ...
                                        'Leaveout', 'on');
                                    % generalized accuracy (CV acc)
                                    accs(r) = 1 - kfoldLoss(CVMdl);                        

                                case 'LOROV'    % leave-one-run-out validation
                                    
                                    runs = arbMBMF_load_var(Exp, 'Session', id, []);
                                    [runSet, lenRuns] = unq_elms(runs);
                                    accs_run = nan(1, length(runSet));
                                    accs_r_run = nan(rRepeat, length(runSet));
                                    pred_lab = nan(size(label_ori));

                                    for run = runSet
                                        [label_test, MVP_test] = ... 
                                            lab_pat_undersample(label_ori(runs == run), MVP_ori(:, runs == run));
                                        [label_train, MVP_train] = ... 
                                            lab_pat_undersample(label_ori(runs ~= run), MVP_ori(:, runs ~= run));
                                        Mdl = fitcecoc(MVP_train, label_train, ...
                                            'Coding', Coding, ...
                                            'Learners', Learners, ...
                                            'ObservationsIn', 'columns');
                                        accs_run(run) = 1 - loss(Mdl, MVP_test', label_test');
                                        test_idx = (runs == run) & (~isnan(label_ori));
                                        pred_lab(test_idx) = predict(Mdl, MVP_ori(:, test_idx)');
                                        % clear Mdl

                                        if rRepeat
                                            accs_rr = nan*ones(rRepeat, 1);
                                            for rr = 1:rRepeat
                                                labTrainPerm = label_train(randperm(length(label_train)));
                                                % disp(labTrainPerm(1:20))
                                                labTestPerm = label_test(randperm(length(label_test)));
                                                MdlR = fitcecoc(MVP_train, labTrainPerm, ...
                                                    'Coding', Coding, ...
                                                    'Learners', Learners, ...
                                                    'ObservationsIn', 'columns');
                                                accs_rr(rr) = 1 - loss(MdlR, MVP_test', labTestPerm');
                                            end
                                            accs_r_run(:, run) = accs_rr; % (accs_r_run: rRepeat, nRun)
                                        else
                                            accs_r_run(:, run) = repmat(1/N_allclass(labi), rRepeat, 1);
                                        end

                                    end % run loop
                                    accs(r) = accs_run * lenRuns / sum(lenRuns); % (1, nRun) x (nRun, 1) / scalar
                                    accs_run_mat(r,:) = accs_run;
                                    accs_r(r, :) = ( accs_r_run * lenRuns / sum(lenRuns) )'; % (RHS) ( (rRepeat, nRun) x (nRun, 1) / scalar )'
                                    % sanity check
                                    % if any(isnan(pred_lab)); error('error in pred_lab'); end
                                    if sum(~isnan(pred_lab)) ~= sum(~isnan(label_ori))
                                        error('(ana_decoding) error in pred_lab')
                                    end
                                    pred_labs(r, :) = pred_lab;
                            end

                        else % CV
                            [CVMdl, Param] = fitcecoc(MVP, label, ...
                                'Coding', Coding, ...
                                'Learners', Learners, ...
                                'ObservationsIn', 'columns', ...
                                'KFold', KFold);
                            % generalized accuracy (CV acc)
                            accs(r) = 1 - kfoldLoss(CVMdl);
                            if contains(temp_name, 'fitcecoc_weights_')
                                betas_fold = nan * ones(size(MVP, 1), KFold);
                                for foldi = 1:KFold
%                                     betas(:,r,foldi) = CVMdl.Trained{foldi}.BinaryLearners{:}.Beta;
                                    betas_fold(:,foldi) = CVMdl.Trained{foldi}.BinaryLearners{:}.Beta;
                                end
%                                 betas_r{r} = betas_fold;
                                beta_foldmeans_r(:,r) = mean(betas_fold, 2);
                            end
                
                        end
                        
                    end % if length(label) < 2 * KFold
                    
                    %                 % save
                    %                 save([save_path '/' roi_type '/' ROI{roii} '/' , ...
                    %                     LABEL{labi} '/' num2str(EVENT(evi)) '/', ...
                    %                     save_name num2str(id) '.mat'], ...
                    %                     'CVMdl', 'run_info', 'save_info')
                    
                    clear CVMdl
                    
                end % Repeat
                
                % save information
                save_info.accs = accs;
                save_info.accs_r = accs_r;
                if (isempty(KFold) || (KFold==0)) && strcmp(CVtype, 'LOROV')
                    save_info.accs_run_mat = accs_run_mat;
                    save_info.pred_labs = pred_labs;
                end
                if contains(temp_name, 'fitcecoc_weights_')
                    save_info.beta_foldmeans_r = beta_foldmeans_r;
                end
                
            elseif contains(temp_name, 'fitcecoc_CVMdl_')
                
                % undersampling
                rng('shuffle');
                [label, MVP] = lab_pat_undersample(label_ori, MVP_ori);
                
                % linear multiclass classification
                [CVMdl, Param] = fitcecoc(MVP, label, ...
                    'Coding', Coding, ...
                    'Learners', Learners, ...
                    'ObservationsIn', 'columns');
                
                % save information
                save_info.CVMdl = CVMdl;
                save_info.Param = Param;
                
                clear CVMdl
                
            elseif contains(temp_name, 'fitcecoc_CCGP_')
                
                if isempty(trained_mdl_info)
                    try
                        %                             mdl_load = load([save_path '/' roi_type '/' ROI{roii} '/' , ...
                        %                                 LABEL{labi} '/' num2str(ev_name) '/', ...
                        %                                 save_name num2str(id) '.mat'], 'save_info');
                        mdl_load = load([save_path '/' roi_type '/' ROI{roii} '/' , ...
                            trained_label '/' num2str(ev_name) '/', ...
                            'fitcecoc_CVMdl_' Exp_sfx num2str(id) '.mat'], 'save_info');
                        mdl = mdl_load.save_info.CVMdl;
                    catch
                        ana_decoding(id, ...
                            'ROI', ROI(roii), ...
                            'LABEL', {trained_label}, ...
                            'EVENT', EVENT(evi), ...
                            'save_name', ['fitcecoc_CVMdl_' Exp_sfx]);
                        %                             mdl_load = load([save_path '/' roi_type '/' ROI{roii} '/' , ...
                        %                                 LABEL{labi} '/' num2str(ev_name) '/', ...
                        %                                 save_name num2str(id) '.mat'], 'save_info');
                        mdl_load = load([save_path '/' roi_type '/' ROI{roii} '/' , ...
                            trained_label '/' num2str(ev_name) '/', ...
                            'fitcecoc_CVMdl_' Exp_sfx num2str(id) '.mat'], 'save_info');
                        mdl = mdl_load.save_info.CVMdl;
                    end
                    
                else
                    
                    if isempty(trained_roi); mdlroi = ROI{roii}; else; mdlroi = trained_roi; end
                    if isempty(trained_label); mdllab = LABEL{labi}; else; mdllab = trained_label; end
                    if isempty(trained_event); mdleve = EVENT(evi); else; mdleve = trained_event; end
                    
                    try
                        mdl_load = load([save_path '/' roi_type '/' mdlroi '/' , ...
                            mdllab '/' num2str(mdleve) '/', ...
                            'fitcecoc_CVMdl_' Exp_sfx num2str(id) '.mat'], 'save_info');
                        mdl = mdl_load.save_info.CVMdl;
                    catch
                        ana_decoding(id, ...
                            'ROI', {mdlroi}, ...
                            'LABEL', {mdllab}, ...
                            'EVENT', mdleve, ...
                            'save_name', ['fitcecoc_CVMdl_' Exp_sfx]);
                        mdl_load = load([save_path '/' roi_type '/' mdlroi '/' , ...
                            mdllab '/' num2str(mdleve) '/', ...
                            'fitcecoc_CVMdl_' Exp_sfx num2str(id) '.mat'], ...
                            'save_info');
                        mdl = mdl_load.save_info.CVMdl;
                    end
                    
                end
                
                CCpredict = predict(mdl, MVP_ori')';
                CCGP = 1 - loss(mdl, MVP_ori', label_ori');
                
                % save information
                
                save_info.CCpredict = CCpredict;
                save_info.CCGP = CCGP;
                
                clear mdl
                
            elseif contains(temp_name, 'fitrlinear_corrcoef')
                
                accs = nan * ones(1, Repeat);
                for r = 1:Repeat

                    % ===== linear regression =====
                    CVMdl = fitrlinear(MVP, label, ...
                        'Learner', Learners, ...
                        'ObservationsIn', 'columns', ...
                        'Regularization', 'ridge', ...
                        'Solver', 'dual', ...
                        'KFold', KFold);
                    % Learners = templateLinear('Learner', 'svm', 'Regularization', 'ridge', 'Solver', 'dual');
                    % generalized accuracy (CV acc)
                    tempp = corrcoef(kfoldPredict(CVMdl), CVMdl.Y);
                    accs(r) = tempp(1,2);
                    
                end % Repeat
                save_info.accs = accs;
                
            else
                error('(ana_decoding) save_name error')
            end
                    
            save_info.label_true = label_ori;

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


