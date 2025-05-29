
function [] = BATCH_ORI_FULL(parID, nSimul, humanSimul, mbmf, saveDirPrefix)

    disp(datetime('now')) % For running time check by the log files

    % t = datetime('now');
    % t = clock;
    LIST_SBJ={'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
        'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};
    list_sbj={LIST_SBJ{[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,23,24]}};
    list_sbj= { list_sbj{parID:parID} };
%     parpool('local', maxcore);
    
    % ===== Save info =====
    save_dir = [saveDirPrefix num2str(nSimul)]; 
    if all(humanSimul)
        fit_or_gen = 1; % parameter fitting to human behavior data
    else
        fit_or_gen = 2; % agent episode generation
    end
    save_name = 'SBJ';

    %% Parameter range setting ('param_init')
    % original parameter description
    % Description : param1 = zero prediction error threshold of SPE
    %               param2 = learning rate for the estimate of absolute reward
    %                           prediction error
    %               param3 = the amplitude of a transition rate function(mb-mf
    %               param4 = the amplitude fo a transition rate function(mf-mb
    %               param5 = inverse softmax temparature
    %               param6 = learning rate of the model based and the model
    %                      free (both of them use same value of learning rate)

    mode.param_BoundL = [0.3 0.01  0.1 0.1 0.01 0.01];
    mode.param_BoundU = [0.8 0.35 10 10 0.5 0.2];


    %% Simulation mode setting
    mode.param_length = 6 ;%size(param_init,2);
    mode.total_simul = 1; % number of simulation iteration
    
    mode.experience_sbj_events = [humanSimul humanSimul]; % subject episode replication vs virtual simulation in eval_ArbitrationRL3.m
    % - experience_sbj_events(1): subject¡¯s data (1) or generated data (0) (phase 1)
    % - experience_sbj_events(2): subject¡¯s data (1) or generated data (0) (phase 2)

    if mbmf>=3 % 3Q
        mbmf = mbmf - 1;
        mode.ori = 2; % following eval_ArbitrationRL5_tauonpsa.m convention
    else
        mode.ori = 1; % following eval_ArbitrationRL5_tauonpsa.m convention
    end
    mode.USE_FWDSARSA_ONLY = round(mbmf); % 0: arbitration, 1: use fwd only, 2: use sarsa only
    mode.tau_unc = round(10 * mod(mbmf, 1)); % 1: uncertainty modulation of tau (mbmf = x.1x)
    mode.lr_unc = mod(100 * mod(mbmf, 1), 10); % 1: uncertainty modulation of lr (mbmf = x.x1)
    if any([mode.tau_unc mode.lr_unc])
        mode.param_length = 8;
        mode.param_BoundL(7) = mode.param_BoundL(5); % inv. temp.(tau) in high uncertainty
        mode.param_BoundL(8) = mode.param_BoundL(6); % lr in high uncertainty
        mode.param_BoundU(7) = mode.param_BoundU(5); % inv. temp.(tau) in high uncertainty
        mode.param_BoundU(8) = mode.param_BoundU(6); % lr in high uncertainty
    end

    mode.USE_BWDupdate_of_FWDmodel = 1; % using BWD update.
    mode.DEBUG_Q_VALUE_CHG = 0;
    mode.simul_process_display = 0;

    mode.out = fit_or_gen; % eval_ArbitrationRL3.m output control
    % - 0: OBS
    % - 1(default): Sum_NegLogLik_val
    % - 99: data_in
    % - 2: sbj

    mode.opt_ArbModel = 1; % 1: naive model(m1_wgt) . 2: posterior model(posterior)
    mode.boundary_12 = 0.01; % beta(1) = 0.01 - boundary condition for mb->mf
    mode.boundary_21 = 0.1; % alpha(1) = 0.1 - boundary condition for mf->mb
    mode.max_iter=200; % max iteration number for parameter fitting


    %% Subject behavior data loading ('data_in')
    % Load from list_sbj.
    maxi=size(list_sbj,2);% because it is for 1 subject. maxi=size(list_sbj,2);
    data_in_tag{1,1}='The columns in data in are subject numbers. Each column in the HIST_behavior_info inside represents a session.';
    % outputval={};

    % storing
    PreBehav = {};
    PreBlck = {};
    MainBehav = {};
    MainBlck = {};
    SBJ = cell(1, maxi);

    for i = 1 : maxi
        TEMP_PRE=load([list_sbj{i} '_pre_1.mat']);
        PreBehav{i}=TEMP_PRE.HIST_behavior_info{1,1};
        PreBlck{i}=TEMP_PRE.HIST_block_condition{1,1};


        tt = dir('functions/analysis_simulation/behavior_record');
        tt = {tt.name};
        maxsess = sum(cell2mat(strfind(tt,[list_sbj{i} '_fmri_']))) - 1;

        try
            for ii = 1 : maxsess
                TEMP_MAIN=load([list_sbj{i} '_fmri_' num2str(ii) '.mat']);
                MainBehav{i,ii}=TEMP_MAIN.HIST_behavior_info0;
                MainBlck{i,ii}=TEMP_MAIN.HIST_block_condition{1,ii};
            end
        catch
            warning('MAX session error');
        end

    end
    % save([pwd '/behavior_record/F_DATA.mat'], 'PreBehav', 'PreBlck','MainBehav', 'MainBlck');
    disp('DATA STORING DONE!');

    SBJtot = cell(nSimul,1);
    for d = 1 : nSimul
        %% main
        % Random parameter initialization
        param_init= zeros(1, mode.param_length);
        rng(d);
        for i = 1 : mode.param_length            
            param_init(i) = rand  * (mode.param_BoundU(i) - mode.param_BoundL(i))  + mode.param_BoundL(i);
        end

        SBJ = ArbBat_Ori(maxi,mode,PreBehav,PreBlck,MainBehav,MainBlck,list_sbj,param_init);
        SBJtot{d,1} = SBJ;
    end

    try
        save([pwd '/' save_dir '/' save_name '_' list_sbj{1} '.mat'], 'SBJtot'); 
    catch 
        mkdir(save_dir); 
        save([pwd '/' save_dir '/' save_name '_' list_sbj{1} '.mat'], 'SBJtot');
    end

end

