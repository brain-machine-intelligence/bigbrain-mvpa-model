
function Variable = arbMBMF_load_var(Exp, var_name, id, tria_cond, varargin)

switch Exp
    case 'Lee2014'
        sbjID = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24};
        % Subject Names
        Sbj_names = {'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
            'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};  % Lee2014
        try
            base_path = '/home/ydsung/A_Research';
            SBJload = load([base_path, ...
                '/Dim_control_PFC_metaRL/M3_2014subs_ori/SBJ_structure_each_exp_BESTcollectionJan29']);
            % Session marks
            % sbj_sess_mark : session change mark for MRI slices
%             sess_mark_struct = load([base_path '/Dim_control_PFC_metaRL/session_marks/sess_mark_Lee2014.mat']);
%             sess_mark_Lee2014 = sess_mark_struct.sess_mark_Lee2014;
        catch
            % Loading sbj data
%             base_path = '\\143.248.30.94\bmlsamba\ydsung/A_Research';
%             SBJload = load([base_path, ...
%                 '/Dim_control_PFC_metaRL/M3_2014subs_ori/SBJ_structure_each_exp_BESTcollectionJan29']);
            SBJload = load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace\SBJ_structure_each_exp_BESTcollectionJan29');
            % Session marks
            % sbj_sess_mark : session change mark for MRI slices
%             sess_mark_struct = load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace\sess_mark_Lee2014');
%             sess_mark_Lee2014 = sess_mark_struct.sess_mark_Lee2014;
        end
        SBJcell = struct2cell(SBJload);
        SBJs = SBJcell{:};
        % Session : 1~2, 1~3, 1~4 or 1~5
        N_sessions = 5*ones(1,length(Sbj_names)); N_sessions(11)=4; N_sessions(16)=3; N_sessions(20)=2; N_sessions(21)=4;     % Lee2014
        
    case 'Heo2018'
        
%         sbjID = 1:63; sbjID = num2cell(sbjID, 1);
%         SBJ3_load = load('C:\Users\User\Desktop\RL_depression\PLOS\revision\SBJ_structure_1027_for_revision.mat');
        SBJ3_load = load('SBJ_structure_1027_for_revision.mat');
        SBJs = SBJ3_load.SBJ3;
        N_sessions = 4*ones(1,63);      
        
        sbjID = {1, 3, 11, 12, 18, 19, 20, 24, 25, 27, 32, 33, 34, 35, 37, 38, 40, 41, 42, 44, 46, 48, 50, 51, 52, 53, 54, 55};
%         Sbj_names = {'subject001_ksy',[],'subject03_keb','subject004_ksj','subject5_lwc','subject6_sjh','subject7_kch','subject8_kjs', ...
%             'subject9_lks','subject10_ssh','subject11_hjy','subject12_khs','subject13_pjb','subject14_syh','subject15_kik','subject16_jja', ...
%             'subject16_lsh','subj18_kjh','subj19_kny','subj20_jjh','subject21_lsl','subject22_jsh','subject23_yhj','subj24_kjy','subj25_kjs', ...
%             'subject26_cjb','subj27_kkm','subject28_ljh','subject29_cyj','subject30','subj31_ssw','subj32_khj','subj33_ohy','subj34_yhw','subj35_ajs', ...
%             [],'subj37_shi','subj38_lsh',[],'subj40_kkr','subj41_jhj', 'sbj42_ljk',[],'sbj44_ksy', 'sbj45_jej','sbj46_ses','sbj47_hej','sbj48_kdy', ...
%             'subject49_ljy','sbj50_kdh','sbj51_akr', 'sbj52_jkw','sbj53_nyk','sbj54_jyh', 'sbj55_ker'};   % Heo2018
%         try
%             base_path = '/home/ydsung/A_Research';
%             SBJload = load([base_path, ...
%                 '/Dim_control_PFC_metaRL/SBJ_structure_2014_yd']);
%         catch
%             % Loading sbj data
%             base_path = '\\143.248.30.94\bmlsamba\ydsung/A_Research';
%             SBJload = load([base_path, ...
%                 '/Dim_control_PFC_metaRL/SBJ_structure_2014_yd']);
%         end
%         SBJcell = struct2cell(SBJload);
%         SBJs = SBJcell{:};
%         N_sessions = 4*ones(1,55);      % Heo2018
%         load([base_path '/Dim_control_PFC_metaRL/session_marks/mark_depression.mat']); sbj_sess_mark = mark_depression(sbj_id,:);   % Heo2018

    case 'Kim2019'
        
        sbjID = 1:24;  
        sbjID([9, 12, 16]) = [];
        sbjID = num2cell(sbjID);
        
        try
%             sbj_struct = ... 
%                 load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace/SBJ_structure_Kim2019.mat');
%             sbj_struct = ... 
%                 load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace/SBJ_structure_each_exp_113tauonpsa_10_vMF_Coin_ydsung21.mat');
            sbj_struct = ... 
                load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace/SBJ_structure_each_exp_113tauonpsa_ydrep_myObj_10_vMF_Coin21.mat');
            SBJs = sbj_struct.SBJ;
        catch
            sbj_struct = ... 
                load('/home/bmlshare/complexity/modelRLsource/result_simul/SBJ_structure.mat');
            SBJs = sbj_struct.SBJ;
        end
        
        N_sessions = [5,5,4,6,5,5,6,6,-1,5,5,-1,5,5,5,-1,6,4,5,6,5,5,5,5];
%         N_sessions = nan(length(SBJs), 1);
%         for si = 1:length(SBJs)
%             N_sessions(si) = length(SBJs{si}.HIST_behavior_info);
%         end
        
    case 'Kim2019 24'
        
        sbjID = 1:24;
        %         sbjID([9, 12, 16]) = [];
        sbjID = num2cell(sbjID);
        
        try
            %             sbj_struct = ...
            %                 load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace/SBJ_structure_Kim2019.mat');
            %             sbj_struct = ...
            %                 load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace/SBJ_structure_each_exp_113tauonpsa_10_vMF_Coin_ydsung21.mat');
%             sbj_struct = ...
%                 load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace/SBJ_structure_each_exp_113tauonpsa_10_vMF_Coin_ydsung.mat');
            sbj_struct = ... 
                load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace/SBJ_structure_each_exp_113tauonpsa_ydrep_myObj_10_vMF_Coin.mat');
            SBJs = sbj_struct.SBJ;
        catch
            sbj_struct = ...
                load('/home/bmlshare/complexity/modelRLsource/result_simul/SBJ_structure.mat');
            SBJs = sbj_struct.SBJ;
        end
        
        %         N_sessions = [5,5,4,6,5,5,6,6,-1,5,5,-1,5,5,5,-1,6,4,5,6,5,5,5,5];
        N_sessions = nan(length(SBJs), 1);
        for si = 1:length(SBJs)
            N_sessions(si) = length(SBJs{si}.HIST_behavior_info);
        end
         
end

if id == 0
    
    mean_target_id = nan * ones(1,length(sbjID));
    for idi = 1:length(sbjID)
        
        temp_tria_cond = [];
        
        
        if ischar(tria_cond)
            
            temp_tria_cond = tria_cond;
%             switch tria_cond
%                 case 'G' % for Lee 2014
%                     GT = arbMBMF_load_var(Exp, 'GoalCond', idi, []); % 1:G, 2:H
%                     temp_tria_cond = ismember(GT, 1);
% %                     fprintf('2')
%                 case 'G uL' % for Lee 2014
%                     context = arbMBMF_load_var(Exp, 'blkCond', idi, []);
%                     temp_tria_cond = (context == 1);
%                 case 'G uH' % for Lee 2014
%                     context = arbMBMF_load_var(Exp, 'blkCond', idi, []);
%                     temp_tria_cond = (context == 2);
%                 case 'uL'
%                     UC = arbMBMF_load_var(Exp, 'UncCond', idi, []); % 1:low, 2:high
%                     temp_tria_cond = (UC~=1);
%                 case 'uH'
%                     UC = arbMBMF_load_var(Exp, 'UncCond', idi, []); % 1:low, 2:high
%                     temp_tria_cond = (UC~=2);
%                 case 'binPMB0'
%                     PMB = arbMBMF_load_var(Exp, 'PMB', idi, []);
%                     binPMB = regr2class(PMB, 2);
%                     temp_tria_cond = (binPMB==0);
% %                     fprintf('test ')
%                 case 'binPMB1'
%                     PMB = arbMBMF_load_var(Exp, 'PMB', idi, []);
%                     binPMB = regr2class(PMB, 2);
%                     temp_tria_cond = (binPMB==1);
%                 case 'binPMB1 uL'
%                     PMB = arbMBMF_load_var(Exp, 'PMB', idi, []); binPMB = regr2class(PMB, 2);
%                     UC = arbMBMF_load_var(Exp, 'UncCond', idi, []); % 1:low, 2:high
%                     temp_tria_cond = (binPMB==1 & UC==1); 
%                 case 'binPMB1 uH'
%                     PMB = arbMBMF_load_var(Exp, 'PMB', idi, []); binPMB = regr2class(PMB, 2);
%                     UC = arbMBMF_load_var(Exp, 'UncCond', idi, []); % 1:low, 2:high
%                     temp_tria_cond = (binPMB==1 & UC==2); 
%                 case 'G uLate'
%                     GT = arbMBMF_load_var(Exp, 'GoalCond', idi, []); % 1:G, 2:H
%                     UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', idi, []); % 0:early, 1:late
%                     temp_tria_cond = (GT==1 & UncPhaseEL==1);
% %                     disp('test')
%                 case 'G uEarly'
%                     GT = arbMBMF_load_var(Exp, 'GoalCond', idi, []); % 1:G, 2:H
%                     UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', idi, []); % 0:early, 1:late
%                     temp_tria_cond = (GT==1 & UncPhaseEL==0);
%                 case 'binPMB1 uLate'
%                     PMB = arbMBMF_load_var(Exp, 'PMB', idi, []); binPMB = regr2class(PMB, 2);
%                     UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', idi, []); % 0:early, 1:late
%                     temp_tria_cond = (binPMB==1 & UncPhaseEL==1);
%                 case 'binPMB1 uEarly'
%                     PMB = arbMBMF_load_var(Exp, 'PMB', idi, []); binPMB = regr2class(PMB, 2);
%                     UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', idi, []); % 0:early, 1:late
%                     temp_tria_cond = (binPMB==1 & UncPhaseEL==0);
%             end
            
        else
            
            % goalCond
            if tria_cond == 1
                % 200410: temp_tria_cond - for G trials
                temp_tria_cond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
                temp_tria_cond = ismember(temp_tria_cond, [1, 2]);
%                 fprintf('t2 ')
            end
            if tria_cond == 2
                % 200417: temp_tria_cond - for H trials
                temp_tria_cond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
                temp_tria_cond = ismember(temp_tria_cond, [3, 4]);
            end
            % uncCond
            if tria_cond == 3
                % 200417: temp_tria_cond - for lUnc trials
                temp_tria_cond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
                temp_tria_cond = ismember(temp_tria_cond, [1, 4]);
            end
            if tria_cond == 4
                % 200417: temp_tria_cond - for hUnc trials
                temp_tria_cond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
                temp_tria_cond = ismember(temp_tria_cond, [2, 3]);
            end
            % blkCond
            if tria_cond == 11
                blkCond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
                temp_tria_cond = ismember(blkCond, 1);
            end
            if tria_cond == 22
                blkCond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
                temp_tria_cond = ismember(blkCond, 2);
            end
            if tria_cond == 33
                blkCond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
                temp_tria_cond = ismember(blkCond, 3);
            end
            if tria_cond == 44
                blkCond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
                temp_tria_cond = ismember(blkCond, 4);
            end
        end
        
        target = arbMBMF_load_var(Exp, var_name, idi, temp_tria_cond);
        if isvector(target)
            mean_target_id(idi) = nanmean(target);
        else
            mean_target_id(idi) = nanmean(nanmean(target));
        end
                
    end
    Variable = mean_target_id;
    
else
    
    switch Exp
        case 'Lee2014'
            SBJ = SBJs{id};
            N_session = N_sessions(sbjID{id});
%             sbj_sess_mark = sess_mark_Lee2014(sbjID{id},:);
            
            % for relMB simplicity
            % relMB
            relMB = SBJ.regressor{1,7}.value(8,:)./...
                sum(SBJ.regressor{1,7}.value(7:9,:), 1);
            relMB(1:3:end) = nan;
%             skips = ceil(0.15 * length(relMB));
            skips = 1;
            relMB(1:skips) = nan;
            relMB = (relMB-nanmin(relMB(skips:end)))/ ...
                (nanmax(relMB(skips:end))-nanmin(relMB(skips:end)));
            relMB = relMB - nanmean(relMB);
            % relMF
            relMF = SBJ.regressor{1,8}.value(8,:);
            relMF(1:3:end) = nan;
%             skips = ceil(0.15 * length(relMF));
            skips = 1;
            relMF(1:skips) = nan;
            relMF = (relMF-nanmin(relMF(skips:end)))/ ...
                (nanmax(relMF(skips:end))-nanmin(relMF(skips:end)));
            relMF = relMF - nanmean(relMF);
            
        case 'Heo2018'
            SBJ = SBJs{sbjID{id}};
            N_session = N_sessions(sbjID{id});
            
        case {'Kim2019', 'Kim2019 24'}
            SBJ = SBJs{id};
            N_session = N_sessions(sbjID{id});
            
            % relMB
            relMB = SBJ.regressor{1,7}.value(8,:)./...
                sum(SBJ.regressor{1,7}.value(7:9,:), 1);
            relMB(1:3:end) = nan;
%             skips = ceil(0.15 * length(relMB));
            skips = 1;
            relMB(1:skips) = nan;
            relMB = (relMB-nanmin(relMB(skips:end)))/ ...
                (nanmax(relMB(skips:end))-nanmin(relMB(skips:end)));
            relMB = relMB - nanmean(relMB);
            % relMF
            relMF = SBJ.regressor{1,8}.value(8,:);
            relMF(1:3:end) = nan;
%             skips = ceil(0.15 * length(relMF));
            skips = 1;
            relMF(1:skips) = nan;
            relMF = (relMF-nanmin(relMF(skips:end)))/ ...
                (nanmax(relMF(skips:end))-nanmin(relMF(skips:end)));
            relMF = relMF - nanmean(relMF);
            
    end
    
    switch var_name
        
        case 0
            reward = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                reward = [reward, sess_behav_info(:,16)'];
            end
            Variable = ones(size(reward));

% -------------------------------------------------------------------------
% exp. condition: session
            
        case 'Null'
            Variable = [];
            for sess = 1:N_session
                temp_sess = sess * ones(1, size(SBJ.HIST_behavior_info{sess}, 1));
                Variable = [Variable, temp_sess];
            end
            Variable(:) = 1;
        
        case 'Session' % session
            Variable = [];
            for sess = 1:N_session
                temp_sess = sess * ones(1, size(SBJ.HIST_behavior_info{sess}, 1));
                Variable = [Variable, temp_sess];
            end
            
        case 'Half'
            temp_Variable = [];
            for sess = 1:N_session
                temp_sess = sess * ones(1, size(SBJ.HIST_behavior_info{sess}, 1));
                temp_Variable = [temp_Variable, temp_sess];
            end
            Variable = temp_Variable;
            Variable(1:floor(length(Variable)/2)) = 1;
            Variable(ceil(length(Variable)/2):end) = 2;
            
        case 'EarlyLate'
            temp_Variable = [];
            for sess = 1:N_session
                temp_sess = sess * ones(1, size(SBJ.HIST_behavior_info{sess}, 1));
                temp_Variable = [temp_Variable, temp_sess];
            end
            Variable = nan(size(temp_Variable));
            Variable(1:floor(length(Variable)/3)) = 1;
            Variable(ceil(2*length(Variable)/3):end) = 2;
            
% -------------------------------------------------------------------------
% exp. condition: stimulus (state)
            
        case 'S2' % 2 3 4 5
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,5);
                Variable = [Variable, sess_lab'];
            end
            
        case 'S3' % 6 7 8 9
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,6);
                Variable = [Variable, sess_lab'];
            end
            % sanity check
            if ~all(ismember(Variable, [6 7 8 9]))
                error('(arbMBMF_load_var) S3 error');
            end
            
        case 'Coin' % 6 7 8 9
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,6);
                Variable = [Variable, sess_lab'];
            end
            switch Exp
                case 'Kim2019'
                    temp = Variable;
                    Variable(ismember(temp, [5 6])) = 1;
                    Variable(ismember(temp, [7 9])) = 2;
                    Variable(ismember(temp, [8 11])) = 3;
                    Variable(ismember(temp, [4 10])) = 4;
            end
            
        case 'S3G'
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,6);
                Variable = [Variable, sess_lab'];
            end
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
            end
            Variable(blk_con==3) = nan;
            Variable(blk_con==4) = nan;
            
        case 10 % S2 in G
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,5);
                Variable = [Variable, sess_lab'];
            end
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
            end
            Variable(blk_con==3) = nan;
            Variable(blk_con==4) = nan;
                    
        case 12 % S2 in H
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,5);
                Variable = [Variable, sess_lab'];
            end
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
            end
            Variable(blk_con==1) = nan;
            Variable(blk_con==2) = nan;
            
        case 13 % S3 in H
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,6);
                Variable = [Variable, sess_lab'];
            end
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
            end
            Variable(blk_con==1) = nan;
            Variable(blk_con==2) = nan;
            
            
% -------------------------------------------------------------------------
% exp. condition: block condition (context)
            
        case 'blkCond' % block condition 1:G, 2:G', 3:H, 4:H'
            blk_con = [];
            for sess=1:length(SBJ.HIST_block_condition)
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            Variable = blk_con;
            
        case 'GoalCond' % block condition [1, 2]:G, [3, 4]:H
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            Variable = (blk_con > 2.5) + 1; % 1:G, 2:H
            
        case 'UncCond'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    Variable = ismember(blk_con, [2 3]) + 1; % 1:Low, 2:High
                case 'Kim2019'
                    Variable = ismember(blk_con, [3 4]) + 1; % 1:Low, 2:High
            end
            
        case 'Session x UC'
            switch Exp
                case 'Lee2014'
                    BHVcell = SBJ.HIST_behavior_info;
                    
                    % Session
                    sess_lens = cellfun(@(x) size(x, 1), BHVcell);
                    Session = cell(1, length(sess_lens));
                    for s = 1:length(sess_lens)
                        Session{s} = s * ones(1, sess_lens(s));
                    end
                    Session = cell2mat(Session); el_sess = unq_elms(Session);
                    
                    % UC
                    BHVs = cell2mat(BHVcell');
                    blkCond = BHVs(:, 3)'; % 1 2 3 4
                    UncCond = ismember(blkCond, [2 3]) + 1; % 1:Low, 2:High
                    el_UC = unq_elms(UncCond);
                    
                    comb_lab = [Session; UncCond];
                    new_el = combvec(el_sess, el_UC);
                    Variable = nan(size(UncCond));
                    for k = 1:size(new_el, 2)
                        Variable(all(comb_lab == new_el(:, k), 1)) = k;
                    end
                    
            end
            
        case 'UncCond at state 4'
            
            switch Exp
                case 'Lee2014'
                    BHVcell = SBJ.HIST_behavior_info;
                    BHVs = cell2mat(BHVcell');
                    blkCond = BHVs(:, 3)'; % 1 2 3 4
                    UncCond = ismember(blkCond, [2 3]) + 1; % 1:Low, 2:High
                    S2 = BHVs(:, 5)'; % 2 3 4 5
                    
                    Variable = UncCond;
                    Variable(~ismember(S2, 4)) = nan;
            end
                        
        case 'CmplxCond'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            Variable = ismember(blk_con, [2 4]) + 1; % 1:low, 2:high
            
        case 'Interaction'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            Variable = nan(size(blk_con));
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    Variable(ismember(blk_con, [1 3])) = 1;
                    Variable(ismember(blk_con, [2 4])) = 2;
                case 'Kim2019'
                    Variable(ismember(blk_con, [1 4])) = 1;
                    Variable(ismember(blk_con, [2 3])) = 2;
            end
            
        % block/goaltype/unc/cmplx/goal switch timing
        case 'blkSwitch'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            Variable = [1 diff(blk_con)~=0];
            
        case 'GoaltypeSwitch'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = [1 diff(GT)~=0];
            
        case 'UncSwitch'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                case 'Kim2019'
                    UC = ismember(blk_con, [1 2]);  % 1: Low, 0: High
            end
            Variable = [1 diff(UC)~=0];
            
        case 'CmplxSwitch'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            CX = ismember(blk_con, [1 3]);  % 1: Low, 0: High
            Variable = [1 diff(CX)~=0];
            
        case 'GoalSwitch'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                    Variable = [];
                    for sess = 1:N_session
                        sess_lab = SBJ.HIST_behavior_info{sess}(:,18);
                        Variable = [Variable, sess_lab'];
                    end
                    
                case 'Kim2019'
                    % 1: 1st outcome has max value, 2: 2nd, 3: 3rd
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal, [], 2);
                        bin_max = (goalVal == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        Variable = [Variable, sess_lab'];
                    end
                    Variable(isnan(Variable)) = -1;
            end
            
            Variable = [1 diff(Variable)~=0];
            
        case 'Goal(6,7,8)Switch'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                    Variable = [];
                    for sess = 1:N_session
                        sess_lab = SBJ.HIST_behavior_info{sess}(:,18);
                        Variable = [Variable, sess_lab'];
                    end
                    blk_con = [];
                    for sess=1:N_session
                        blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
                    end
                    
                case 'Kim2019'
                    % 1: 1st outcome has max value, 2: 2nd, 3: 3rd
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal, [], 2);
                        bin_max = (goalVal == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        Variable = [Variable, sess_lab'];
                    end
                    Variable(isnan(Variable)) = -1;
                    blk_con = [];
                    for sess=1:N_session
                        blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
                    end
                    
            end
            
            Variable = [1 diff(Variable)~=0];
            Variable(~ismember(blk_con, [1 2 ])) = nan;
            
        case 'Goal(6,7,8)xUCSwitch'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                    Variable = [];
                    for sess = 1:N_session
                        sess_lab = SBJ.HIST_behavior_info{sess}(:,18);
                        Variable = [Variable, sess_lab'];
                    end
                    blk_con = [];
                    for sess=1:N_session
                        blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
                    end
                    
                case 'Kim2019'
                    % 1: 1st outcome has max value, 2: 2nd, 3: 3rd
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal, [], 2);
                        bin_max = (goalVal == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        Variable = [Variable, sess_lab'];
                    end
                    Variable(isnan(Variable)) = -1;
                    blk_con = [];
                    for sess=1:N_session
                        blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
                    end
                    
            end
            
            Variable = [1 diff(Variable)~=0];
            blk_con_sw = [1 diff(blk_con)~=0];
            Variable = 1 * (Variable | blk_con_sw);
            Variable(~ismember(blk_con, [1 2])) = nan;
            
        % phase (e.g. early/late) after block change
        case 'blkPhase'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
%             Variable = class_series2phase_series(blk_con);
            Variable = class_series2phase_series(blk_con, 'PhaseNumber', 3);
            
        case 'GoalPhase'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                    Variable = [];
                    for sess = 1:N_session
                        sess_lab = SBJ.HIST_behavior_info{sess}(:,18);
                        Variable = [Variable, sess_lab'];
                    end
                    
                case 'Kim2019'
                    % 1: 1st outcome has max value, 2: 2nd, 3: 3rd
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal, [], 2);
                        bin_max = (goalVal == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        Variable = [Variable, sess_lab'];
                    end
                    Variable(isnan(Variable)) = -1;
            end
            
            Variable = class_series2phase_series(Variable, 'PhaseNumber', 3);
            
        case 'GoalPhase2'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                    Variable = [];
                    for sess = 1:N_session
                        sess_lab = SBJ.HIST_behavior_info{sess}(:,18);
                        Variable = [Variable, sess_lab'];
                    end
                    
                case 'Kim2019'
                    % 1: 1st outcome has max value, 2: 2nd, 3: 3rd
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal, [], 2);
                        bin_max = (goalVal == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        Variable = [Variable, sess_lab'];
                    end
                    Variable(isnan(Variable)) = -1;
            end
            
            Variable = class_series2phase_series(Variable, 'PhaseNumber', 2);
        
        case 'GoaltypePhase'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = class_series2phase_series(GT, 'PhaseNumber', 3);
            
        case 'GoaltypePhase2'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = class_series2phase_series(GT, 'PhaseNumber', 2);
            
        case 'SpecGoalPhase2'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = class_series2phase_series(GT, 'PhaseNumber', 2);
            Variable(GT==0) = nan;
            
        case 'SpecGoalPhase3'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = class_series2phase_series(GT, 'PhaseNumber', 3);
            Variable(GT==0) = nan;
            
        case 'FlexGoalPhase2'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = class_series2phase_series(GT, 'PhaseNumber', 2);
            Variable(GT==1) = nan;
            
        case 'UncPhase'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                case 'Kim2019'
                    UC = ismember(blk_con, [1 2]);  % 1: Low, 0: High
            end
            Variable = class_series2phase_series(UC, 'PhaseNumber', 3);
            
        case 'UncPhase2'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                case 'Kim2019'
                    UC = ismember(blk_con, [1 2]);  % 1: Low, 0: High
            end
            Variable = class_series2phase_series(UC, 'PhaseNumber', 2);
        
        case 'UncPhaseEL' % Early (0) & Late (1) phase only
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                case 'Kim2019'
                    UC = ismember(blk_con, [1 2]);  % 1: Low, 0: High
            end
            Variable = class_series2phase_series(UC, 'PhaseNumber', 3);
            Variable(Variable==2) = nan;
            Variable(Variable==1) = 0; Variable(Variable==3) = 1;
            
        case 'CmplxPhase'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            CX = ismember(blk_con, [1 3]);  % 1: Low, 0: High
            Variable = class_series2phase_series(CX, 'PhaseNumber', 3);
            
        case 'CmplxPhase2'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            CX = ismember(blk_con, [1 3]);  % 1: Low, 0: High
            Variable = class_series2phase_series(CX, 'PhaseNumber', 2);

% -------------------------------------------------------------------------
% exp. condition: goal

        case 'Goal' % max value goal = argmax_{goal} (value)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                    Variable = [];
                    for sess = 1:N_session
                        sess_lab = SBJ.HIST_behavior_info{sess}(:,18);
                        Variable = [Variable, sess_lab'];
                    end
                case 'Kim2019'
                    % 1: 1st outcome has max value, 2: 2nd, 3: 3rd
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal, [], 2);
                        bin_max = (goalVal == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        Variable = [Variable, sess_lab'];
                    end
            end
            
        case 'GoalValEntropy' % entropy of normalized goal currency
            switch Exp
                case 'Kim2019'
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
                        normVal = goalVal ./ sum(goalVal, 2);
                        ent_sess = nansum(- normVal .* log2(normVal), 2);
                        Variable = [Variable, ent_sess'];
                    end
            end
            
        case 'binGoalValEntropy' % entropy of normalized goal currency
            switch Exp
                case 'Kim2019'
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
                        normVal = goalVal ./ sum(goalVal, 2);
                        ent_sess = nansum(- normVal .* log2(normVal), 2);
                        Variable = [Variable, ent_sess'];
                    end
                    Variable = regr2class(Variable, 2);
            end
            
        case 'AchGoalS1RandA'  % = argmax_{goal} (achievability x value) in stage 1, random action
            switch Exp
                case 'Kim2019'
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20); % goal values
                        goalAch = nan(size(goalVal, 1), 3); % goal achievabilities
                        context = SBJ.HIST_block_condition{sess}(2,:); % block conditions
                        S2 = SBJ.HIST_behavior_info{sess}(:,5);
                        for si = 1:length(S2)
                            d = goal_distr_complex(context(si), 1, ...
                                'TransitionStep', 2);
                            goalAch(si, :) = d(1:3);
                        end
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal.*goalAch, [], 2);
                        bin_max = (goalVal.*goalAch == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        Variable = [Variable, sess_lab'];
                    end
            end
            
        case 'AchGoalS1OptA'  % = argmax_{goal} (achievability x value) in stage 1, optimal action
            switch Exp
                case 'Kim2019'
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20); % goal values
                        goalAch = nan(size(goalVal, 1), 3); % goal achievabilities
                        context = SBJ.HIST_block_condition{sess}(2,:); % block conditions
                        S2 = SBJ.HIST_behavior_info{sess}(:,5);
                        for si = 1:length(S2)
                            d = goal_distr_complex(context(si), 1, ...
                                'Policy', 'optimal', ...
                                'TransitionStep', 2, ...
                                'GoalValues', goalVal(si, :));
                            goalAch(si, :) = d(1:3);
                        end
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal.*goalAch, [], 2);
                        bin_max = (goalVal.*goalAch == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        Variable = [Variable, sess_lab'];
                    end
            end
            
        case 'AchGoalS2RandA'  % = argmax_{goal} (achievability x value) in stage 2, random action
            switch Exp
                case 'Kim2019'
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20); % goal values
                        goalAch = nan(size(goalVal, 1), 3); % goal achievabilities
                        context = SBJ.HIST_block_condition{sess}(2,:); % block conditions
                        S2 = SBJ.HIST_behavior_info{sess}(:,5);
                        for si = 1:length(S2)
                            d = goal_distr_complex(context(si), S2(si));
                            goalAch(si, :) = d(1:3);
                        end
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal.*goalAch, [], 2);
                        bin_max = (goalVal.*goalAch == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        Variable = [Variable, sess_lab'];
                    end
            end
            
        case 'AchGoalS2OptA'  % = argmax_{goal} (achievability x value) in stage 2, optimal action
            switch Exp
                case 'Kim2019'
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20); % goal values
                        goalAch = nan(size(goalVal, 1), 3); % goal achievabilities
                        context = SBJ.HIST_block_condition{sess}(2,:); % block conditions
                        S2 = SBJ.HIST_behavior_info{sess}(:,5);
                        for si = 1:length(S2)
                            d = goal_distr_complex(context(si), S2(si), ...
                                'Policy', 'optimal', ...
                                'GoalValues', goalVal(si, :));
                            goalAch(si, :) = d(1:3);
                        end
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal.*goalAch, [], 2);
                        bin_max = (goalVal.*goalAch == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        Variable = [Variable, sess_lab'];
                    end
            end
            
        case 'Goal(6,7,8)'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                    Variable = [];
                    for sess = 1:N_session
                        sess_lab = SBJ.HIST_behavior_info{sess}(:,18);
                        Variable = [Variable, sess_lab'];
                    end
            end
            Variable(Variable == -1) = nan;
            
%         case 'GoalSwitch'
%             switch Exp
%                 case {'Lee2014', 'Heo2018'}
%                     % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
%                     G = [];
%                     for sess = 1:N_session
%                         sess_lab = SBJ.HIST_behavior_info{sess}(:,18);
%                         G = [G, sess_lab'];
%                     end
%                 case 'Kim2019'
%                     % 1: 1st outcome has max value, 2: 2nd, 3: 3rd
%                     G = [];
%                     for sess = 1:N_session
%                         goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
%                         sess_lab = nan(size(goalVal, 1), 1);
%                         maxgoal = max(goalVal, [], 2);
%                         bin_max = (goalVal == maxgoal);
%                         sess_lab(bin_max(:,1)) = 1;
%                         sess_lab(bin_max(:,2)) = 2;
%                         sess_lab(bin_max(:,3)) = 3;
%                         sess_lab(sum(bin_max,2)>1) = nan;
%                         G = [G, sess_lab'];
%                     end
%             end
%             Variable = [1 diff(G)~=0];
            
        case 'GoalValues' % mainly for Kim 2019
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % not implemented
                case 'Kim2019'
%                     % 1: 1st outcome has max value, 2: 2nd, 3: 3rd
%                     Variable = [];
%                     for sess = 1:N_session
%                         goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
%                         sess_lab = nan(size(goalVal, 1), 1);
%                         maxgoal = max(goalVal, [], 2);
%                         bin_max = (goalVal == maxgoal);
%                         sess_lab(bin_max(:,1)) = 1;
%                         sess_lab(bin_max(:,2)) = 2;
%                         sess_lab(bin_max(:,3)) = 3;
%                         sess_lab(sum(bin_max,2)>1) = nan;
%                         Variable = [Variable, sess_lab'];
%                     end
                    
                    % not implemented
            end
            
% -------------------------------------------------------------------------
% exp. condition: episode (S1-A1-S2-A2-S3)
        case 'episode'
            HISTbhv = SBJ.HIST_behavior_info;
            episodes = [];
            for ss = 1:N_session
                episodes = [episodes; num2str(HISTbhv{ss}(:,4:8))];
            end
            
            Variable = nan(1, size(episodes,1));
            epi_set = unique(episodes, 'rows');
            for e = 1:size(epi_set,1)
                Variable(all(episodes==epi_set(e,:),2)) = e;
            end
            
% -------------------------------------------------------------------------
% behavior: action
            
        case 'A1' % A1 1:left, 2:right?
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,7);
                Variable = [Variable, sess_lab'];
            end
            
        case 'A2' % A2 1:left, 2:right?
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,8);
                Variable = [Variable, sess_lab'];
            end
            
        case 'A2 at state 4'
            
            switch Exp
                case 'Lee2014'
                    BHVcell = SBJ.HIST_behavior_info;
                    BHVs = cell2mat(BHVcell');
                    S2 = BHVs(:, 5)'; % 2 3 4 5
                    A2 = BHVs(:, 8)'; % 1 2
                    
                    Variable = A2;
                    Variable(~ismember(S2, 4)) = nan;
            end
            
% -------------------------------------------------------------------------
% behavior: performance
            
        case 'Hit'  % Hit
            reward = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                reward = [reward, sess_behav_info(:,16)'];
            end
            Variable = zeros(size(reward));
            Variable(reward > 0) = 1;
            
        case 'R'  % Reward
            reward = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                reward = [reward, sess_behav_info(:,16)'];
            end
            Variable = reward;
%             Variable = reward/max(reward);
            
% -------------------------------------------------------------------------
% behavior: choice optimality
            
        case 'ChoOpt_old'  % Choice optimality (Lee2014)
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
                
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                sess_blkcon = sess_behav_info(:,3)';
                S2 = [S2, sess_behav_info(:,5)'];
                A2 = [A2, sess_behav_info(:,8)'];
                goal_state = [goal_state, sess_behav_info(:,18)'];
                
                if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0         % Sanity check
                    error('block condition information error'); end
            end
            
            % Finding choice optimality
            is_opt_choice = zeros(size(goal_state));
            
            % specific/low_unc/40goal
            idx_GL40_opt = [find(blk_con==1 & goal_state==6 & S2 == 4 & A2 == 2), ...
                find(blk_con==1 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/low_unc/20goal
            idx_GL20_opt = [find(blk_con==1 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==1 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/low_unc/10goal
            idx_GL10_opt = [find(blk_con==1 & goal_state==8 & S2 == 2 & A2 == 2), ...
                find(blk_con==1 & goal_state==8 & S2 == 3 & A2 == 1)];
            % specific/high_unc/40goal
            idx_GH40_opt = [find(blk_con==2 & goal_state==6 & S2 == 4 & A2 > 0), ...
                find(blk_con==2 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/high_unc/20goal
            idx_GH20_opt = [find(blk_con==2 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==2 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/high_unc/10goal
            idx_GH10_opt = [find(blk_con==2 & goal_state==8 & S2 == 2 & A2 > 0), ...
                find(blk_con==2 & goal_state==8 & S2 == 3 & A2 == 1)];
            % flexible/high_unc/uni_goal
            idx_FH_opt = [find(blk_con==3 & goal_state==-1 & S2 == 4 & A2 == 1), ...
                find(blk_con==3 & goal_state==-1 & S2 == 5 & A2 == 2)];
            % flexible/low_unc/uni_goal
            idx_FL_opt = [find(blk_con==4 & goal_state==-1 & S2 == 4 & A2 == 2), ...
                find(blk_con==4 & goal_state==-1 & S2 == 5 & A2 == 1)];
            
            is_opt_choice(idx_GL40_opt) = 1; is_opt_choice(idx_GH40_opt) = 1;
            is_opt_choice(idx_GL20_opt) = 1; is_opt_choice(idx_GH20_opt) = 1;
            is_opt_choice(idx_GL10_opt) = 1; is_opt_choice(idx_GH10_opt) = 1;
            is_opt_choice(idx_FH_opt) = 1; is_opt_choice(idx_FL_opt) = 1;
            
            Variable = is_opt_choice;
            
        case 'ChoOpt' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
                    
                case 'Kim2019'
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal, blk, 'likelihood', varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal, blk, 'likelihood', varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
            Variable = mean([opt_likelihood1; opt_likelihood2], 1); % [2 x N_trial]
            
        case 'ChoOptG' % **** different definition between Lee 2014 vs Kim 2019
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    blk_con = [];
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        blk_con = [blk_con, temp_blk_con(2,:)];
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    context = ismember(blk_con, [1 2]);
                    
                case 'Kim2019'
                    PMB = SBJ.regressor{1,9}.value(7,:); % weigtM1 value (Lee2014, Kim2019)
                    PMB = (PMB(1:3:end)+PMB(2:3:end)+PMB(3:3:end))/3;
                    binPMB = regr2class(PMB, 2);
                    context = ismember(binPMB, 1);
                    
                    blk_con = [];
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        blk_con = [blk_con, temp_blk_con(2,:)];
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal, blk, 'likelihood');
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal, blk, 'likelihood');
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
            Variable = mean([opt_likelihood1; opt_likelihood2], 1); % [2 x N_trial]
            Variable(~context) = nan;
            
        case 'ChoOptH' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    blk_con = [];
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        blk_con = [blk_con, temp_blk_con(2,:)];
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
                    
                case 'Kim2019'
                    disp('not implemented yet')
            end
            
            Variable = mean([opt_likelihood1; opt_likelihood2], 1); % [2 x N_trial]
            Variable(~ismember(blk_con, [3 4])) = nan;
            
        case 'ChoOpt1n2' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
                    
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:length(SBJ.HIST_block_condition)
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal, blk, 'likelihood', varargin{:});
                            % L1 = choice_opt_complex(blk, state1, action1, goal, varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal, blk, 'likelihood', varargin{:});
                            % L2 = choice_opt_complex(blk, state2, action2, goal, varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
        case 'ChoOpt1n2TransEst' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
                    
                case {'Kim2019', 'Kim2019 24'}
                    temp_T = SBJ.('S1_state_fwd_T'); % N_state x N_action x N_state x N_trial
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    t_cnt = 0;
                    for sess=1:length(SBJ.HIST_block_condition)
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            t_cnt = t_cnt + 1;
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = choice_opt_complex(blk, state1, action1, goal, ... 
                                'TransitionProb', permute(temp_T(:,:,:,t_cnt), [4 3 1 2]), varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = choice_opt_complex(blk, state2, action2, goal, ... 
                                'TransitionProb', permute(temp_T(:,:,:,t_cnt), [4 3 1 2]), varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
                        
        case 'ChoOptA1R1'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % not implemented yet
                    
                case 'Kim2019'
                    opt_likelihood1 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            L1 = opt_normalization_1st2ndstage_complex(state1, 1, ... % A1 = R1
                                goal, blk, 'likelihood', varargin{:}); 
                            
                            opt_L1_sess(t_sess) = L1;                            
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                                                
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood1; % [1 x N_trial]
            end
            
        case 'ChoOptA1L1'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % not implemented yet
                    
                case 'Kim2019'
                    opt_likelihood1 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            L1 = opt_normalization_1st2ndstage_complex(state1, 2, ... % A1 = L1
                                goal, blk, 'likelihood', varargin{:}); 
                            
                            opt_L1_sess(t_sess) = L1;                            
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                                                
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood1; % [1 x N_trial]
            end
            
        case 'ChoOptA1LeftBias'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % not implemented yet
                    
                case 'Kim2019'
                    opt_likelihood1 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            L1_L = opt_normalization_1st2ndstage_complex(state1, 2, ... % A1 = L1
                                goal, blk, 'likelihood', varargin{:}); 
                            L1_R = opt_normalization_1st2ndstage_complex(state1, 1, ... % A1 = R1
                                goal, blk, 'likelihood', varargin{:}); 
                            
                            opt_L1_sess(t_sess) = L1_L - L1_R;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                                                
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood1; % [1 x N_trial]
            end
            
        case 'ChoOptA2LeftBias' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % not implemented yet
                    
                case 'Kim2019'
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal, blk, 'leftbias', varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood2; % [1 x N_trial]
            end
                        
        case 'ChoOptBin1n2' % choice optimality (Kim2021)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            if L1==.5; L1=.5; elseif L1>.5; L1=1; else; L1=0; end % discretization, Kim 2021
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            if L2==.5; L2=.5; elseif L2>.5; L2=1; else; L2=0; end % discretization, Kim 2021
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
                    
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal, blk, 'binary', varargin{:}); % discretization, Kim 2021
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal, blk, 'binary', varargin{:}); % discretization, Kim 2021
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
                    
        case 'ChoOptBin1' % choice optimality (Kim2021)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = [];
                    
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            if L1==.5; L1=.5; elseif L1>.5; L1=1; else; L1=0; end % discretization, Kim 2021
                            opt_L1_sess(t_sess) = L1;
                            
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = opt_likelihood1; % [1 x N_trial]
                    
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood1 = [];
                    
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal, blk, 'binary', varargin{:}); % discretization, Kim 2021
                            
                            opt_L1_sess(t_sess) = L1;
                            
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                                                
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood1; % [1 x N_trial]
            end
            
        case 'ChoOptBin2' % choice optimality (Kim2021)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            if L2==.5; L2=.5; elseif L2>.5; L2=1; else; L2=0; end % discretization, Kim 2021
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = opt_likelihood2; % [1 x N_trial]
                    
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal, blk, 'binary', varargin{:}); % discretization, Kim 2021
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood2; % [1 x N_trial]
            end
            
        case 'ChoOptBinNew1n2' % strictly binary
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            if L1>=.5; L1=1; else; L1=0; end % discretization, Kim 2021
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            if L2>=.5; L2=1; else; L2=0; end % discretization, Kim 2021
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
                    
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal, blk, 'binary', varargin{:}); % discretization, Kim 2021
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal, blk, 'binary', varargin{:}); % discretization, Kim 2021
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
        case 'ChoOptBinNew1' % strictly binary
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            if L1>=.5; L1=1; else; L1=0; end % discretization, Kim 2021
                            opt_L1_sess(t_sess) = L1;
                            
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = opt_likelihood1; % [1 x N_trial]
                    
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood1 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal, blk, 'binary', varargin{:}); % discretization, Kim 2021
                            
                            opt_L1_sess(t_sess) = L1;
                            
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood1; % [1 x N_trial]
            end
            
        case 'ChoOptBinNew2' % strictly binary
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            if L2>=.5; L2=1; else; L2=0; end % discretization, Kim 2021
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = opt_likelihood2; % [1 x N_trial]
                    
                case {'Kim2019', 'Kim2019 24'}
                    
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal, blk, 'binary', varargin{:}); % discretization, Kim 2021
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood2; % [1 x N_trial]
            end
            
        case 'ChoOpt1n2MaxGoal' % choice optimality (Kim2019)
            switch Exp
                
                case {'Kim2019', 'Kim2019 24'}
%                     max_goal = arbMBMF_load_var(Exp, 'Goal', id, []);
%                     ach_goal = arbMBMF_load_var(Exp, 'AchGoalS2RandA', id, []);
                    
                    t_cnt = 0; % trial count
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:length(SBJ.HIST_block_condition)                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            t_cnt = t_cnt + 1;
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            max_goal = max(goal);
                            new_goal = zeros(size(goal));
%                             new_goal(max_goal(t_cnt)) = goal(max_goal(t_cnt));
                            new_goal(goal==max_goal) = goal(goal==max_goal);

                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, new_goal, blk, 'likelihood', varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, new_goal, blk, 'likelihood', varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
        case 'ChoOpt1n2MaxGoalTransEst' % choice optimality (Kim2019)
            switch Exp
                
                case {'Kim2019', 'Kim2019 24'}
                    temp_T = SBJ.('S1_state_fwd_T'); % N_state x N_action x N_state x N_trial
                    t_cnt = 0; % trial count
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:length(SBJ.HIST_block_condition)                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            t_cnt = t_cnt + 1;
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            max_goal = max(goal);
                            new_goal = zeros(size(goal));
%                             new_goal(max_goal(t_cnt)) = goal(max_goal(t_cnt));
                            new_goal(goal==max_goal) = goal(goal==max_goal);

                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = choice_opt_complex(blk, state1, action1, new_goal, ... 
                                'TransitionProb', permute(temp_T(:,:,:,t_cnt), [4 3 1 2]), varargin{:});

                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = choice_opt_complex(blk, state2, action2, new_goal, ... 
                                'TransitionProb', permute(temp_T(:,:,:,t_cnt), [4 3 1 2]), varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
        case 'ChoOpt1n2MaxGoalLeftBias'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % not implemented yet
                    
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            max_goal = max(goal);
                            new_goal = zeros(size(goal));
%                             new_goal(max_goal(t_cnt)) = goal(max_goal(t_cnt));
                            new_goal(goal==max_goal) = goal(goal==max_goal);

                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, new_goal, blk, 'leftbias', varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, new_goal, blk, 'leftbias', varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
        case 'ChoOptMaxGoalA1LeftBias'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % not implemented yet
                    
                case 'Kim2019'
                    opt_likelihood1 = [];
                    for sess=1:N_session                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            max_goal = max(goal);
                            new_goal = zeros(size(goal));
%                             new_goal(max_goal(t_cnt)) = goal(max_goal(t_cnt));
                            new_goal(goal==max_goal) = goal(goal==max_goal);

                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, new_goal, blk, 'leftbias', varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                            
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood1; % [1 x N_trial]
            end
            
        case 'ChoOptMaxGoalA2LeftBias'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % not implemented yet
                    
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood2 = [];
                    for sess=1:N_session                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            max_goal = max(goal);
                            new_goal = zeros(size(goal));
%                             new_goal(max_goal(t_cnt)) = goal(max_goal(t_cnt));
                            new_goal(goal==max_goal) = goal(goal==max_goal);

                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, new_goal, blk, 'leftbias', varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood2; % [1 x N_trial]
            end
            
        case 'ChoOpt1n2AchGoalRan' 
            switch Exp
                
                case {'Kim2019', 'Kim2019 24'}
%                     max_goal = arbMBMF_load_var(Exp, 'Goal', id, []);
%                     ach_goal_s1 = arbMBMF_load_var(Exp, 'AchGoalS1RandA', id, []);
%                     ach_goal_s2 = arbMBMF_load_var(Exp, 'AchGoalS2RandA', id, []);
                    
                    t_cnt = 0; % trial count
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            t_cnt = t_cnt + 1;
                            blk = temp_blk_con(2, t_sess);
                            S2 = sess_behav_info(t_sess, 5);
                            goal = sess_behav_info(t_sess, 18:20);
                            ach1 = goal_distr_complex(blk, 1, 'TransitionStep', 2); % AchGoal 'Rand'
                            ach2 = goal_distr_complex(blk, S2);
                            maxgoal1 = max(ach1(1:3)'.*goal);
                            maxgoal2 = max(ach2(1:3)'.*goal);
                            goal1 = zeros(size(goal)); 
                            goal2 = zeros(size(goal));
                            
%                             goal1(ach_goal_s1(t_cnt)) = goal(ach_goal_s1(t_cnt));
                            goal1(ach1(1:3)'.*goal == maxgoal1) = goal(ach1(1:3)'.*goal == maxgoal1);
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal1, blk, 'likelihood', varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                            
%                             goal2(ach_goal_s2(t_cnt)) = goal(ach_goal_s2(t_cnt));
                            goal2(ach2(1:3)'.*goal == maxgoal2) = goal(ach2(1:3)'.*goal == maxgoal2);
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal2, blk, 'likelihood', varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
        case 'ChoOpt1n2AchGoalRanTransEst' 
            switch Exp
                
                case {'Kim2019', 'Kim2019 24'}
%                     max_goal = arbMBMF_load_var(Exp, 'Goal', id, []);
%                     ach_goal_s1 = arbMBMF_load_var(Exp, 'AchGoalS1RandA', id, []);
%                     ach_goal_s2 = arbMBMF_load_var(Exp, 'AchGoalS2RandA', id, []);
                    
                    temp_T = SBJ.('S1_state_fwd_T'); % N_state x N_action x N_state x N_trial
                    t_cnt = 0; % trial count
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            t_cnt = t_cnt + 1;
                            blk = temp_blk_con(2, t_sess);
                            S2 = sess_behav_info(t_sess, 5);
                            goal = sess_behav_info(t_sess, 18:20);
                            ach1 = goal_distr_complex(blk, 1, 'TransitionStep', 2); % AchGoal 'Rand'
                            ach2 = goal_distr_complex(blk, S2);
                            maxgoal1 = max(ach1(1:3)'.*goal);
                            maxgoal2 = max(ach2(1:3)'.*goal);
                            goal1 = zeros(size(goal)); 
                            goal2 = zeros(size(goal));
                            
%                             goal1(ach_goal_s1(t_cnt)) = goal(ach_goal_s1(t_cnt));
                            goal1(ach1(1:3)'.*goal == maxgoal1) = goal(ach1(1:3)'.*goal == maxgoal1);
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = choice_opt_complex(blk, state1, action1, goal1, ... 
                                'TransitionProb', permute(temp_T(:,:,:,t_cnt), [4 3 1 2]), varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                            
%                             goal2(ach_goal_s2(t_cnt)) = goal(ach_goal_s2(t_cnt));
                            goal2(ach2(1:3)'.*goal == maxgoal2) = goal(ach2(1:3)'.*goal == maxgoal2);
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = choice_opt_complex(blk, state2, action2, goal2, ... 
                                'TransitionProb', permute(temp_T(:,:,:,t_cnt), [4 3 1 2]), varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
        case 'ChoOptAchGoalRanA1LeftBias' 
            switch Exp
                
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood1 = [];
                    for sess=1:N_session                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            ach1 = goal_distr_complex(blk, 1, 'TransitionStep', 2); % AchGoal 'Rand'
                            maxgoal1 = max(ach1(1:3)'.*goal);
                            goal1 = zeros(size(goal)); 
                                                        
%                             goal1(ach_goal_s1(t_cnt)) = goal(ach_goal_s1(t_cnt));
                            goal1(ach1(1:3)'.*goal == maxgoal1) = goal(ach1(1:3)'.*goal == maxgoal1);
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal1, blk, 'leftbias', varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                            
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                                                
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood1; % [1 x N_trial]
            end
            
        case 'ChoOptAchGoalRanA2LeftBias' 
            switch Exp
                
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood2 = [];
                    for sess=1:N_session                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            S2 = sess_behav_info(t_sess, 5);
                            goal = sess_behav_info(t_sess, 18:20);
                            ach2 = goal_distr_complex(blk, S2);
                            maxgoal2 = max(ach2(1:3)'.*goal);
                            goal2 = zeros(size(goal));
                            
%                             goal2(ach_goal_s2(t_cnt)) = goal(ach_goal_s2(t_cnt));
                            goal2(ach2(1:3)'.*goal == maxgoal2) = goal(ach2(1:3)'.*goal == maxgoal2);
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal2, blk, 'leftbias', varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood2; % [1 x N_trial]
            end
        
        case 'ChoOpt1n2AchGoalOpt' 
            switch Exp
                
                case {'Kim2019', 'Kim2019 24'}
%                     max_goal = arbMBMF_load_var(Exp, 'Goal', id, []);
%                     ach_goal_s1 = arbMBMF_load_var(Exp, 'AchGoalS1RandA', id, []);
%                     ach_goal_s2 = arbMBMF_load_var(Exp, 'AchGoalS2RandA', id, []);
                    
                    t_cnt = 0; % trial count
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            t_cnt = t_cnt + 1;
                            blk = temp_blk_con(2, t_sess);
                            S2 = sess_behav_info(t_sess, 5);
                            goal = sess_behav_info(t_sess, 18:20);
                            ach1 = goal_distr_complex(blk, 1, ...
                                'TransitionStep', 2, ...
                                'Policy', 'optimal', ...
                                'GoalValues', goal); % AchGoal 'Opt'
                            ach2 = goal_distr_complex(blk, S2, ...
                                'Policy', 'optimal', ...
                                'GoalValues', goal); % AchGoal 'Opt'
                            maxgoal1 = max(ach1(1:3)'.*goal);
                            maxgoal2 = max(ach2(1:3)'.*goal);
                            goal1 = zeros(size(goal)); 
                            goal2 = zeros(size(goal));
                            
%                             goal1(ach_goal_s1(t_cnt)) = goal(ach_goal_s1(t_cnt));
                            goal1(ach1(1:3)'.*goal == maxgoal1) = goal(ach1(1:3)'.*goal == maxgoal1);
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal1, blk, 'likelihood', varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                            
%                             goal2(ach_goal_s2(t_cnt)) = goal(ach_goal_s2(t_cnt));
                            goal2(ach2(1:3)'.*goal == maxgoal2) = goal(ach2(1:3)'.*goal == maxgoal2);
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal2, blk, 'likelihood', varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
        case 'ChoOpt1n2AchGoalOptTransEst' 
            switch Exp
                
                case {'Kim2019', 'Kim2019 24'}
%                     max_goal = arbMBMF_load_var(Exp, 'Goal', id, []);
%                     ach_goal_s1 = arbMBMF_load_var(Exp, 'AchGoalS1RandA', id, []);
%                     ach_goal_s2 = arbMBMF_load_var(Exp, 'AchGoalS2RandA', id, []);
                    
                    temp_T = SBJ.('S1_state_fwd_T'); % N_state x N_action x N_state x N_trial
                    t_cnt = 0; % trial count
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            t_cnt = t_cnt + 1;
                            blk = temp_blk_con(2, t_sess);
                            S2 = sess_behav_info(t_sess, 5);
                            goal = sess_behav_info(t_sess, 18:20);
                            ach1 = goal_distr_complex(blk, 1, ...
                                'TransitionStep', 2, ...
                                'Policy', 'optimal', ...
                                'GoalValues', goal); % AchGoal 'Opt'
                            ach2 = goal_distr_complex(blk, S2, ...
                                'Policy', 'optimal', ...
                                'GoalValues', goal); % AchGoal 'Opt'
                            maxgoal1 = max(ach1(1:3)'.*goal);
                            maxgoal2 = max(ach2(1:3)'.*goal);
                            goal1 = zeros(size(goal)); 
                            goal2 = zeros(size(goal));
                            
%                             goal1(ach_goal_s1(t_cnt)) = goal(ach_goal_s1(t_cnt));
                            goal1(ach1(1:3)'.*goal == maxgoal1) = goal(ach1(1:3)'.*goal == maxgoal1);
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = choice_opt_complex(blk, state1, action1, goal1, ... 
                                'TransitionProb', permute(temp_T(:,:,:,t_cnt), [4 3 1 2]), varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                            
%                             goal2(ach_goal_s2(t_cnt)) = goal(ach_goal_s2(t_cnt));
                            goal2(ach2(1:3)'.*goal == maxgoal2) = goal(ach2(1:3)'.*goal == maxgoal2);
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = choice_opt_complex(blk, state2, action2, goal2, ... 
                                'TransitionProb', permute(temp_T(:,:,:,t_cnt), [4 3 1 2]), varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
        case 'ChoOptAchGoalOptA1LeftBias' 
            switch Exp
                
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood1 = [];
                    for sess=1:N_session                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            ach1 = goal_distr_complex(blk, 1, ...
                                'TransitionStep', 2, ...
                                'Policy', 'optimal', ...
                                'GoalValues', goal); % AchGoal 'Opt'
                            maxgoal1 = max(ach1(1:3)'.*goal);
                            goal1 = zeros(size(goal)); 
                            
%                             goal1(ach_goal_s1(t_cnt)) = goal(ach_goal_s1(t_cnt));
                            goal1(ach1(1:3)'.*goal == maxgoal1) = goal(ach1(1:3)'.*goal == maxgoal1);
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal1, blk, 'leftbias', varargin{:});
                            
                            opt_L1_sess(t_sess) = L1;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood1; % [1 x N_trial]
            end
            
        case 'ChoOptAchGoalOptA2LeftBias' 
            switch Exp
                
                case {'Kim2019', 'Kim2019 24'}
                    opt_likelihood2 = [];
                    for sess=1:N_session                        
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            S2 = sess_behav_info(t_sess, 5);
                            goal = sess_behav_info(t_sess, 18:20);
                            ach2 = goal_distr_complex(blk, S2, ...
                                'Policy', 'optimal', ...
                                'GoalValues', goal); % AchGoal 'Opt'
                            maxgoal2 = max(ach2(1:3)'.*goal);
                            goal2 = zeros(size(goal));
                            
%                             goal2(ach_goal_s2(t_cnt)) = goal(ach_goal_s2(t_cnt));
                            goal2(ach2(1:3)'.*goal == maxgoal2) = goal(ach2(1:3)'.*goal == maxgoal2);
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal2, blk, 'leftbias', varargin{:});
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = opt_likelihood2; % [1 x N_trial]
            end
            
        case 'ChoOpt1n2(uEarly)' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    blk_con = [];
                    for sess=1:N_session
                        blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                    UncPhase = class_series2phase_series(UC, 'PhaseNumber', 3);
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
                    Variable(:, ~ismember(UncPhase, 1)) = nan;
                    
                case {'Kim2019', 'Kim2019 24'}
            end
            
        case 'ChoOpt1(uEarly)' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = []; blk_con = [];
                    for sess=1:N_session
                        blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            
                            opt_L1_sess(t_sess) = L1;
                            
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                    UncPhase = class_series2phase_series(UC, 'PhaseNumber', 3);
                    Variable = opt_likelihood1;
                    Variable(:, ~ismember(UncPhase, 1)) = nan;
                    
                case {'Kim2019', 'Kim2019 24'}
                    
            end
            
        case 'ChoOpt2(uEarly)' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood2 = [];
                    blk_con = [];
                    for sess=1:N_session
                        blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                    UncPhase = class_series2phase_series(UC, 'PhaseNumber', 3);
                    Variable = opt_likelihood2; 
                    Variable(:, ~ismember(UncPhase, 1)) = nan;
                    
                case 'Kim2019'
            end
            
        case 'ChoOpt1n2(uLate)' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    blk_con = [];
                    for sess=1:N_session
                        blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                    UncPhase = class_series2phase_series(UC, 'PhaseNumber', 3);
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
                    Variable(:, ~ismember(UncPhase, 3)) = nan;
                    
                case 'Kim2019'
            end
            
        case 'ChoOpt1(uLate)' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = []; blk_con = [];
                    for sess=1:N_session
                        blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            
                            opt_L1_sess(t_sess) = L1;
                            
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                    UncPhase = class_series2phase_series(UC, 'PhaseNumber', 3);
                    Variable = opt_likelihood1;
                    Variable(:, ~ismember(UncPhase, 3)) = nan;
                    
                case 'Kim2019'
                    
            end
            
        case 'ChoOpt2(uLate)' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood2 = [];
                    blk_con = [];
                    for sess=1:N_session
                        blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                    UncPhase = class_series2phase_series(UC, 'PhaseNumber', 3);
                    Variable = opt_likelihood2; 
                    Variable(:, ~ismember(UncPhase, 3)) = nan;
                    
                case 'Kim2019'
            end
            
% -------------------------------------------------------------------------
% behavior: choice consistency
        case 'ChoConsist1n2'
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % ===== state action history =====
                    BHVs = [];  % for A1, S2, A2, G save
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        BHVs = [BHVs, BHVmat([7 5 8 18],:)];
                    end
                    A1 = BHVs(1,:);
                    S2 = BHVs(2,:);
                    A2 = BHVs(3,:);
                    G = BHVs(4,:);
                    
                    %             % ===== when tria_cond given =====
                    %             if ~isempty(tria_cond)
                    %                 if length(tria_cond)~=A1; error('tria_cond error'); end
                    %                 A1 = A1(tria_cond);
                    %                 S2 = S2(tria_cond);
                    %                 A2 = A2(tria_cond);
                    %             end
                    
                    ChoCon1 = NaN(1,length(A1));
                    ChoCon2 = NaN(1,length(A2));
                    Goal_set = [6 7 8 -1];
                    for gi = 1:length(Goal_set)
                        
                        g_idx = (G == Goal_set(gi));
                        
                        % ===== stage 1 choice consistency =====
                        A1_temp = A1(g_idx);
                        ChoCon1_temp = [0 diff(A1_temp)==0];
                        ChoCon1(g_idx) = ChoCon1_temp;
                        
                        % ===== stage 2 choice consistency =====
                        S2_temp = S2(g_idx);
                        A2_temp = A2(g_idx);
                        ChoCon2_temp = NaN(1,length(A2_temp));
                        S2_set = [2 3 4 5];
                        for si = 1:length(S2_set)
                            s2_idx = (S2_temp == S2_set(si));
                            if any(s2_idx)
                                a2 = A2_temp(s2_idx);
                                temp_consist2 = [0 diff(a2)==0];
                                ChoCon2_temp(s2_idx) = temp_consist2;
                            end
                        end
                        ChoCon2(g_idx) = ChoCon2_temp;
                        
                    end
                    
                case 'Kim2019'
                    error('(arbMBMF_load_var) ChoConsist1n2 was not implemented for Kim 2019')
                    
            end
            
            % sanity check
            if any(isnan(ChoCon2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with 2 rows =====
            Variable = [ChoCon1; ChoCon2];
            
        case 'ChoConsist1'
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % ===== state action history =====
                    BHVs = [];  % for A1, S2, A2, G save
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        BHVs = [BHVs, BHVmat([7 5 8 18],:)];
                    end
                    A1 = BHVs(1,:);
                    G = BHVs(4,:);
                    
                    ChoCon1 = NaN(1,length(A1));
                    Goal_set = [6 7 8 -1];
                    for gi = 1:length(Goal_set)
                        
                        g_idx = (G == Goal_set(gi));
                        
                        % ===== stage 1 choice consistency =====
                        A1_temp = A1(g_idx);
                        ChoCon1_temp = [0 diff(A1_temp)==0];
                        ChoCon1(g_idx) = ChoCon1_temp;
                        
                    end
                    
                case 'Kim2019'
                    error('(arbMBMF_load_var) ChoConsist1n2 was not implemented for Kim 2019')
                    
            end
            
            % ===== choice consistency row vector =====
            Variable = ChoCon1;
            
        case 'ChoConsist2'
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % ===== state action history =====
                    BHVs = [];  % for A1, S2, A2, G save
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        BHVs = [BHVs, BHVmat([7 5 8 18],:)];
                    end
                    S2 = BHVs(2,:);
                    A2 = BHVs(3,:);
                    G = BHVs(4,:);
                                        
                    ChoCon2 = NaN(1,length(A2));
                    Goal_set = [6 7 8 -1];
                    for gi = 1:length(Goal_set)
                        
                        g_idx = (G == Goal_set(gi));
                        
                        % ===== stage 2 choice consistency =====
                        S2_temp = S2(g_idx);
                        A2_temp = A2(g_idx);
                        ChoCon2_temp = NaN(1,length(A2_temp));
                        S2_set = [2 3 4 5];
                        for si = 1:length(S2_set)
                            s2_idx = (S2_temp == S2_set(si));
                            if any(s2_idx)
                                a2 = A2_temp(s2_idx);
                                temp_consist2 = [0 diff(a2)==0];
                                ChoCon2_temp(s2_idx) = temp_consist2;
                            end
                        end
                        ChoCon2(g_idx) = ChoCon2_temp;
                        
                    end
                    
                case 'Kim2019'
                    error('(arbMBMF_load_var) ChoConsist1n2 was not implemented for Kim 2019')
                    
            end
            
            % sanity check
            if any(isnan(ChoCon2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with 2 rows =====
            Variable = ChoCon2;
            
        case 'StrictChoCon'
            
            % ===== state action history =====
            A1 = []; S2 = []; A2 = []; G = [];
            for sess = 1:N_session
                BHVmat = SBJ.HIST_behavior_info{sess}';
                
                % sanity check
                if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                
                % S1 is always the same                
                A1 = [A1, BHVmat(7,:)];
                S2 = [S2, BHVmat(5,:)];
                A2 = [A2, BHVmat(8,:)];
                G = [G, BHVmat(18,:)];
            end
            
%             % ===== when tria_cond given =====
%             if ~isempty(tria_cond)
%                 if length(tria_cond)~=A1; error('tria_cond error'); end
%                 A1 = A1(tria_cond);
%                 S2 = S2(tria_cond);
%                 A2 = A2(tria_cond);
%             end
            
            ChoCon1 = NaN(1,length(A1));
            ChoCon2 = NaN(1,length(A2));
            Goal_set = [6 7 8 -1];
            for gi = 1:length(Goal_set)
                
                g_idx = (G == Goal_set(gi));
                
                % ===== stage 1 choice consistency =====
                A1_temp = A1(g_idx);
                ChoCon1_temp = [0 diff(A1_temp)==0];
                ChoCon1(g_idx) = ChoCon1_temp;
                
                % ===== stage 2 choice consistency =====
                S2_temp = S2(g_idx);
                A2_temp = A2(g_idx);
                ChoCon2_temp = NaN(1,length(A2_temp));
                S2_set = [2 3 4 5];
                for si = 1:length(S2_set)
                    s2_idx = (S2_temp == S2_set(si));
                    if any(s2_idx)
                        a2 = A2_temp(s2_idx);
                        temp_consist2 = [0 diff(a2)==0];
                        ChoCon2_temp(s2_idx) = temp_consist2;
                    end
                end
                ChoCon2(g_idx) = ChoCon2_temp;
                
            end
            
            % sanity check
            if any(isnan(ChoCon2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with 2 rows =====
            Variable = all([ChoCon1; ChoCon2], 1) + 0;
            
        case 'ChoCon ignoring context'
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    
                    % ===== state action history =====
                    A1 = []; S2 = []; A2 = []; % G = [];
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        
                        % S1 is always the same
                        A1 = [A1, BHVmat(7,:)];
                        S2 = [S2, BHVmat(5,:)];
                        A2 = [A2, BHVmat(8,:)];
                        %                 G = [G, BHVmat(18,:)];
                    end
                    
                    %             % ===== when tria_cond given =====
                    %             if ~isempty(tria_cond)
                    %                 if length(tria_cond)~=A1; error('tria_cond error'); end
                    %                 A1 = A1(tria_cond);
                    %                 S2 = S2(tria_cond);
                    %                 A2 = A2(tria_cond);
                    %             end
                    
                    % ===== stage 1 choice consistency =====
                    ChoCon1 = [0 diff(A1)==0];
                    
                    % ===== stage 2 choice consistency =====
                    ChoCon2 = NaN(1,length(A2));
                    S2_set = [2 3 4 5];
                    for si = 1:length(S2_set)
                        s2_idx = (S2 == S2_set(si));
                        if any(s2_idx)
                            a2_s2 = A2(s2_idx);
                            ChoCon2(s2_idx) = [0 diff(a2_s2)==0];
                        end
                    end
                    
                case 'Kim2019'
                    
                    % ===== state action history =====
                    A1 = []; S2 = []; A2 = []; % G = [];
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3]));       error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2 3 4]));       error('A2 error'); end
                        
                        % S1 is always the same
                        A1 = [A1, BHVmat(7,:)];
                        S2 = [S2, BHVmat(5,:)];
                        A2 = [A2, BHVmat(8,:)];
                        
                    end
                    
                    % ===== stage 1 choice consistency =====
                    ChoCon1 = [0 diff(A1)==0];
                    
                    % ===== stage 2 choice consistency =====
                    ChoCon2 = NaN(1,length(A2));
                    S2_set = [2 3];
                    for si = 1:length(S2_set)
                        s2_idx = (S2 == S2_set(si));
                        if any(s2_idx)
                            a2_s2 = A2(s2_idx);
                            ChoCon2(s2_idx) = [0 diff(a2_s2)==0];
                        end
                    end
                    
            end
            
            % sanity check
            if any(isnan(ChoCon2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with 2 rows =====
            Variable = [ChoCon1; ChoCon2];
            
        case 'ChoiceSwitch1n2' % based on ChoCon ignoring context
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    
                    % ===== state action history =====
                    A1 = []; S2 = []; A2 = []; % G = [];
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        
                        % S1 is always the same
                        A1 = [A1, BHVmat(7,:)];
                        S2 = [S2, BHVmat(5,:)];
                        A2 = [A2, BHVmat(8,:)];
                        %                 G = [G, BHVmat(18,:)];
                    end
                    
                    % ===== stage 1 choice switch =====
                    ChoSw1 = 1 - [0 diff(A1)==0];
%                     % ===== stage 1 choice consistency =====
%                     ChoCon1 = [0 diff(A1)==0];
                    
                    % ===== stage 2 choice consistency =====
                    ChoSw2 = NaN(1,length(A2));
                    S2_set = [2 3 4 5];
                    for si = 1:length(S2_set)
                        s2_idx = (S2 == S2_set(si));
                        if any(s2_idx)
                            a2_s2 = A2(s2_idx);
                            ChoSw2(s2_idx) = 1 - [0 diff(a2_s2)==0];
                        end
                    end
                    
%                     % ===== stage 2 choice consistency =====
%                     ChoCon2 = NaN(1,length(A2));
%                     S2_set = [2 3 4 5];
%                     for si = 1:length(S2_set)
%                         s2_idx = (S2 == S2_set(si));
%                         if any(s2_idx)
%                             a2_s2 = A2(s2_idx);
%                             ChoCon2(s2_idx) = [0 diff(a2_s2)==0];
%                         end
%                     end
                    
                case 'Kim2019'
                    
%                     % ===== state action history =====
%                     A1 = []; S2 = []; A2 = []; % G = [];
%                     for sess = 1:N_session
%                         BHVmat = SBJ.HIST_behavior_info{sess}';
%                         
%                         % sanity check
%                         if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
%                         if any(~ismember(BHVmat(5,:),[2 3]));       error('S2 error'); end
%                         if any(~ismember(BHVmat(8,:),[1 2 3 4]));       error('A2 error'); end
%                         
%                         % S1 is always the same
%                         A1 = [A1, BHVmat(7,:)];
%                         S2 = [S2, BHVmat(5,:)];
%                         A2 = [A2, BHVmat(8,:)];
%                         
%                     end
%                     
%                     % ===== stage 1 choice consistency =====
%                     ChoCon1 = [0 diff(A1)==0];
%                     
%                     % ===== stage 2 choice consistency =====
%                     ChoCon2 = NaN(1,length(A2));
%                     S2_set = [2 3];
%                     for si = 1:length(S2_set)
%                         s2_idx = (S2 == S2_set(si));
%                         if any(s2_idx)
%                             a2_s2 = A2(s2_idx);
%                             ChoCon2(s2_idx) = [0 diff(a2_s2)==0];
%                         end
%                     end
                    
            end
            
            % sanity check
            if any(isnan(ChoSw2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with 2 rows =====
            Variable = [ChoSw1; ChoSw2];
            
        case 'ChoiceSwitch1n2(GoalSwitch)' % based on ChoCon ignoring context
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    
                    % ===== state action history =====
                    A1 = []; S2 = []; A2 = []; G = [];
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        
                        % S1 is always the same
                        A1 = [A1, BHVmat(7,:)];
                        S2 = [S2, BHVmat(5,:)];
                        A2 = [A2, BHVmat(8,:)];
                        G = [G, BHVmat(18,:)];
                    end
                    
                    GS = [1 diff(G)~=0]; % goal switch
                    
                    % ===== stage 1 choice switch =====
                    ChoSw1 = 1 - [0 diff(A1)==0];
                    % ===== stage 2 choice consistency =====
                    ChoSw2 = NaN(1,length(A2));
                    S2_set = [2 3 4 5];
                    for si = 1:length(S2_set)
                        s2_idx = (S2 == S2_set(si));
                        if any(s2_idx)
                            a2_s2 = A2(s2_idx);
                            ChoSw2(s2_idx) = 1 - [0 diff(a2_s2)==0];
                        end
                    end
                                        
                case 'Kim2019'
                                        
            end
            
            % sanity check
            if any(isnan(ChoSw2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with 2 rows =====
            Variable = [ChoSw1; ChoSw2];
            Variable(:, ~ismember(GS, 1)) = nan;
            
        case 'ChoiceSwitch1n2(GoalStay)' % based on ChoCon ignoring context
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    
                    % ===== state action history =====
                    A1 = []; S2 = []; A2 = []; G = [];
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        
                        % S1 is always the same
                        A1 = [A1, BHVmat(7,:)];
                        S2 = [S2, BHVmat(5,:)];
                        A2 = [A2, BHVmat(8,:)];
                        G = [G, BHVmat(18,:)];
                    end
                    
                    GS = [1 diff(G)~=0]; % goal switch
                    
                    % ===== stage 1 choice switch =====
                    ChoSw1 = 1 - [0 diff(A1)==0];
                    % ===== stage 2 choice consistency =====
                    ChoSw2 = NaN(1,length(A2));
                    S2_set = [2 3 4 5];
                    for si = 1:length(S2_set)
                        s2_idx = (S2 == S2_set(si));
                        if any(s2_idx)
                            a2_s2 = A2(s2_idx);
                            ChoSw2(s2_idx) = 1 - [0 diff(a2_s2)==0];
                        end
                    end
                                        
                case 'Kim2019'
                                        
            end
            
            % sanity check
            if any(isnan(ChoSw2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with 2 rows =====
            Variable = [ChoSw1; ChoSw2];
            Variable(:, ~ismember(GS, 0)) = nan; % only different from ChoiceSwitch1n2(GoalSwitch)
            
        case 'ChoiceSwitch1' % based on ChoCon ignoring context
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    
                    % ===== state action history =====
                    A1 = []; S2 = []; A2 = []; % G = [];
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        
                        % S1 is always the same
                        A1 = [A1, BHVmat(7,:)];
                        S2 = [S2, BHVmat(5,:)];
                        A2 = [A2, BHVmat(8,:)];
                        %                 G = [G, BHVmat(18,:)];
                    end
                    
                    % ===== stage 1 choice switch =====
                    ChoSw1 = 1 - [0 diff(A1)==0];
%                     % ===== stage 1 choice consistency =====
%                     ChoCon1 = [0 diff(A1)==0];
                                        
                case 'Kim2019'
                    
            end
            % ===== choice consistency with a single row =====
            Variable = ChoSw1;
            
        case 'ChoiceSwitch1(GoalSwitch)' % based on ChoCon ignoring context
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    
                    % ===== state action history =====
                    A1 = []; S2 = []; A2 = []; G = [];
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        
                        % S1 is always the same
                        A1 = [A1, BHVmat(7,:)];
                        S2 = [S2, BHVmat(5,:)];
                        A2 = [A2, BHVmat(8,:)];
                        G = [G, BHVmat(18,:)];
                    end
                    GS = [1 diff(G)~=0]; % goal switch
                    
                    % ===== stage 1 choice switch =====
                    ChoSw1 = 1 - [0 diff(A1)==0];
                                        
                case 'Kim2019'
                    
            end
            % ===== choice consistency with a single row =====
            Variable = ChoSw1;
            Variable(:, ~ismember(GS, 1)) = nan; 
            
        case 'ChoiceSwitch1(GoalStay)' % based on ChoCon ignoring context
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    
                    % ===== state action history =====
                    A1 = []; S2 = []; A2 = []; G = [];
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        
                        % S1 is always the same
                        A1 = [A1, BHVmat(7,:)];
                        S2 = [S2, BHVmat(5,:)];
                        A2 = [A2, BHVmat(8,:)];
                        G = [G, BHVmat(18,:)];
                    end
                    GS = [1 diff(G)~=0]; % goal switch
                    
                    % ===== stage 1 choice switch =====
                    ChoSw1 = 1 - [0 diff(A1)==0];
                                        
                case 'Kim2019'
                    
            end
            % ===== choice consistency with a single row =====
            Variable = ChoSw1;
            Variable(:, ~ismember(GS, 0)) = nan; % only different from ChoiceSwitch1n2(GoalSwitch)
            
        case 'ChoiceSwitch2' % based on ChoCon ignoring context
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    
                    % ===== state action history =====
                    A1 = []; S2 = []; A2 = []; % G = [];
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        
                        % S1 is always the same
                        A1 = [A1, BHVmat(7,:)];
                        S2 = [S2, BHVmat(5,:)];
                        A2 = [A2, BHVmat(8,:)];
                        %                 G = [G, BHVmat(18,:)];
                    end
                    
                    % ===== stage 2 choice consistency =====
                    ChoSw2 = NaN(1,length(A2));
                    S2_set = [2 3 4 5];
                    for si = 1:length(S2_set)
                        s2_idx = (S2 == S2_set(si));
                        if any(s2_idx)
                            a2_s2 = A2(s2_idx);
                            ChoSw2(s2_idx) = 1 - [0 diff(a2_s2)==0];
                        end
                    end
                    
%                     % ===== stage 2 choice consistency =====
%                     ChoCon2 = NaN(1,length(A2));
%                     S2_set = [2 3 4 5];
%                     for si = 1:length(S2_set)
%                         s2_idx = (S2 == S2_set(si));
%                         if any(s2_idx)
%                             a2_s2 = A2(s2_idx);
%                             ChoCon2(s2_idx) = [0 diff(a2_s2)==0];
%                         end
%                     end
                    
                case 'Kim2019'
                    
            end
            
            % sanity check
            if any(isnan(ChoSw2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with a single row =====
            Variable = ChoSw2;
            
        case 'ChoiceSwitch2(GoalSwitch)' % based on ChoCon ignoring context
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    
                    % ===== state action history =====
                    A1 = []; S2 = []; A2 = []; G = [];
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        
                        % S1 is always the same
                        A1 = [A1, BHVmat(7,:)];
                        S2 = [S2, BHVmat(5,:)];
                        A2 = [A2, BHVmat(8,:)];
                        G = [G, BHVmat(18,:)];
                    end
                    GS = [1 diff(G)~=0]; % goal switch
                    
                    % ===== stage 2 choice consistency =====
                    ChoSw2 = NaN(1,length(A2));
                    S2_set = [2 3 4 5];
                    for si = 1:length(S2_set)
                        s2_idx = (S2 == S2_set(si));
                        if any(s2_idx)
                            a2_s2 = A2(s2_idx);
                            ChoSw2(s2_idx) = 1 - [0 diff(a2_s2)==0];
                        end
                    end
                    
                case 'Kim2019'
                    
            end
            
            % sanity check
            if any(isnan(ChoSw2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with a single row =====
            Variable = ChoSw2;
            Variable(:, ~ismember(GS, 1)) = nan; 
            
        case 'ChoiceSwitch2(GoalStay)' % based on ChoCon ignoring context
            
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    
                    % ===== state action history =====
                    A1 = []; S2 = []; A2 = []; G = [];
                    for sess = 1:N_session
                        BHVmat = SBJ.HIST_behavior_info{sess}';
                        
                        % sanity check
                        if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                        if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                        if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                        if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                        
                        % S1 is always the same
                        A1 = [A1, BHVmat(7,:)];
                        S2 = [S2, BHVmat(5,:)];
                        A2 = [A2, BHVmat(8,:)];
                        G = [G, BHVmat(18,:)];
                    end
                    GS = [1 diff(G)~=0]; % goal switch
                    
                    % ===== stage 2 choice consistency =====
                    ChoSw2 = NaN(1,length(A2));
                    S2_set = [2 3 4 5];
                    for si = 1:length(S2_set)
                        s2_idx = (S2 == S2_set(si));
                        if any(s2_idx)
                            a2_s2 = A2(s2_idx);
                            ChoSw2(s2_idx) = 1 - [0 diff(a2_s2)==0];
                        end
                    end
                    
                case 'Kim2019'
                    
            end
            
            % sanity check
            if any(isnan(ChoSw2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with a single row =====
            Variable = ChoSw2;
            Variable(:, ~ismember(GS, 0)) = nan; % only different from ChoiceSwitch1n2(GoalSwitch)
            
% -------------------------------------------------------------------------
% behavior: response time

        case 'RTA1'
            RTA1 = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};
                RTA1 = [RTA1, sess_behav_info(:,9)'];
            end
            Variable = RTA1;
            
        case 'RTA2'
            RTA2 = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};
                RTA2 = [RTA2, sess_behav_info(:,10)'];
            end
            Variable = RTA2;
            
        case 'RToptA1'
            RTA1 = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};
                RTA1 = [RTA1, sess_behav_info(:,9)'];
            end
            
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
                
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                sess_blkcon = sess_behav_info(:,3)';
                S2 = [S2, sess_behav_info(:,5)'];
                A2 = [A2, sess_behav_info(:,8)'];
                goal_state = [goal_state, sess_behav_info(:,18)'];
                
                if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0         % Sanity check
                    error('block condition information error'); end
            end
            
            % Finding choice optimality
            is_opt_choice = zeros(size(goal_state));
            
            % specific/low_unc/40goal
            idx_GL40_opt = [find(blk_con==1 & goal_state==6 & S2 == 4 & A2 == 2), ...
                find(blk_con==1 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/low_unc/20goal
            idx_GL20_opt = [find(blk_con==1 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==1 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/low_unc/10goal
            idx_GL10_opt = [find(blk_con==1 & goal_state==8 & S2 == 2 & A2 == 2), ...
                find(blk_con==1 & goal_state==8 & S2 == 3 & A2 == 1)];
            % specific/high_unc/40goal
            idx_GH40_opt = [find(blk_con==2 & goal_state==6 & S2 == 4 & A2 > 0), ...
                find(blk_con==2 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/high_unc/20goal
            idx_GH20_opt = [find(blk_con==2 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==2 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/high_unc/10goal
            idx_GH10_opt = [find(blk_con==2 & goal_state==8 & S2 == 2 & A2 > 0), ...
                find(blk_con==2 & goal_state==8 & S2 == 3 & A2 == 1)];
            % flexible/high_unc/uni_goal
            idx_FH_opt = [find(blk_con==3 & goal_state==-1 & S2 == 4 & A2 == 1), ...
                find(blk_con==3 & goal_state==-1 & S2 == 5 & A2 == 2)];
            % flexible/low_unc/uni_goal
            idx_FL_opt = [find(blk_con==4 & goal_state==-1 & S2 == 4 & A2 == 2), ...
                find(blk_con==4 & goal_state==-1 & S2 == 5 & A2 == 1)];
            
            is_opt_choice(idx_GL40_opt) = 1; is_opt_choice(idx_GH40_opt) = 1;
            is_opt_choice(idx_GL20_opt) = 1; is_opt_choice(idx_GH20_opt) = 1;
            is_opt_choice(idx_GL10_opt) = 1; is_opt_choice(idx_GH10_opt) = 1;
            is_opt_choice(idx_FH_opt) = 1; is_opt_choice(idx_FL_opt) = 1;
            
            Variable = RTA1;
            Variable(~is_opt_choice) = nan;
            
        case 'RToptA2'
            RTA2 = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};
                RTA2 = [RTA2, sess_behav_info(:,10)'];
            end
            
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
                
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                sess_blkcon = sess_behav_info(:,3)';
                S2 = [S2, sess_behav_info(:,5)'];
                A2 = [A2, sess_behav_info(:,8)'];
                goal_state = [goal_state, sess_behav_info(:,18)'];
                
                if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0         % Sanity check
                    error('block condition information error'); end
            end
            
            % Finding choice optimality
            is_opt_choice = zeros(size(goal_state));
            
            % specific/low_unc/40goal
            idx_GL40_opt = [find(blk_con==1 & goal_state==6 & S2 == 4 & A2 == 2), ...
                find(blk_con==1 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/low_unc/20goal
            idx_GL20_opt = [find(blk_con==1 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==1 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/low_unc/10goal
            idx_GL10_opt = [find(blk_con==1 & goal_state==8 & S2 == 2 & A2 == 2), ...
                find(blk_con==1 & goal_state==8 & S2 == 3 & A2 == 1)];
            % specific/high_unc/40goal
            idx_GH40_opt = [find(blk_con==2 & goal_state==6 & S2 == 4 & A2 > 0), ...
                find(blk_con==2 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/high_unc/20goal
            idx_GH20_opt = [find(blk_con==2 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==2 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/high_unc/10goal
            idx_GH10_opt = [find(blk_con==2 & goal_state==8 & S2 == 2 & A2 > 0), ...
                find(blk_con==2 & goal_state==8 & S2 == 3 & A2 == 1)];
            % flexible/high_unc/uni_goal
            idx_FH_opt = [find(blk_con==3 & goal_state==-1 & S2 == 4 & A2 == 1), ...
                find(blk_con==3 & goal_state==-1 & S2 == 5 & A2 == 2)];
            % flexible/low_unc/uni_goal
            idx_FL_opt = [find(blk_con==4 & goal_state==-1 & S2 == 4 & A2 == 2), ...
                find(blk_con==4 & goal_state==-1 & S2 == 5 & A2 == 1)];
            
            is_opt_choice(idx_GL40_opt) = 1; is_opt_choice(idx_GH40_opt) = 1;
            is_opt_choice(idx_GL20_opt) = 1; is_opt_choice(idx_GH20_opt) = 1;
            is_opt_choice(idx_GL10_opt) = 1; is_opt_choice(idx_GH10_opt) = 1;
            is_opt_choice(idx_FH_opt) = 1; is_opt_choice(idx_FL_opt) = 1;
            
            Variable = RTA2;
            Variable(~is_opt_choice) = nan;
                        
        case 'EarlyLateDeltaChoOptAfterBlkChange'
            % Len: # of blkChange + 1
            
        case 'TimeUntilFirstOptCho'
            % Len: # of blkChange + 1
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
                
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                sess_blkcon = sess_behav_info(:,3)';
                S2 = [S2, sess_behav_info(:,5)'];
                A2 = [A2, sess_behav_info(:,8)'];
                goal_state = [goal_state, sess_behav_info(:,18)'];
                
                if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0         % Sanity check
                    error('block condition information error'); end
            end
            
            % Finding choice optimality
            is_opt_choice = zeros(size(goal_state));
            
            % specific/low_unc/40goal
            idx_GL40_opt = [find(blk_con==1 & goal_state==6 & S2 == 4 & A2 == 2), ...
                find(blk_con==1 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/low_unc/20goal
            idx_GL20_opt = [find(blk_con==1 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==1 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/low_unc/10goal
            idx_GL10_opt = [find(blk_con==1 & goal_state==8 & S2 == 2 & A2 == 2), ...
                find(blk_con==1 & goal_state==8 & S2 == 3 & A2 == 1)];
            % specific/high_unc/40goal
            idx_GH40_opt = [find(blk_con==2 & goal_state==6 & S2 == 4 & A2 > 0), ...
                find(blk_con==2 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/high_unc/20goal
            idx_GH20_opt = [find(blk_con==2 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==2 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/high_unc/10goal
            idx_GH10_opt = [find(blk_con==2 & goal_state==8 & S2 == 2 & A2 > 0), ...
                find(blk_con==2 & goal_state==8 & S2 == 3 & A2 == 1)];
            % flexible/high_unc/uni_goal
            idx_FH_opt = [find(blk_con==3 & goal_state==-1 & S2 == 4 & A2 == 1), ...
                find(blk_con==3 & goal_state==-1 & S2 == 5 & A2 == 2)];
            % flexible/low_unc/uni_goal
            idx_FL_opt = [find(blk_con==4 & goal_state==-1 & S2 == 4 & A2 == 2), ...
                find(blk_con==4 & goal_state==-1 & S2 == 5 & A2 == 1)];
            
            is_opt_choice(idx_GL40_opt) = 1; is_opt_choice(idx_GH40_opt) = 1;
            is_opt_choice(idx_GL20_opt) = 1; is_opt_choice(idx_GH20_opt) = 1;
            is_opt_choice(idx_GL10_opt) = 1; is_opt_choice(idx_GH10_opt) = 1;
            is_opt_choice(idx_FH_opt) = 1; is_opt_choice(idx_FL_opt) = 1;
            
            % block change timing
            chng_idx = find([1 diff(blk_con)~=0]~=0);
            T_FirstOpt = nan * ones(size(chng_idx)); % LEN = # of blkChange
            for i = 1:length(chng_idx)
                changed = chng_idx(i);
                opt_after_chng = is_opt_choice(changed:end);
                T_FirstOpt(i) = find(opt_after_chng==1,1,'first');
            end
            Variable = T_FirstOpt;

            
% -------------------------------------------------------------------------
% RL variable: state prediction error

        case 'SPE2' % SPE in stage 2
            SPE = SBJ.regressor{1,1}.value(7,:);
            Variable = SPE(2:3:end);
        case 'SPE3' % SPE in stage 3
            SPE = SBJ.regressor{1,1}.value(7,:);
            Variable = SPE(3:3:end);
        case 'SPE' % average SPE in stage 2 & 3
            SPE = SBJ.regressor{1,1}.value(7,:);
            Variable = SPE(2:3:end) + SPE(3:3:end) / 2;

% -------------------------------------------------------------------------
% RL variable: reward prediction error

        case 'RPE2' % SPE in stage 2
            RPE = SBJ.regressor{1,2}.value(7,:);
            Variable = RPE(2:3:end);
        case 'RPE3' % SPE in stage 3
            RPE = SBJ.regressor{1,2}.value(7,:);
            Variable = RPE(3:3:end);
        case 'RPE' % average SPE in stage 2 & 3
            RPE = SBJ.regressor{1,2}.value(7,:);
            Variable = RPE(2:3:end) + RPE(3:3:end) / 2;

% -------------------------------------------------------------------------
% RL variable: reliability

        case 'relMB2' % relMB in stage 2
            relMB = SBJ.regressor{1,7}.value(8,:)./...
                sum(SBJ.regressor{1,7}.value(7:9,:), 1);
            Variable = relMB(2:3:end);
        case 'relMB3' % relMB in stage 3
            relMB = SBJ.regressor{1,7}.value(8,:)./...
                sum(SBJ.regressor{1,7}.value(7:9,:), 1);
            Variable = relMB(3:3:end);
        case 'relMB' % average relMB in stage 2 & 3
            relMB = SBJ.regressor{1,7}.value(8,:)./...
                sum(SBJ.regressor{1,7}.value(7:9,:), 1);
            Variable = relMB(2:3:end) + relMB(3:3:end) / 2;
        case 'relMB_G' % average relMB in stage 2 & 3
            relMB = SBJ.regressor{1,7}.value(8,:)./...
                sum(SBJ.regressor{1,7}.value(7:9,:), 1);
            Variable = relMB(2:3:end) + relMB(3:3:end) / 2;
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case 'Lee2014'
                    Variable(ismember(blk_con, [3 4])) = nan;
            end            
        
        case 'relMF2' % relMF in stage 2
            relMF = SBJ.regressor{1,8}.value(8,:);
            Variable = relMF(2:3:end);
        case 'relMF3' % relMF in stage 3
            relMF = SBJ.regressor{1,8}.value(8,:);
            Variable = relMF(3:3:end);
        case 'relMF' % average relMF in stage 2 & 3
            relMF = SBJ.regressor{1,8}.value(8,:);
            Variable = relMF(2:3:end) + relMF(3:3:end) / 2;
            
        case 'relMAX2' % relMAX in stage 2
            switch Exp
                case 'Lee2014'
                    relMAX = nanmax([relMB; relMF], [], 1);
                    Variable = relMAX(2:3:end);
                case 'Kim2019'
                    relMAX = SBJ.regressor{1,22}.value(8,:); % MAXinvFano12
                    Variable = relMAX(2:3:end);
            end            
        case 'relMAX3' % relMAX in stage 3
            switch Exp
                case 'Lee2014'
                    relMAX = nanmax([relMB; relMF], [], 1);
                    Variable = relMAX(3:3:end);
                case 'Kim2019'
                    relMAX = SBJ.regressor{1,22}.value(8,:); % MAXinvFano12
                    Variable = relMAX(3:3:end);
            end            
        case 'relMAX' % average relMAX in stage 2 & 3
            switch Exp
                case 'Lee2014'
                    relMAX = nanmax([relMB; relMF], [], 1);
                    Variable = relMAX(2:3:end) + relMAX(3:3:end) / 2;
                case 'Kim2019'
                    relMAX = SBJ.regressor{1,22}.value(8,:); % MAXinvFano12
                    Variable = relMAX(2:3:end) + relMAX(3:3:end) / 2;
            end             
        case 'relMAX_G'
            relMAX = nanmax([relMB; relMF], [], 1);
            Variable = relMAX(2:3:end) + relMAX(3:3:end) / 2;
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case 'Lee2014'
                    Variable(ismember(blk_con, [3 4])) = nan;
            end
            
        case 'relDiff2'
            relDiff = relMB - relMF;
            Variable = relDiff(2:3:end);
        case 'relDiff3'
            relDiff = relMB - relMF;
            Variable = relDiff(3:3:end);
        case 'relDiff'
            relDiff = relMB - relMF;
            Variable = relDiff(2:3:end) + relDiff(3:3:end) / 2;
        case 'relDiff_G'
            relDiff = relMB - relMF;
            Variable = relDiff(2:3:end) + relDiff(3:3:end) / 2;
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case 'Lee2014'
                    Variable(ismember(blk_con, [3 4])) = nan;
            end
            
        case 'relComp2'
            relComp = (relMB > relMF) + 0;
            Variable = relComp(2:3:end);
        case 'relComp3'
            relComp = (relMB > relMF) + 0;
            Variable = relComp(3:3:end);
        case 'relComp'
            relComp = (relMB > relMF) + 0;
            Variable = (relComp(2:3:end) + relComp(3:3:end)) / 2;
            
% -------------------------------------------------------------------------
% RL variable: P_MB

        case 'PMB2'  % PMB
            PMB = SBJ.regressor{1,9}.value(7,:); % weigtM1 value (Lee2014, Kim2019)
            Variable = PMB(2:3:end);
        case 'PMB3'  % PMB
            PMB = SBJ.regressor{1,9}.value(7,:); % weigtM1 value (Lee2014, Kim2019)
            Variable = PMB(3:3:end);
        case 'PMB'  % PMB
            PMB = SBJ.regressor{1,9}.value(7,:); % weigtM1 value (Lee2014, Kim2019)
            Variable = (PMB(1:3:end)+PMB(2:3:end)+PMB(3:3:end))/3;
        case 'PMB_G'
            PMB = SBJ.regressor{1,9}.value(7,:);
            Variable = (PMB(1:3:end)+PMB(2:3:end)+PMB(3:3:end))/3;
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case 'Lee2014'
                    Variable(ismember(blk_con, [3 4])) = nan;
            end

% -------------------------------------------------------------------------
% RL variable: state-action value
            
        case 'Qarb1' % Qarb(chosen - unchosen) in stage 1 (Lee2014, Kim2019)
            Qarb = SBJ.regressor{1,13}.value(7,:);
            Variable = Qarb(1:3:end);
        case 'Qarb2' % Qarb(chosen - unchosen) in stage 2 (Lee2014, Kim2019)
            Qarb = SBJ.regressor{1,13}.value(7,:);
            Variable = Qarb(2:3:end);
        case 'Qarb' % average Qarb(chosen - unchosen) in stage 1 & 2 (Lee2014, Kim2019)
            Qarb = SBJ.regressor{1,13}.value(7,:);
            Variable = Qarb(1:3:end) + Qarb(2:3:end) / 2;
            
% -------------------------------------------------------------------------
% RL variable: choice probability (likelihood)

        case 'A1NLL' % negative log likelihood of chosen action in ther 1st stage
            NLL = SBJ.model_error{1,3}.value(7,:);
            Variable = NLL(1:3:end);
        case 'A2NLL' % negative log likelihood of chosen action in ther 2nd stage
            NLL = SBJ.model_error{1,3}.value(7,:);
            Variable = NLL(2:3:end);
        case 'A1NLLfwd' % negative log likelihood of chosen action in ther 1st stage
            NLL = SBJ.model_error{1,1}.value(7,:);
            Variable = NLL(1:3:end);
        case 'A2NLLfwd' % negative log likelihood of chosen action in ther 2nd stage
            NLL = SBJ.model_error{1,1}.value(7,:);
            Variable = NLL(2:3:end);
        case 'A1NLLsarsa' % negative log likelihood of chosen action in ther 1st stage
            NLL = SBJ.model_error{1,2}.value(7,:);
            Variable = NLL(1:3:end);
        case 'A2NLLsarsa' % negative log likelihood of chosen action in ther 2nd stage
            NLL = SBJ.model_error{1,2}.value(7,:);
            Variable = NLL(2:3:end);
            
% -------------------------------------------------------------------------
% RL parameters

        case 'ThrSPE' % zero-SPE threshold
            Variable = SBJ.model_BayesArb.param(1); % Lee 2014 & Kim 2019
        case 'lrRPE' % learning rate for |RPE|
            Variable = SBJ.model_BayesArb.param(2); % Lee 2014 & Kim 2019
        case 'MB2MF' % transition rate (MB->MF)
            Variable = SBJ.model_BayesArb.param(3); % Lee 2014 & Kim 2019
        case 'MF2MB' % transition rate (MB->MF)
            Variable = SBJ.model_BayesArb.param(5); % Lee 2014 & Kim 2019
        case 'invTemp' % inverse softmax temperature
            Variable = SBJ.model_BayesArb.param(7); % Lee 2014 & Kim 2019
        case 'lr'
            Variable = SBJ.model_BayesArb.param(8); % Lee 2014 supp
                        
    end
    
% -------------------------------------------------------------------------

    % ****** when id ~= 0 ******

    if ~isempty(tria_cond) && ~isscalar(Variable)
        if ischar(tria_cond)
            switch tria_cond
                case 'S1'
                    session = arbMBMF_load_var(Exp, 'Session', id, []);
                    Variable(:, ~ismember(session, 1)) = nan;
                case 'S2'
                    session = arbMBMF_load_var(Exp, 'Session', id, []);
                    Variable(:, ~ismember(session, 2)) = nan;
                case 'S3'
                    session = arbMBMF_load_var(Exp, 'Session', id, []);
                    Variable(:, ~ismember(session, 3)) = nan;
                case 'S4'
                    session = arbMBMF_load_var(Exp, 'Session', id, []);
                    Variable(:, ~ismember(session, 4)) = nan;
                case 'S5'
                    session = arbMBMF_load_var(Exp, 'Session', id, []);
                    Variable(:, ~ismember(session, 5)) = nan;
                case 'state2'
                    BHVcell = SBJ.HIST_behavior_info; % [1 x N_session]
                    BHVs = cell2mat(BHVcell'); % [N_trial x N_measures]
                    S2 = BHVs(:, 5)'; % 2 3 4 5, row vector
                    Variable(:, ~ismember(S2, 2)) = nan;
                case 'state3'
                    BHVcell = SBJ.HIST_behavior_info; % [1 x N_session]
                    BHVs = cell2mat(BHVcell'); % [N_trial x N_measures]
                    S2 = BHVs(:, 5)'; % 2 3 4 5, row vector
                    Variable(:, ~ismember(S2, 3)) = nan;
                case 'state4'
                    BHVcell = SBJ.HIST_behavior_info; % [1 x N_session]
                    BHVs = cell2mat(BHVcell'); % [N_trial x N_measures]
                    S2 = BHVs(:, 5)'; % 2 3 4 5, row vector
                    Variable(:, ~ismember(S2, 4)) = nan;
                case 'state5'
                    BHVcell = SBJ.HIST_behavior_info; % [1 x N_session]
                    BHVs = cell2mat(BHVcell'); % [N_trial x N_measures]
                    S2 = BHVs(:, 5)'; % 2 3 4 5, row vector
                    Variable(:, ~ismember(S2, 5)) = nan;
                case 'G' % for Lee 2014
                    BlkMat = cell2mat(SBJ.HIST_block_condition);
                    context = BlkMat(2,:);
                    Variable(:, ~ismember(context, [1 2])) = nan;
                case 'H'
                    BlkMat = cell2mat(SBJ.HIST_block_condition);
                    context = BlkMat(2,:);
%                     GT = arbMBMF_load_var(Exp, 'GoalCond', id, []); % 1:G, 2:H
                    Variable(:, ~ismember(context, [3 4])) = nan;
                case 'H + Red'
                    switch Exp
                        case {'Lee2014', 'Heo2018'}
                            % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                            Goal = arbMBMF_load_var(Exp, 'Goal', id, []); 
                            Variable(:, ~ismember(Goal, [6, -1])) = nan;
                        otherwise
                            error('(arbMBMF_load_var.m) not implemented for Kim2019')
                    end
                case 'Goal(6,-1)'
                    switch Exp
                        case {'Lee2014', 'Heo2018'}
                            % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                            Goal = arbMBMF_load_var(Exp, 'Goal', id, []); 
                            Variable(:, ~ismember(Goal, [6, -1])) = nan;
                        otherwise
                            error('(arbMBMF_load_var.m) not implemented for Kim2019')
                    end
                case 'Goal(6,8)'
                    switch Exp
                        case {'Lee2014', 'Heo2018'}
                            % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                            Goal = arbMBMF_load_var(Exp, 'Goal', id, []); 
                            Variable(:, ~ismember(Goal, [6, 8])) = nan;
                        otherwise
                            error('(arbMBMF_load_var.m) not implemented for Kim2019')
                    end
                case 'Goal6'
                    switch Exp
                        case {'Lee2014', 'Heo2018'}
                            % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                            Goal = arbMBMF_load_var(Exp, 'Goal', id, []); 
                            Variable(:, ~ismember(Goal, 6)) = nan;
                        otherwise
                            error('(arbMBMF_load_var.m) not implemented for Kim2019')
                    end
                case 'Goal7'
                    switch Exp
                        case {'Lee2014', 'Heo2018'}
                            % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                            Goal = arbMBMF_load_var(Exp, 'Goal', id, []); 
                            Variable(:, ~ismember(Goal, 7)) = nan;
                        otherwise
                            error('(arbMBMF_load_var.m) not implemented for Kim2019')
                    end
                case 'Goal8'
                    switch Exp
                        case {'Lee2014', 'Heo2018'}
                            % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                            Goal = arbMBMF_load_var(Exp, 'Goal', id, []); 
                            Variable(:, ~ismember(Goal, 8)) = nan;
                        otherwise
                            error('(arbMBMF_load_var.m) not implemented for Kim2019')
                    end
                case 'G uL' % for Lee 2014
                    BlkMat = cell2mat(SBJ.HIST_block_condition);
                    context = BlkMat(2,:);
                    Variable(:, context~=1) = nan;
                case 'G uH' % for Lee 2014
                    BlkMat = cell2mat(SBJ.HIST_block_condition);
                    context = BlkMat(2,:);
                    Variable(:, context~=2) = nan;
                case 'uL'
                    UC = arbMBMF_load_var(Exp, 'UncCond', id, []); % 1:low, 2:high
                    Variable(:, UC~=1) = nan;
                case 'uH'
                    UC = arbMBMF_load_var(Exp, 'UncCond', id, []); % 1:low, 2:high
                    Variable(:, UC~=2) = nan;
                case 'binPMB0'
                    PMB = arbMBMF_load_var(Exp, 'PMB', id, []);
                    binPMB = regr2class(PMB, 2);
                    Variable(:, binPMB~=0) = nan;
%                     disp('test')
                case 'binPMB1'
                    PMB = arbMBMF_load_var(Exp, 'PMB', id, []);
                    binPMB = regr2class(PMB, 2);
                    Variable(:, binPMB~=1) = nan;
                case 'binPMB1 uL'
                    PMB = arbMBMF_load_var(Exp, 'PMB', id, []); binPMB = regr2class(PMB, 2);
                    UC = arbMBMF_load_var(Exp, 'UncCond', id, []); % 1:low, 2:high
                    Variable(:, binPMB~=1) = nan; Variable(:, UC~=1) = nan; 
                case 'binPMB1 uH'
                    PMB = arbMBMF_load_var(Exp, 'PMB', id, []); binPMB = regr2class(PMB, 2);
                    UC = arbMBMF_load_var(Exp, 'UncCond', id, []); % 1:low, 2:high
                    Variable(:, binPMB~=1) = nan; Variable(:, UC~=2) = nan; 
                case 'uLate'
                    UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, []); % 0:early, 1:late
                    Variable(:, UncPhaseEL~=1) = nan;
                case 'uEarly'
                    UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, []); % 0:early, 1:late
                    Variable(:, UncPhaseEL~=0) = nan;
                case 'G uLate'
                    GT = arbMBMF_load_var(Exp, 'GoalCond', id, []); % 1:G, 2:H
                    UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, []); % 0:early, 1:late
                    Variable(:, GT~=1) = nan; Variable(:, UncPhaseEL~=1) = nan;
                case 'G uEarly'
                    GT = arbMBMF_load_var(Exp, 'GoalCond', id, []); % 1:G, 2:H
                    UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, []); % 0:early, 1:late
                    Variable(:, GT~=1) = nan; Variable(:, UncPhaseEL~=0) = nan;
                case 'H uLate'
                    GT = arbMBMF_load_var(Exp, 'GoalCond', id, []); % 1:G, 2:H
                    UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, []); % 0:early, 1:late
                    Variable(:, GT~=2) = nan; Variable(:, UncPhaseEL~=1) = nan;
                case 'H uEarly'
                    GT = arbMBMF_load_var(Exp, 'GoalCond', id, []); % 1:G, 2:H
                    UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, []); % 0:early, 1:late
                    Variable(:, GT~=2) = nan; Variable(:, UncPhaseEL~=0) = nan;
                case 'binPMB1 uLate'
                    PMB = arbMBMF_load_var(Exp, 'PMB', id, []); binPMB = regr2class(PMB, 2);
                    UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, []); % 0:early, 1:late
                    Variable(:, binPMB~=1) = nan; Variable(:, UncPhaseEL~=1) = nan;
                case 'binPMB1 uEarly'
                    PMB = arbMBMF_load_var(Exp, 'PMB', id, []); binPMB = regr2class(PMB, 2);
                    UncPhaseEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, []); % 0:early, 1:late
                    Variable(:, binPMB~=1) = nan; Variable(:, UncPhaseEL~=0) = nan;
                case 'Goal(6,7,8)xUCSwitch'
                    Switch = arbMBMF_load_var(Exp, 'Goal(6,7,8)xUCSwitch', id, []);
                    Variable(:, Switch~=1) = nan;
                case 'Goal(6,7,8)xUCStay'
                    Switch = arbMBMF_load_var(Exp, 'Goal(6,7,8)xUCSwitch', id, []);
                    Variable(:, Switch~=0) = nan;
            end
        elseif isscalar(tria_cond)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    switch tria_cond
                        case 1
                            blkCond = arbMBMF_load_var(Exp, 'blkCond', id, []);
                            Variable(:, ismember(blkCond, [3, 4])) = nan;
%                             fprintf('t1 ')
                        case 2
                            blkCond = arbMBMF_load_var(Exp, 'blkCond', id, []);
                            Variable(:, ismember(blkCond, [1, 2])) = nan;
                    end
                case 'Kim2019'
                    disp('(arbMBMF_load_var) not implemented yet');
            end
        else
            if isvector(Variable)
                if any(isnan(tria_cond))
                    Variable(isnan(tria_cond)) = nan;
                else
                    Variable(tria_cond==0) = nan;
                end
            else
                if any(isnan(tria_cond))
                    Variable(:, isnan(tria_cond)) = nan;
                else
                    Variable(:, tria_cond==0) = nan;
                end
            end
        end
    end
    
end

