function [SBJ]=ArbBat_Ori(maxi,mode,PreBehav,PreBlck,MainBehav,MainBlck,list_sbj,param_init)

mode.param_init_ori = param_init;
for i = 1 : maxi
    fprintf('###   SUB_NUM: [%d / %d]\n',i,maxi);
    fprintf('### OPT_ITER : [%d]\n',mode.max_iter);
    disp('############################################');
    disp('############################################');
    disp('############################################');
    
    data_in{1,1}.map_type=1;

    % pre save
    data_pre=load([list_sbj{i} '_pre_1.mat']);
    data_in{1,1}.HIST_behavior_info_pre{1,1}=PreBehav{i};
    data_in{1,1}.HIST_block_condition_pre{1,1}=PreBlck{i};

    % main save
    tt = dir('functions/analysis_simulation/behavior_record');
    tt = {tt.name};
    maxsess = sum(cell2mat(strfind(tt,[list_sbj{i} '_fmri_']))) - 1;
    for ii = 1 : maxsess
        data_in{1,1}.HIST_behavior_info{1,ii} = MainBehav{i,ii};
        data_in{1,1}.HIST_block_condition{1,ii} = MainBlck{i,ii};
    end
        
    %         % NO opt process
    %         outvalval = eval_ArbitrationRL6c(param_init, data_in, mode);
    %         outputval{d}=[outputval{d} outvalval];
    
    
    % optimization part
    myFunc_bu = @(x) eval_ArbitrationRL3c(x, data_in, mode);
    disp(['    ***** subject number : [' num2str(i) '], subject name : [' list_sbj{i} ']' ]);
    
    if mode.out == 1 % parameter fitting to human behavior data
        
        [model_BayesArb.param, model_BayesArb.val] = fminsearchbnd(myFunc_bu, ...
            param_init, mode.param_BoundL, mode.param_BoundU, optimset('Display', 'final','MaxIter',mode.max_iter)); % X0,LB,UB
        mode_temp = mode;
        mode_temp.out = 2;
        model_BayesArb.episode = eval_ArbitrationRL3c(model_BayesArb.param, data_in, mode_temp);

    elseif mode.out == 2 % agent episode generation

        model_BayesArb.episode = eval_ArbitrationRL3c(param_init, data_in, mode);

    elseif mode.out == 99 % regressor generation with fixed parameters

        data_in = eval_ArbitrationRL3c(param_init, data_in, mode);    

    end
    model_BayesArb.mode = mode;
    % for Storing
    SBJ{1,i} = data_in{1,1};
    SBJ{1,i}.model_BayesArb = model_BayesArb;
    disp('############################################');
    disp('############################################');
    disp('############################################');
end
end