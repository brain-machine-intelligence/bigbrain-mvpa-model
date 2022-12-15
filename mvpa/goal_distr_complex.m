function d = goal_distr_complex(context, state, varargin)

N_c = 4;    % N_context
N_s = 11;   % N_state
N_a = 4;    % N_action

c = context;
s = state;

goal1_id = [5 6]; % silver coin
goal2_id = [7 9]; % red coin
goal3_id = [8 11]; % blue coin
goal4_id = [4 10]; % nothing
d = zeros(4,1);

options = struct('Policy', 'random', ...
    'TransitionStep', 1, ...
    'GoalValues', []);
% varargin validity
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(goal_distr_complex) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(goal_distr_complex) %s is not a recognized parameter name', pair{1})
    end
end

pi = zeros(N_c, N_s, N_a);    % pi_c(s,a), context-dependent policy

switch options.Policy
    case 'random' % random policy
        pi(:,1,1:2) = .5;
        pi(1,2:3,1:2) = .5; % random policy among L1 & R1
        pi(2,2:3,1:4) = .25; % random policy among L1 ~ R2
        pi(3,2:3,1:2) = .5; % random policy among L1 & R1
        pi(4,2:3,1:4) = .25; % random policy among L1 ~ R2
    case 'optimal' % optimal policy based on the given goal values
        if isempty(options.GoalValues); error('(goal_distr_complex) no goal values'); end
        goalVals = options.GoalValues;
end

%%%%%%%% T: context-dependent state-action-state transition matrix %%%%%%%%
T = zeros(N_c, N_s', N_s, N_a); % T_c(s'|s,a), a = L1, R1, L2, R2 respectively

p = .9;
q = .1;
r = .5;

% ----- lowUC, lowCX ----- 
% T(s'|s,L1)
             %1 2 3 4 5 6 7 8 9 0 1
T(1,:,:,1) = [0 0 0 0 0 0 0 0 0 0 0;  % 1
              q 0 0 0 0 0 0 0 0 0 0;  % 2
              p 0 0 0 0 0 0 0 0 0 0;  % 3
              0 p 0 0 0 0 0 0 0 0 0;  % 4 
              0 0 0 0 0 0 0 0 0 0 0;  % 5
              0 0 p 0 0 0 0 0 0 0 0;  % 6
              0 0 0 0 0 0 0 0 0 0 0;  % 7
              0 q 0 0 0 0 0 0 0 0 0;  % 8
              0 0 0 0 0 0 0 0 0 0 0;  % 9
              0 0 q 0 0 0 0 0 0 0 0;  % 10
              0 0 0 0 0 0 0 0 0 0 0]; % 11

% T(s'|s,R1)
             %1 2 3 4 5 6 7 8 9 0 1
T(1,:,:,2) = [0 0 0 0 0 0 0 0 0 0 0;  % 1
              p 0 0 0 0 0 0 0 0 0 0;  % 2
              q 0 0 0 0 0 0 0 0 0 0;  % 3
              0 0 0 0 0 0 0 0 0 0 0;  % 4
              0 p 0 0 0 0 0 0 0 0 0;  % 5
              0 0 0 0 0 0 0 0 0 0 0;  % 6
              0 0 p 0 0 0 0 0 0 0 0;  % 7
              0 0 0 0 0 0 0 0 0 0 0;  % 8
              0 q 0 0 0 0 0 0 0 0 0;  % 9
              0 0 0 0 0 0 0 0 0 0 0;  % 10
              0 0 q 0 0 0 0 0 0 0 0]; % 11

% ----- lowUC, highCX ----- 
% T(s'|s,L1)
T(2,:,:,1) = zeros(11);
T(2,3,1,1) = p; % T(3|1,L1)
T(2,2,1,1) = q; % T(2|1,L1)
T(2,4,2,1) = p; % T(4|2,L1)
T(2,8,2,1) = q; % T(8|2,L1)
T(2,6,3,1) = p; 
T(2,10,3,1) = q; 

% T(s'|s,R1)
T(2,:,:,2) = zeros(11);
T(2,2,1,2) = p; 
T(2,3,1,2) = q; 
T(2,5,2,2) = p; 
T(2,9,2,2) = q; 
T(2,7,3,2) = p; 
T(2,11,3,2) = q; 

% T(s'|s,L2)
T(2,:,:,3) = zeros(11);
T(2,4,2,3) = q; 
T(2,8,2,3) = p; 
T(2,6,3,3) = q; 
T(2,10,3,3) = p; 

% T(s'|s,R2)
T(2,:,:,4) = zeros(11);
T(2,5,2,4) = q; 
T(2,9,2,4) = p; 
T(2,7,3,4) = q; 
T(2,11,3,4) = p; 

% ----- highUC, lowCX ----- 
% T(s'|s,L1)
T(3,:,:,1) = zeros(11);
T(3,3,1,1) = r; 
T(3,2,1,1) = r; 
T(3,4,2,1) = r; 
T(3,8,2,1) = r; 
T(3,6,3,1) = r; 
T(3,10,3,1) = r; 

% T(s'|s,R1)
T(3,:,:,2) = zeros(11);
T(3,2,1,2) = r; 
T(3,3,1,2) = r; 
T(3,5,2,2) = r; 
T(3,9,2,2) = r; 
T(3,7,3,2) = r; 
T(3,11,3,2) = r; 

% ----- lowUC, highCX ----- 
% T(s'|s,L1)
T(4,:,:,1) = zeros(11);
T(4,3,1,1) = r; 
T(4,2,1,1) = r; 
T(4,4,2,1) = r; 
T(4,8,2,1) = r; 
T(4,6,3,1) = r; 
T(4,10,3,1) = r; 

% T(s'|s,R1)
T(4,:,:,2) = zeros(11);
T(4,2,1,2) = r; 
T(4,3,1,2) = r; 
T(4,5,2,2) = r; 
T(4,9,2,2) = r; 
T(4,7,3,2) = r; 
T(4,11,3,2) = r; 

% T(s'|s,L2)
T(4,:,:,3) = zeros(11);
T(4,4,2,3) = r; 
T(4,8,2,3) = r; 
T(4,6,3,3) = r; 
T(4,10,3,3) = r; 

% T(s'|s,R2)
T(4,:,:,4) = zeros(11);
T(4,5,2,4) = r; 
T(4,9,2,4) = r; 
T(4,7,3,4) = r; 
T(4,11,3,4) = r; 

% =====================================================

% T_c(s'|s,pi) = sum_a [pi(a|s) * T_c(s'|s,a)]
switch options.Policy
    case 'random'
        temp_pi = squeeze(pi(c,:,:));
    case 'optimal'
        r_s3_vec = zeros(N_s, 1); % column
        r_s3_vec(goal1_id) = goalVals(1);
        r_s3_vec(goal2_id) = goalVals(2);
        r_s3_vec(goal3_id) = goalVals(3);
        r_s3_vec(goal4_id) = 0;
        
        temp_pi = zeros(N_s, N_a);
        for si = 2:3
            % stage 2 - Q(s_fix, a) vector [N_a x 1]
            q_a2_vec = r_s3_vec' * squeeze(T(c, :, si, :)); % [1 x N_s] * [N_s' x N_a]
            q_a2_vec = q_a2_vec';
            max_a2_vec = (q_a2_vec == max(q_a2_vec));
            temp_pi(si, :) = max_a2_vec / sum(max_a2_vec);
        end    
        
        % T_c(s'|s,pi) = sum_a [pi(a|s) * T_c(s'|s,a)]
        T_pi = squeeze(T(c,:,:,:)) .* repmat(permute(reshape(temp_pi, N_s, N_a, 1), [3 1 2]), N_s', 1, 1); % N_s' x N_s x N_a
        T_pi = squeeze(sum(T_pi, 3)); % N_s' x N_s
        
        % stage 2 - V(s) vector [N_s x 1]
        % V(s2) = SUM_a2 {pi(s2, a2) * Q(s2, a2)} = SUM_a2 pi(s2, a2) * [SUM_s3 T_c(s2,a2,s3) {R(s3)}]
        % V(s2) = SUM_s3 { T_c(s3|s2,pi) * R(s3) }
        v_s2_vec = r_s3_vec' * T_pi; % [1 x N_s'] * [N_s' x N_s]
        v_s2_vec = v_s2_vec';
        
        % stage 1 - Q(s_fix, a) vector [N_a x 1]
        % Q(s1, a1) = SUM_s2 T_c(s1,a1,s2) {V(s2)}
        q_a1_vec = v_s2_vec' * squeeze(T(c, :, state, :)); % [1 x N_s] * [N_s' x N_a]        
        max_a1_vec = (q_a1_vec == max(q_a1_vec));
        temp_pi(1, :) = max_a1_vec / sum(max_a1_vec);
end
T_pi = squeeze(T(c,:,:,:)) .* repmat(permute(reshape(temp_pi, N_s, N_a, 1), [3 1 2]), N_s', 1, 1); % N_s' x N_s x N_a
T_pi = squeeze(sum(T_pi, 3));

% after n step transitions
T_pi_step = T_pi ^ options.TransitionStep;
d(1) = sum(T_pi_step(goal1_id, s));
d(2) = sum(T_pi_step(goal2_id, s));
d(3) = sum(T_pi_step(goal3_id, s));
d(4) = sum(T_pi_step(goal4_id, s));

% sanity check
if sum(d) - 1 > 0.0001
    disp(context)
    disp(state)
    disp(d)
    error('(goal_distr_complex) goal distribution sum-to-one violated'); 
end


