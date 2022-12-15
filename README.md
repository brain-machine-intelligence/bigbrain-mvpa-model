# BML_BigBrain
### Building multi-modal datasets to predict individual differences in  decision-making and psychiatric disorders with behavioral modeling and AI technology development
 - This work has supported by the National Research Foundation of Korea(NRF) grant funded by the Korea government(MSIT)(No. 2021M3E5D2A0102249311).

## Task
### 2-stage Markov decision task
![image](https://user-images.githubusercontent.com/31230723/207741819-dc551411-f1b9-4d62-b1ec-241a2b3e94ae.png)
(A) Sequential two-choice Markov decision task [2]. We used 20 participants’ data. A participant starts from the initial state (S1) and makes two sequential choices (A1, A2) to get a reward in the goal state (S3). The state-action-state transitions are made with the probability p. Each fractal image indicates a distinct state ID. (B) There are two types of contexts: goal (specific or flexible goal condition) and uncertainty (high or low uncertainty level associated with different state-action-state probability (p) values). In the specific goal condition, subjects receive a reward only if the collected coin color matches the color of the coin box, whereas all colors of coins are redeemable in the flexible goal condition.

## Simulation & model fitting
* main functions: BATCH_ORI_FULL.m & ArbMat_Ori.m

## Multi-voxel pattern analysis

### SPM preprocessing (recommended)
 - Slice-time correction
 - Motion correction
 - (Coregistration)
 - Spatial normalization

### Multi-voxel pattern (ROI)
* main function: arbMBMF_boldpat.m

    - To read the neuroimaging files and save/load them as a multi-dimensional array
    - requirements
        - nii file directory path
        - SBJ_structure.mat (recorded timepoint)
    - input
        - experiment name
        - ROI name
        - participant ID
        - varargin
            - z-scoring / percent signal
            - preprocessing prefix
    - output
        - bold activity array (N_voxel x N_event x N_trial)

### Shattering dimensionality analysis [3, 4]. 
 To assess the amount of information associated with the task variables in multi-voxel patterns of brain regions, we computed the shattering dimensionality (SD) by averaging test accuracies of all linear support vector machines (SVM) trained to classify goal and uncertainty conditions. Accordingly, the SD quantiﬁes the separability of neural embeddings associated with each task condition. Since it is based on binary classiﬁers, the chance level is 0.5.

* main function: arbMBMF_shattering.m

    - 5 steps
        1. variable setting 
        2. label loading 
            - arbMBMF_load_var.m
        3. BOLD pattern loading 
            - arbMBMF_boldpat.m
        4. linear shattering 
            - linear_shattering.m
        5. result saving
    - requirements: 
        - parameter setting in the script
            - experiment name
            - ROI
            - label
            - result save path
    - input
        - participant ID
    - output
        - run_info, save_info (save them to the designated path)

### Cross-condition generalization analysis [4]
 To investigate the eﬀect of one task variable (A) on another variable(B)'s representation, we performed another SD analysis, in which SVMs were trained in one condition of variable 'A' and tested in the other conditions of it to decode variable 'B'. The average test accuracy is called the CCGP score.

## References
[1] O'Doherty, J. P., Lee, S. W., & McNamee, D. The structure of reinforcement-learning mechanisms in the human brain, Current Opinion in Behavioral Sciences 1, 94-100 (2015).

[2] Lee, S. W., Shimojo, S., & O’Doherty, J. P. Neural computations underlying arbitration between model-based and model-free learning. Neuron, 81(3), 687-699 (2014).

[3] Rigotti, M., Barak, O., Warden, M. R., Wang, X. J., Daw, N. D., Miller, E. K., & Fusi, S. The importance of mixed selectivity in complex cognitive tasks. Nature, 497(7451), 585-590 (2013).

[4] Bernardi, S., Benna, M. K., Rigotti, M., Munuera, J., Fusi, S., & Salzman, C. D. The geometry of abstraction in the hippocampus and prefrontal cortex. Cell, 183(4), 954-967 (2020).


## Data collection plan

* 2021년 12월 fMRI 촬영 계획

KBSI 오창센터 3T Human MRI System (https://use.kbsi.re.kr/html/user/equip/equip_view.php?e_code=AR53&cPage=1&f_name=&sel_major=&f_equ=MRI&f_branch=%EC%98%A4%EC%B0%BD&f_dept=&page=special)

![fMRI_schedule_2021Dec](https://user-images.githubusercontent.com/31230723/145715887-18b8b29b-8956-4d7c-876b-25687e4dd519.png)
