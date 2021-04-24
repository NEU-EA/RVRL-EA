function RVRLEA(Global)
% <algorithm> <R> 
% alpha ---   2 --- The parameter controlling the rate of change of penalty

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    alpha1 = Global.ParameterSet(2);
    gamma=0.80;              % discount factor
    alpha=0.50;              % learning rate
    epsilon=0.9;
    k1 = 0.5;
    k2 = 0.5;
    cumReward=0;             % counter to calculate accumulated reward
    exploitCount=0;
    exploreCount=0;

    %% Generate the reference points and random population
    [V0,Global.N] = UniformPoint(Global.N,Global.M);
    [q,q1,R]=InitMatrix(Global.N);
    Population    = Global.Initialization();
    V             = [V0;V0];%rand(Global.N,Global.M)
    V_old = V(Global.N+1:end,:);
    
    %% 实验分析用
    v1=[];
    a1=[];
    s1=[];
    q11=cell(1,1200);
    tt=1;
   

    %% Optimization
    while Global.NotTermination(Population)

       temp = V(Global.N+1:end,:);
        %% Reinforcement learning  
         % 开始的状态
         % 附近有一个非支配个体，state=1,附近无个体，state=2，附近有多个非支配个体，state=3，附近一个支配个体，state=4
         % 附近多个支配个体，state=5，附近多个混合个体，state=6
         [state,new_ind] = getState(Population.objs,V(Global.N+1:end,:));
         previous_ind    = new_ind;
         Action = [];
         for i=1:Global.N  %n个参考向量
             if sum(sum(abs(q1{i}-q{i})))<0.001 && sum(sum(q{i} >0)) && epsilon<0.0001 %&& i==Global.N %该向量的Q表收敛，值稳定
                if state(i)==4
                    x=find(R{i}(state(i),:)>=-50);         % find possible action of this state
                    if size(x,1)>0
                        [~,qmax]=(max(q{i}(state(i),x(1:end)))); % check for action with highest Q value
                        x1 = x(qmax);  % set action with highest Q value as next state
                        Action = horzcat(Action,x1);
                        % 采取动作
                        switch x1
                            case 4
                                V(Global.N+i,:) =  V(Global.N+i,:);%不变
                            case 1
                                V(Global.N+i,:) = Regeneration(V(Global.N+1:end,:),Population.objs);
                                %                                 V(Global.N+i,:) = rand(1,Global.M).*max(Population.objs,[],1);%均匀产生
                            case 3
                                V(Global.N+i,:) = V_old(i,:);%返回上次的位置
                            case 2
                                V(Global.N+i,:) = V(i,:);%返回均匀位置
                        end
                    end
                else
                  Action = horzcat(Action,1);  
                end
            else
                % 两个阶段 探索和开发     
                % select any action from this state using ?-greedy
                x=find(R{i}(state(i),:)>=-50);         % find possible action of this state
                if size(x,1)>0
                    r=rand; % get a uniform random number between 0-1
                    
                    % choose either explore or exploit
                    if r>=epsilon   % exploit
                        [~,qmax]=(max(q{i}(state(i),x(1:end)))); % check for action with highest Q value
                        x1 = x(qmax);  % set action with highest Q value as next state
                        % 采取动作
                        switch x1
                            case 4
                               V(Global.N+i,:) =  V(Global.N+i,:);%不变
                            case 1
                                V(Global.N+i,:) = Regeneration(V(Global.N+1:end,:),Population.objs);
%                                 V(Global.N+i,:) = rand(1,Global.M).*max(Population.objs,[],1);%均匀产生
                            case 3
                                V(Global.N+i,:) = V_old(i,:);%返回上次的位置
                            case 2
                                V(Global.N+i,:) = V(i,:);%返回均匀位置
                        end
                        Action = horzcat(Action,x1);  %记录每个向量采取的动作
                        
                        exploitCount=exploitCount+1;
                        display('exploit');
                        
                    else        % explore
                        x1=RandomPermutation(x);   % randomize the possible action
                        x1=x1(1);                  % select an action (only the first element of random sequence)
                        % 采取动作
                        switch x1
                            case 4
                               V(Global.N+i,:) =  V(Global.N+i,:);%不变
                            case 1
                                V(Global.N+i,:) = Regeneration(V(Global.N+1:end,:),Population.objs);
%                                 V(Global.N+i,:) = rand(1,Global.M).*max(Population.objs,[],1);%均匀产生
                            case 3
                                V(Global.N+i,:) = V_old(i,:);%返回上次的位置
                            case 2
                                V(Global.N+i,:) = V(i,:);%返回均匀位置
                        end
                        Action = horzcat(Action,x1);  %记录每个向量采取的动作
                        exploreCount=exploreCount+1;
                        display('explore');
                    end   
                end
             end
            %% 保存某条向量采取的动作和状态，实验分析用
            if i==91
              a1=[a1,x1];
              s1=[s1,state(i)];
              v1=[v1;V(Global.N+i,:)];
            end
         end
          
       %%  更新Q-table
         old_q    =   q;      
        [state2,new_ind] = getState(Population.objs,V(Global.N+1:end,:));
        epsilon=epsilon*0.99;% decrease epsilon
        for i = 1:Global.N
            if sum(sum(abs(q1{i}-q{i})))<0.001 && sum(sum(q{i} >0)) && epsilon<0.0001
                i
            else
                if i==91
                    q11{tt}=q{i};%记录Q表的变化过程
                    tt=tt+1;
                end
                R{i}(state(i),Action(i))=k1/(exp(new_ind(i)-1)-0.9)+k2*exp(previous_ind(i)-new_ind(i)); % Update reward
                %                 R{i}(state(i),Action(i))=exp(new_ind(i)-previous_ind(i)); % Update reward
                cumReward=cumReward+q{i}(state(i),Action(i)); %keep track of cumulative reward for graph

                x2 = find(R{i}(:,Action(i))>=-50);   % find possible steps from next step
                qMax=(max(q{i}(x2(1:end),Action(i)))); % extract qmax from all possible next states
                q{i}(state(i),Action(i))= q{i}(state(i),Action(i))+alpha*((R{i}(state(i),Action(i))+gamma*qMax)-q{i}(state(i),Action(i)));    % Temporal Difference Error
                %state(i)=x1;    % set state to next state
            end
        end
        q1=old_q;

        
        %%  环境选择
        MatingPool = randi(length(Population),1,Global.N);
        Offspring  = GA(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],V,(Global.gen/Global.maxgen)^alpha1);
         
        %% 最后
        if Global.evaluated >= Global.evaluation
            Population = Truncation(Population,Global.N);
            folder = fullfile('Experiment\Data3');
            [~,~]  = mkdir(folder);
            save(fullfile(folder,sprintf('%s_%s_M%d.mat',func2str(Global.algorithm),class(Global.problem),Global.M)),'a1','s1','v1','q11');
        end
        state = state2;
        V_old = temp;

    end
end