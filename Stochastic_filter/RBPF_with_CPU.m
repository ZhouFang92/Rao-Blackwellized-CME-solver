function pi_RBPF=RBPF_with_CPU...                 
                 (TY,Y,k,phi,FT,num_species,reactant_vector,product_vector,...
                  follower_subsystems,...
                  reaction_decomposition_LF,reaction_decomposition_F,...
                  species_decomposition_LF,...
                  IniDistribution_Mean,IniDistribution_Variance,...
                  sample_size,truncated_size,propensity,h,sigma)
              
% This algorithm solves the filtering problem using the Rao-Blackwellized
% particle filter. 

%Preparation 

[rowI,colI,reactionI,followerX]...
         =prepare_for_the_evolution_matrix(num_species,truncated_size,...
                                           reactant_vector,product_vector,...
                                           follower_subsystems,...
                                           reaction_decomposition_LF,...
                                           reaction_decomposition_F);
[rowI_Jump,colI_Jump,reactionI_Jump,followerX_Jump,...
           leader_subsystems_reactions,Leader_subsystems]=...
          prepare_for_the_jump_update(num_species,truncated_size,...
                                           reactant_vector,product_vector,...
                                           follower_subsystems,...
                                           reaction_decomposition_F,...
                                           species_decomposition_LF); 
                                       
A_update_or_not=update_the_evolution_matrix_or_not(num_species,...
                            propensity,reaction_decomposition_F,...
                            follower_subsystems,Leader_subsystems,...
                            species_decomposition_LF);
                        
observation_decomposition_F=observations_belonging_to_each_follower_subsystem...
                            (num_species,h,follower_subsystems);
followerX_correction=construct_state_for_follower_subsystems...
                              (num_species,truncated_size,...
                               follower_subsystems);


% Initialize particles
X=Initial_samples(IniDistribution_Mean,IniDistribution_Variance,...
                  sample_size);
P=cell(size(follower_subsystems,2),sample_size);
W=ones(1,sample_size)/sample_size;

% Initial marginal distribution 
% prepare the marginal distribution
for i=1:size(truncated_size,1)
    Marginal{i}=zeros(truncated_size(i),1);
    for j=1:truncated_size(i)
        x=j-1; % the value of the state
        if x>=IniDistribution_Mean(i)-IniDistribution_Variance(i)
           x=x-IniDistribution_Mean(i)+IniDistribution_Variance(i);
           x=double(x);
           Marginal{i}(j)= poisspdf(x,IniDistribution_Variance(i));
        end
    end
end

p=cell(1,size(follower_subsystems,2));
for i=1:size(follower_subsystems,2) % the i-th follower subsystem is considered.
    species_follower_subsystem=follower_subsystems{i};
    truncated_size_follower_subsystem=truncated_size(species_follower_subsystem);
    p{i}=zeros(prod(truncated_size_follower_subsystem),1);
    for LI=1:prod(truncated_size_follower_subsystem) % linear index
        MI=Linear_index_2_matrix_index(truncated_size_follower_subsystem,LI);
        value=1;
        for j=1:max(size(MI)) 
            value=value*Marginal{species_follower_subsystem(j)}(MI(j));
        end
        p{i}(LI)=value;
    end
    P(i,:)=p(1,i);
%     for j=1:sample_size
%         P{i,j}=p{i}; 
%     end
end



% Iteration
pi_RBPF=cell(3,max(size(TY)));
for i=1:max(size(TY))
    % determine dt
    if i==1
        dt=TY(1);
    else
        dt=TY(i)-TY(i-1);
    end
    
    % predication and correction
    y=Y(:,i);
    parfor j=1:sample_size
      [X_temp,p]=RB_simulator...
                (X(:,j),P(:,j),dt,...
                  k,phi,reactant_vector,product_vector,...
                  follower_subsystems,truncated_size,...
                  reaction_decomposition_F,species_decomposition_LF,...
                  rowI,colI,reactionI,followerX,...
                  rowI_Jump,colI_Jump,reactionI_Jump,followerX_Jump,...
                  leader_subsystems_reactions,Leader_subsystems,...
                  A_update_or_not);
       X(:,j)=X_temp;
       [W(j),P(:,j)]=Correction_step_RBPF(W(j),X(:,j),p,y,...
                                    follower_subsystems,...
                                    h,sigma,...
                                    observation_decomposition_F,...
                                    followerX_correction);
    end
     
    % save
    pi_RBPF{1,i}=X;
    pi_RBPF{2,i}=P;
    pi_RBPF{3,i}=W/sum(W);
    
    % resampling
    [X,P,W]=resampling_step_RB(X,P,W);
    X=sample_z(X,P,follower_subsystems,followerX_correction);
    1;
end

