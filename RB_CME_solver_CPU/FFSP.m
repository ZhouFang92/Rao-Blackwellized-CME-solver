function p=FFSP(t,X,p,k,phi,...
                follower_subsystems,truncated_size,...
                reaction_decomposition_F,species_decomposition_LF,...
                rowI,colI,reactionI,followerX,...
                rowI_Jump,colI_Jump,reactionI_Jump,followerX_Jump,...
                leader_subsystems_reactions,Leader_subsystems,...
                A_update_or_not)

% the FFSP method, (t,X) is the process, p is the initial distribution


% initialization 
CurrentX=X(:,1);
t_last=t(1);
leader_system_species=find(species_decomposition_LF==0);
index_leader_last=0; % the last leader level-reaction

% Construct the matrix 
% A=cell(1,size(follower_subsystems,2));

for i=1:size(follower_subsystems,2)
    A{i}=Construct_evolution_matrix_FFSP_optimized(CurrentX,k,phi,...
                                             follower_subsystems{i},...
                                             rowI{i},colI{i},reactionI{i},followerX{i});
end
 

for i=2:size(t,2)
    
    %find the leader system 
    dX=X(leader_system_species,i)-X(leader_system_species,i-1);
    if all(dX==0)==1 %it is not a leader level reaction
       continue;   
    end
    for j=1:size(Leader_subsystems,2)
        if all(Leader_subsystems(:,j)==dX)==1
            index_leader=j;
            break; 
        end
    end  
    % identify which matrix (of the subsystem) needs to be updated. 
    if index_leader_last>0
       A_evolution_update_check=A_update_or_not(:,index_leader_last);
    else
       A_evolution_update_check=ones(size(follower_subsystems,2),1);
    end

    % evolve p to the next time point 
     dt=t(i)-t_last;
     [p,A]=evolution_of_marginal_conditional_distribution(A,p,dt,CurrentX,k,phi,...
             follower_subsystems,rowI,colI,reactionI,followerX,A_evolution_update_check);
       
    % update according to the jump event 
     p=jump_update(p,index_leader,CurrentX,k,phi,truncated_size,...
                        follower_subsystems,leader_subsystems_reactions,...
                        reaction_decomposition_F,...
                        rowI_Jump,colI_Jump,reactionI_Jump,followerX_Jump);
    
    % update the time and state
    t_last=t(i);
    CurrentX=X(:,i);
    index_leader_last=index_leader;
end

% evolve to the final time
dt=t(end)-t_last;
if dt>0
   if index_leader_last>0
       A_evolution_update_check=A_update_or_not(:,index_leader_last);
    else
       A_evolution_update_check=ones(size(follower_subsystems,2),1);
    end
    [p,A]=evolution_of_marginal_conditional_distribution(A,p,dt,CurrentX,k,phi,...
             follower_subsystems,rowI,colI,reactionI,followerX,A_evolution_update_check);
end

% normalization
for i=1:max(size(p))
    p{i}=p{i}/sum(p{i}); 
end

    