function p=jump_update(p,index_leader,X,k,phi,truncated_size,...
                       follower_subsystems,leader_subsystems_reactions,...
                       reaction_decomposition_F,...
                       rowI_Jump,colI_Jump,reactionI_Jump,followerX_Jump)
       
% here mu is the reaction of the jump, p is the probability to be updated

% find the follower subsystem to be updated
index_follower=...
    reaction_decomposition_F(leader_subsystems_reactions{index_leader});

index_follower=max(index_follower);

if index_follower==0
   return; % no follower subsystem needs to be update 
end


% find the leader system accociated with mu
% for i=1:size(leader_subsystems_reactions,2)
%     if ismember(mu,leader_subsystems_reactions{i})
%         index_leader=i;
%         break;
%     end
% end

% construct the matrix
A=Construct_matrix_for_Jumps(X,k,phi,truncated_size,...
                                      follower_subsystems{index_follower},...
                                      rowI_Jump{index_leader},colI_Jump{index_leader},...
                                      reactionI_Jump{index_leader},followerX_Jump{index_leader});
                                  
% update the probability
q=A*p{index_follower};
%sumq=sum(q);
p{index_follower}=q/sum(q);