function  [rowI_Jump,colI_Jump,reactionI_Jump,followerX_Jump,...
           leader_subsystems_reactions,Leader_subsystems]=...
          prepare_for_the_jump_update(num_species,truncated_size,...
                                           reactant_vector,product_vector,...
                                           follower_subsystems,...
                                           reaction_decomposition_F,...
                                           species_decomposition_LF)
                                       
% this function shows the leader_susbsystems, the reactions contained in
% eahc of them, and the state of the follower subsystem for each site of
% the probability vector;

num_reactions=size(reactant_vector,2);
reaction_vector=product_vector-reactant_vector;

[Leader_subsystems,leader_subsystems_reactions]=determine_leader_level_reactions(reactant_vector,...
                                product_vector,...
                                species_decomposition_LF);
                            
 rowI_Jump=cell(1,size(Leader_subsystems,2));
 colI_Jump=cell(size(rowI_Jump));
 reactionI_Jump=cell(size(rowI_Jump));
 followerX_Jump=cell(size(rowI_Jump));
 
for i=1:size(Leader_subsystems,2) % the i-th leader subsystem is considered
    
    index_follower=reaction_decomposition_F(leader_subsystems_reactions{i});
    
    index_follower=max(index_follower);
    
    if index_follower==0
       continue; % this leader system does not belongs to any follower subsystem 
    end
    
    truncated_size_follower_subsystem=...
                    truncated_size(follower_subsystems{index_follower});
    size_of_the_matrix=prod(truncated_size_follower_subsystem);
    
    rowI_Jump{i}=zeros(1,2*(num_reactions)*size_of_the_matrix); 
    colI_Jump{i}=zeros(size(rowI_Jump{i})); % the cololumn index of non-zero elements
    reactionI_Jump{i}=zeros(size(rowI_Jump{i}));  % the index of the reaction
    followerX_Jump{i}=zeros(num_species,size(rowI_Jump{i},2));
    index=0;
     
    for LI=1:size_of_the_matrix
         MI=Linear_index_2_matrix_index(truncated_size_follower_subsystem,LI);
         if size(MI,1)==1
            MI=MI'; % make MI a column vector
         end
         CurrentX=zeros(num_species,1);
         CurrentX(follower_subsystems{index_follower})=MI-1; % the state correspond to the LI-th row
         
         for reaction_index=1:max(size(leader_subsystems_reactions{i}))
             j=leader_subsystems_reactions{i}(reaction_index);
             % update the off-diagonal element
             %the state associated with the off-diagonal element
             X_off=zeros(num_species,1);
             X_off(follower_subsystems{index_follower})...
                 =CurrentX(follower_subsystems{index_follower})-reaction_vector(follower_subsystems{index_follower},j); 
             if min(X_off)>=0 &&...
                all(X_off(follower_subsystems{index_follower})+1<=truncated_size_follower_subsystem)==1 
                % the X_off is within the truncated space 
                LI_off=Matrix_index_2_linear_index...
                    (truncated_size_follower_subsystem,X_off(follower_subsystems{index_follower})+1);
                index=index+1;
                rowI_Jump{i}(index)=LI;
                colI_Jump{i}(index)=LI_off;
                reactionI_Jump{i}(index)=j;
                followerX_Jump{i}(:,index)=X_off;
             end
             
         end
         % the LI-th row has been updated
    end
    % the matrix for the i-th leader subsystem has been updated 
     rowI_Jump{i}=rowI_Jump{i}(1:index); 
     colI_Jump{i}=colI_Jump{i}(1:index); 
     reactionI_Jump{i}=reactionI_Jump{i}(1:index); 
     followerX_Jump{i}=followerX_Jump{i}(:,1:index);
end