function [rowI,colI,reactionI,followerX]...
         =prepare_for_the_evolution_matrix(num_species,truncated_size,...
                                           reactant_vector,product_vector,...
                                           follower_subsystems,...
                                           reaction_decomposition_LF,...
                                           reaction_decomposition_F)

% This function shows the column and row indices of the non-zero elements,
% the accociated reaction indices, and the state of the follower system.

 num_reactions=size(reaction_decomposition_F,2);
 reaction_vector=product_vector-reactant_vector;
 rowI=cell(1,size(follower_subsystems,2));
 colI=cell(1,size(follower_subsystems,2));
 reactionI=cell(1,size(follower_subsystems,2));
 followerX=cell(1,size(follower_subsystems,2));
 
 for i=1:size(follower_subsystems,2) % compute the matrix for the i-th follower subsystem
     truncated_size_follower_subsystem=truncated_size(follower_subsystems{i});
     size_of_the_matrix=prod(truncated_size_follower_subsystem);
     
     rowI{i}=zeros(1,2*(num_reactions)*size_of_the_matrix); 
     colI{i}=zeros(size(rowI{i})); % the cololumn index of non-zero elements
     reactionI{i}=zeros(size(rowI{i}));  % the index of the reaction
     followerX{i}=zeros(num_species,size(rowI{i},2));
     index=0;
     
     for LI=1:size_of_the_matrix % LI is the linear index
         % the state of the follower system
         MI=Linear_index_2_matrix_index(truncated_size_follower_subsystem,LI);
         if size(MI,1)==1
            MI=MI'; % make MI a column vector
         end
         CurrentX=zeros(num_species,1);
         CurrentX(follower_subsystems{i})=MI-1;  
         
         for j=1:num_reactions
             if reaction_decomposition_F(j)~=i
               continue;     % this reaction does not belongs to this system.          
             end
             % update the diagonal element
             index=index+1;
             rowI{i}(index)=LI;
             colI{i}(index)=LI;
             reactionI{i}(index)=j;
             followerX{i}(:,index)=CurrentX;
             % update the off-diagonal element
             if reaction_decomposition_LF(j)==0
                continue;      % Off-diagonal element does not contain leader level reactions.  
             end
             X_off=zeros(num_species,1);
             X_off(follower_subsystems{i})...
                 =CurrentX(follower_subsystems{i})-reaction_vector(follower_subsystems{i},j); %the state associated with the off-diagonal element
             if min(X_off)>=0 &&...
                all(X_off(follower_subsystems{i})+1<=truncated_size_follower_subsystem)==1 
                % the X_off is within the truncated space
                LI_off=Matrix_index_2_linear_index...
                    (truncated_size_follower_subsystem,X_off(follower_subsystems{i})+1);
                index=index+1;
                rowI{i}(index)=LI;
                colI{i}(index)=LI_off;
                reactionI{i}(index)=j;
                followerX{i}(:,index)=X_off;
              end 
         end 
         % the LI-th row has been updated
     end
     % the matrix for the i-th follower subsystem has been updated 
     rowI{i}=rowI{i}(1:index); 
     colI{i}=colI{i}(1:index); 
     reactionI{i}=reactionI{i}(1:index); 
     followerX{i}=followerX{i}(:,1:index);
 end