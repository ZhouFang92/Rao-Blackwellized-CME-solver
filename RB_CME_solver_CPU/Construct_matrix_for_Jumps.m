function A=Construct_matrix_for_Jumps(X,k,phi,truncated_size,...
                                      species_follower_subsystem,...
                                      rowI_Jump,colI_Jump,reactionI_Jump,followerX_Jump)
                                  
% This function gives the matrix for the update at Jump events

truncated_size_follower_subsystem=truncated_size(species_follower_subsystem);
size_of_the_matrix=prod(truncated_size_follower_subsystem);

value=zeros(size(rowI_Jump)); % the value of the non-zero matrix
     
state=X.*ones(size(X,1),size(value,2));
state(species_follower_subsystem,:)=followerX_Jump(species_follower_subsystem,:);
     
value=propensity_for_RB(state,k,phi,reactionI_Jump);
     
A=sparse(rowI_Jump,colI_Jump,value,size_of_the_matrix,size_of_the_matrix);