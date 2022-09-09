function A=Construct_evolution_matrix_FFSP_optimized(X,k,phi,...
                                             species_follower_subsystem,...
                                             rowI,colI,reactionI,followerX)
                                       
 % This function generates the evolution matrix for the FFSP method. Here,
 % A{i} is the matrix for the i-th follower subsystem.
 
 % This algorithm has been optimized compared to the previous one.
 
 %for i=1:size(rowI,2)  % the i-th follower subsystem is considered
     value=zeros(size(rowI)); % the value of the non-zero matrix
     
     state=X.*ones(size(X,1),size(value,2));
     state(species_follower_subsystem,:)=followerX(species_follower_subsystem,:);
     
     value=propensity_for_RB(state,k,phi,reactionI);
     value=value-2*value.*(rowI==colI);
     
     A=sparse(rowI,colI,value);

 %end