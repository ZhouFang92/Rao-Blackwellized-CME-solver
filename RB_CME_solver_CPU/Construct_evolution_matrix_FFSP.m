function A=Construct_evolution_matrix_FFSP(X,k,phi,truncated_size,...
                                           reactant_vector,product_vector,...
                                           follower_subsystems,...
                                           reaction_decomposition_LF,...
                                           reaction_decomposition_F)
                                       
 % This function generates the evolution matrix for the FFSP method. Here,
 % A{i} is the matrix for the i-th follower subsystem. 
 
 num_reactions=size(reaction_decomposition_F,2);
 reaction_vector=product_vector-reactant_vector;
 
 for i=1:size(follower_subsystems,2) % compute the matrix for the i-th follower subsystem
     truncated_size_follower_subsystem=truncated_size(follower_subsystems{i});
     size_of_the_matrix=prod(truncated_size_follower_subsystem);
     %A{i}=zeros(size_of_the_matrix);
     
     rowI=zeros(1,(num_reactions+1)*prod(truncated_size_follower_subsystem)); 
     colI=zeros(size(rowI)); % the cololumn index of non-zero elements
     value=zeros(size(rowI)); % the value of non-zero elements.
     index=0;                     % the number of non-zero elements. 
     
     % compute elements
     for LI=1:size_of_the_matrix % LI is the linear index
         MI=Linear_index_2_matrix_index(truncated_size_follower_subsystem,LI);
         if size(MI,1)==1
            MI=MI'; % make MI a column vector
         end
         CurrentX=X;
         CurrentX(follower_subsystems{i})=MI-1; % the state for the LI-th element of the probability vector
         AD=0; % the value of diagonal index.
         
         % jst reaction
         for j=1:num_reactions
            if reaction_decomposition_F(j)~=i
               continue;     % this reaction does not belongs to this system.          
            end
            a=propensity_function(j,CurrentX,k,phi); %propensity at this point
            AD=AD-a; % update the diagonal element
            % update the off-diagonal element
            if reaction_decomposition_LF(j)==0
               continue;      % Off-diagonal element does not contain leader level reactions.  
            end
            X_off=CurrentX-reaction_vector(:,j); %the state associated with the off-diagonal element
            if min(X_off)>=0 &&...
               all(X_off(follower_subsystems{i})+1<=truncated_size_follower_subsystem)==1 
                % the X_off is within the truncated space
                a_off=propensity_function(j,X_off,k,phi); %propensity associated with x_off;
                LI_off=Matrix_index_2_linear_index...
                    (truncated_size_follower_subsystem,X_off(follower_subsystems{i})+1);
                index=index+1;
                rowI(index)=LI;
                colI(index)=LI_off;
                value(index)=a_off;
                %A{i}(LI,LI_off)=a_off;
            end   
         end
    
         % diagonal elements
         index=index+1;
         rowI(index)=LI;
         colI(index)=LI;
         value(index)=AD;
         %A{i}(LI,LI)=AD; 
     end
     
     % Finalize the matrix for the i-th follower subsystem
     A{i}=sparse(rowI(1:index),colI(1:index),value(1:index));
 end