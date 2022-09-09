function p=A_times_P_FFSP(follower_index,p0,...,
                          X,k,phi,truncated_size,...
                          reactant_vector,product_vector,...
                          follower_subsystems,...
                          reaction_decomposition_LF,...
                          reaction_decomposition_F)

% this function returns A*p0 for the FFSP method, where A is the evolution
% matrix for the corresponding follower subsystem, p0 is its distribution.


truncated_size_follower_subsystem=truncated_size(follower_subsystems{follower_index});
num_reactions=size(reaction_decomposition_F,2);
reaction_vector=product_vector-reactant_vector;

p=zeros(size(p0));

for LI=1:prod(truncated_size_follower_subsystem) % LI is the linear index
    MI=Linear_index_2_matrix_index(truncated_size_follower_subsystem,LI);
    if size(MI,1)==1
        MI=MI'; % make MI a column vector
    end
    CurrentX=X;
    CurrentX(follower_subsystems{follower_index})=MI-1; % the state for the LI-th element of the probability vector
    
    for j=1:num_reactions
        if reaction_decomposition_F(j)~=follower_index
            continue;     % this reaction does not belongs to this system.          
        end
        % the outflux
        a=propensity_function(j,CurrentX,k,phi); %propensity at this point
        p(LI)=p(LI)-a*p0(LI);
        % the outflux
        if reaction_decomposition_LF(j)==0
           continue;      % Off-diagonal element does not contain leader level reactions.  
        end
        X_off=CurrentX-reaction_vector(:,j); %the state associated with the off-diagonal element
        if min(X_off)>=0 &&...
           all(X_off(follower_subsystems{follower_index})+1<=truncated_size_follower_subsystem)==1 
           % the X_off is within the truncated space
              a_off=propensity_function(j,X_off,k,phi); %propensity associated with x_off;
              LI_off=Matrix_index_2_linear_index...
                    (truncated_size_follower_subsystem,X_off(follower_subsystems{follower_index})+1);
              p(LI)=p(LI)-a_off*p0(LI_off);
        end        
    end
    
end

