function A=Prepare_A_matrix(k,phi,reactant_vector,product_vector,...
                          truncated_size)
                      
num_reactions=size(reactant_vector,2);
reaction_vector=product_vector-reactant_vector;

% the row index of non-zero elements
rowI=zeros(1,(num_reactions+1)*prod(truncated_size)); 
colI=zeros(size(rowI)); % the cololumn index of non-zero elements
value=zeros(size(rowI)); % the value of non-zero elements.
i=0;                     % the number of non-zero elements. 

for LI=1:prod(truncated_size) %LI is the linear index of the matrix
    X=Linear_index_2_matrix_index(truncated_size,LI)-1; %the state associated with the point.
    X=double(X);
    if size(X,1)==1
       X=X'; % make X a column vector
    end
    AD=0; % the value of diagonal index.
    
    % jst reaction
    for j=1:num_reactions
        a=propensity_function(j,X,k,phi); %propensity at this point
        AD=AD-a; % update the diagonal element
        % update the off-diagonal element
        X_off=X-reaction_vector(:,j); %the state associated with the off-diagonal element
        if min(X_off)>=0 && all(X_off+1<=truncated_size)==1 
          % the X_off is within the truncated space
          a_off=propensity_function(j,X_off,k,phi); %propensity associated with x_off;
          LI_off= Matrix_index_2_linear_index(truncated_size,X_off+1);
          i=i+1;
          rowI(i)=LI;
          colI(i)=LI_off;
          value(i)=a_off;
        end
    end
    
    % diagonal elements
    i=i+1;
    rowI(i)=LI;
    colI(i)=LI;
    value(i)=AD;
end

A=sparse(rowI(1:i),colI(1:i),value(1:i));

%save ./FSP_CME_solver/matrix_A A