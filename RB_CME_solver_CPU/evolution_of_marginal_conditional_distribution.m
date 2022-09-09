function [p,A0]=evolution_of_marginal_conditional_distribution(A0,p,FT,X,k,phi,...
           follower_subsystems,rowI,colI,reactionI,followerX,A_evolution_update_check)


% this function compute the evolution of the marginal conditional
% distribution for a time period of T

epsilon=10^-2;%10^-1.5; % tolerance of the error 


for i=1:size(follower_subsystems,2) % the i-th follower subsystem is considered
    
    % prepare the matrix
    if A_evolution_update_check(i)==1 % the matrix needs to be updated
       A=Construct_evolution_matrix_FFSP_optimized(X,k,phi,...
                                             follower_subsystems{i},...
                                             rowI{i},colI{i},reactionI{i},followerX{i});
       A0{i}=A;                                  
    else
       A=A0{i};
    end

    
    if size(A,1)<=30
       A=full(A); 
    end
                                         
    % compute the result
    q=p{i};
    t=0;
    dt=2/norm(A,1);    % the step size is a trick
    while t<FT
      dt=min(dt,FT-t);
      n=1;
      add=A*q*dt/n; % compute the coefficient for the first order term
      q=q+add;
      
%      for n=2:5        % expan to the xxx-th order
%          add=A*add*dt/n; 
%          q=q+add;
%      end
      
      error=norm(add);
      while error>epsilon 
           n=n+1;
           add=A*add*dt/n;
           q=q+add;
           error=norm(add);
      end
    
      q=max(q,0);
      
      t=t+dt;
    end
    
    
    % update the result
    p{i}=q;
    
    
end