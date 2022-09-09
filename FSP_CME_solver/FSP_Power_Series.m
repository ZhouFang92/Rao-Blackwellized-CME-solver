function d_FSP=FSP_Power_Series(k,phi,reactant_vector,product_vector,...
                                FT,...
                                IniDistribution_Mean,...
                                IniDistribution_Variance,...
                                truncated_size)
                            
% This function solves the CME using the FSP method. Specifically, the
% matrix expoentials are approximated by the Taylor expansion. 

epsilon=10^-15; % tolerance of the error 

%Initial distribution
d0=Initial_distribution(IniDistribution_Mean,...
                                IniDistribution_Variance,...
                                truncated_size);
d0=reshape(d0,[],1); % transform it into a column vector;
A=Prepare_A_matrix(k,phi,reactant_vector,product_vector,...
                          truncated_size);
%load ./FSP_CME_solver./matrix_A;

                      
% Check if we have GPU or not 
if gpuDeviceCount>=1
   d0=gpuArray(d0);
   A=gpuArray(A);
end
 
% Solve the CME
d_FSP=d0;
t=0;
dt=2/norm(A,1);    % the step size is a trick
dt=gather(dt);
%tic
while t<FT
      dt=min(dt,FT-t);
      i=1;
      add=A*d_FSP*dt/i; % compute the coefficient for the first order term
      d_FSP=d_FSP+add;
      
      for i=2:5        % expan to the xxx-th order
          add=A*add*dt/i; 
          d_FSP=d_FSP+add;
      end
      
      error=norm(add);
      while gather(error)>epsilon 
           i=i+1;
           add=A*add*dt/i;
           d_FSP=d_FSP+add;
           error=norm(add);
       end
      
       d_FSP=max(d_FSP,0);
    
       t=t+dt;
end
%toc;

% finalization 
if gpuDeviceCount>=1
   d_FSP=gather(d_FSP);
end

if size(truncated_size,1)>1
   d_FSP=reshape(d_FSP,truncated_size');
end
                            