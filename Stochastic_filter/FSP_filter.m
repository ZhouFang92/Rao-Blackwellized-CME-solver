function pi_FSP=FSP_filter(TY,Y,k,phi,reactant_vector,product_vector,...
                                IniDistribution_Mean,...
                                IniDistribution_Variance,...
                                truncated_size,....
                                h,sigma)
                            
% this algorithm solves the filtering problem using the FSP method. 

% Preparation 
d=Initial_distribution(IniDistribution_Mean,...
                                IniDistribution_Variance,...
                                truncated_size);
d=reshape(d,[],1); % transform it into a column vector;
A=Prepare_A_matrix(k,phi,reactant_vector,product_vector,...
                          truncated_size);                      
if gpuDeviceCount>=1  % use gpu or not
       A=gpuArray(A);
end
pi_FSP=cell(1,max(size(TY))); 
                      
% Iteration 
for i=1:max(size(TY))
    % determine dt
    if i==1
        dt=TY(1);
    else
        dt=TY(i)-TY(i-1);
    end
    
    if gpuDeviceCount>=1  % use gpu or not
       d=gpuArray(d);
    end
    
    % prediction step
    d=FSP_Power_series_with_given_A(dt,d,A,truncated_size);
    
    % Correction step
    Likelihood=likelihood_matrix(Y(:,i),h,sigma,truncated_size);
    d=d.*Likelihood;
    d=d/sum(d,'all');
    
    % store the result
    pi_FSP{i}=d;
    
    if size(truncated_size,1)>1
       d=reshape(d,truncated_size'); % reshape d for the next run
    end
    
end
                            