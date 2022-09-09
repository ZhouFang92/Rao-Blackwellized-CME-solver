function d_MC=Monte_Carlo_CME_solver...
              (k,phi,FT,reactant_vector,product_vector,...
                 IniDistribution_Mean,IniDistribution_Variance,...
                 sample_size,truncated_size)
             
             
 Maximum_size_for_simultaneous_simulations=10^6;
 unsimulated_samples=sample_size;    
 
 %CPU solver
 if gpuDeviceCount==0 % no gpu device
      d_MC=zeros(truncated_size');
      if size(truncated_size,1)==1
         d_MC=zeros(truncated_size,1);
      end
      
       while unsimulated_samples>0
          size_for_this_iteration=min(unsimulated_samples,...
                               Maximum_size_for_simultaneous_simulations);
          unsimulated_samples=unsimulated_samples-size_for_this_iteration;

       % simulate the process
          [TX,X]=next_reaction_method_multi_trajectories_no_GPU...
                (k,phi,FT,reactant_vector,product_vector,...
                 IniDistribution_Mean,IniDistribution_Variance,...
                 size_for_this_iteration);
           d=Distribution_Monte_Carlo(X,truncated_size);
           d_MC=d_MC+d*size_for_this_iteration/sample_size;
        end
     return;
 end
             
 % use GPU for the simulation
   
 d_MC=zeros(truncated_size','single','gpuArray');
 if size(truncated_size,1)==1
    d_MC=zeros(truncated_size,1,'single','gpuArray');
 end
 
 while unsimulated_samples>0
       size_for_this_iteration=min(unsimulated_samples,...
                               Maximum_size_for_simultaneous_simulations);
       unsimulated_samples=unsimulated_samples-size_for_this_iteration;

       % simulate the process
       [TX,X]=next_reaction_method_multi_trajectories...
                (k,phi,FT,reactant_vector,product_vector,...
                 IniDistribution_Mean,IniDistribution_Variance,...
                 size_for_this_iteration);
        d=Distribution_Monte_Carlo(X,truncated_size);
        d=single(gpuArray(d));
        d_MC=d_MC+d*size_for_this_iteration/sample_size;
 end
 
 d_MC=gather(d_MC);
 
 
 