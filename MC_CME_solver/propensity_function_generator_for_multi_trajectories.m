function propensity_function_generator_for_multi_trajectories(propensity,name)

% This function generates propensities for multi-trajectory simulations

fid=fopen('./MC_CME_solver/propensity_function_multi_trajectories.m','wt');

% The title of the function
fprintf(fid, 'function a=propensity_function_multi_trajectories(X,k,phi,one_vector)\n'); 

fprintf(fid, '%% The propensity function of %s for multi-trajectory simulations\n\n',...
        name); % show the name of the network

%fprintf(fid, "one_vector=ones(1,size(X,2),'single','gpuArray');\n\n"); % a row vector.
    
fprintf(fid, "a=zeros(%d,size(X,2),'single','gpuArray');\n\n",size(propensity,2));
for i=1:size(propensity,2) 
    s=translate_propensities_for_multi_trajectory_simulation(propensity(i));
    fprintf(fid, 'a(%d,:)=%s;\n',i, s);
end

fclose(fid);