function propensity_generator_for_RB(propensity,name)

% This function generates propensities for RB simulations

fid=fopen('./RB_CME_solver_CPU/propensity_for_RB.m','wt');

% The title of the function
fprintf(fid, 'function a=propensity_for_RB(X,k,phi,reactionI)\n'); 

fprintf(fid, '%% The propensity function of %s for RB simulations\n\n',...
        name); % show the name of the network
    
fprintf(fid, "a=zeros(1,size(X,2));\n\n");
for i=1:size(propensity,2) 
    s=translate_propensities_for_multi_trajectory_simulation(propensity(i));
    fprintf(fid, 'a=a+(%s).*(%d==reactionI);\n', s,i);
end

fclose(fid);