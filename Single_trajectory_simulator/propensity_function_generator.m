function propensity_function_generator(propensity,name)

% This function generates an m file which can calculate the propensity of
% the network system given the state and the reaction index.


fid=fopen('./Single_trajectory_simulator/propensity_function.m','wt');

fprintf(fid, 'function a=propensity_function(i,X,k,phi)\n'); % The title of the function

fprintf(fid, '%% The propensity function of %s\n\n', name); % show the name of the network

for i=1:size(propensity,2)
    fprintf(fid, 'if i==%d\n', i); 
    fprintf(fid, 'a=%s;\n', propensity(i));
    fprintf(fid, 'end\n\n');
end

fclose(fid);