function h_function_generator(h,name)

% This function generates an m file which can calculate the observation of
% the network system given the state.

fid=fopen('./Single_trajectory_simulator/h_function.m','wt');

fprintf(fid, 'function a=h_function(h,X)\n'); % The title of the function

fprintf(fid, '%% The h function of %s\n\n', name); % show the name of the network

% Case 1: no observation channel
fprintf(fid, 'if size(h,2)==0\n'); 
fprintf(fid, 'a=0;\nreturn;\n');
fprintf(fid, 'end\n\n');

%Case 2: there are observation channels
for i=1:size(h,2)
    fprintf(fid, 'a(%d,1)=%s;\n', i, h(i));
end

fclose(fid);