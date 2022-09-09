function h_function_vector_form_generator(h,name)

% This function generates an m file which can calculate the h function of
% the network system given the state.

fid=fopen('./Stochastic_filter/h_function_vector_form.m','wt');

fprintf(fid, 'function a=h_function_vector_form(h,X)\n'); % The title of the function

fprintf(fid, '%% The h function of %s in the vector form\n\n', name); % show the name of the network

% Case 1: no observation channel
fprintf(fid, 'if size(h,2)==0\n'); 
fprintf(fid, 'a=[];\nreturn;\n');
fprintf(fid, 'end\n\n');

%Case 2: there are observation channels
for i=1:size(h,2)
    h(i)=translate_propensities_for_multi_trajectory_simulation...
           (h(i)); % translate
    fprintf(fid, 'a(%d,:)=%s;\n', i, h(i));
end

fclose(fid);