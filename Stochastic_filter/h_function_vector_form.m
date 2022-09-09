function a=h_function_vector_form(h,X)
% The h function of Genetic_toggle_switch in the vector form

if size(h,2)==0
a=[];
return;
end

a(1,:)=X(3,:);
