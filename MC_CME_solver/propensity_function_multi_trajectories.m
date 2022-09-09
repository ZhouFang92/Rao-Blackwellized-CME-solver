function a=propensity_function_multi_trajectories(X,k,phi,one_vector)
% The propensity function of Linear_metabolite_pathway for multi-trajectory simulations

a=zeros(26,size(X,2),'single','gpuArray');

a(1,:)=k(1);
a(2,:)=k(2).*X(1,:);
a(3,:)=k(3).*X(1,:);
a(4,:)=k(4);
a(5,:)=k(5).*X(2,:);
a(6,:)=k(6).*X(2,:);
a(7,:)=k(7);
a(8,:)=k(8).*X(3,:);
a(9,:)=k(9).*X(3,:);
a(10,:)=k(10);
a(11,:)=k(11).*X(4,:);
a(12,:)=k(12).*X(4,:);
a(13,:)=k(13);
a(14,:)=k(14).*X(5,:);
a(15,:)=k(15).*X(5,:);
a(16,:)=k(16);
a(17,:)=k(17).*X(6,:);
a(18,:)=k(18).*X(6,:);
a(19,:)=k(19);
a(20,:)=k(20).*X(7,:);
a(21,:)=k(21).*X(7,:);
a(22,:)=k(22);
a(23,:)=k(23).*X(8,:);
a(24,:)=k(24).*X(8,:);
a(25,:)=k(25);
a(26,:)=k(26).*X(9,:);
