function d=plot_RB_follower_plus_leader(i,j,X,P,truncated_size,follower_subsystems,...
                        species_decomposition_LF)
% i is the follower, j is the leader

sample_size=size(X,2); 

d=zeros(truncated_size(i),truncated_size(j));
                    
% indentify which subsystem do they belongs to 
for index=1:size(follower_subsystems,2)
    if ismember(i,follower_subsystems{index}) ==1
       follower_index1=index;
       break;
    end
end

for sample_index=1:sample_size
    if X(j,sample_index)+1>truncated_size(j)
       continue; % it is outside the truncated space;
    end
    for LI1=1:size(P{follower_index1,1},1)  % tranverse all the probability vector
           truncated_size_follower_subsystem=truncated_size(follower_subsystems{follower_index1});
           MI_i=Linear_index_2_matrix_index(truncated_size_follower_subsystem,LI1);
           a=find(follower_subsystems{follower_index1}==i);
           index_i=MI_i(a);
           value1=P{follower_index1,sample_index}(LI1);
           d(index_i,X(j,sample_index)+1)=...
               d(index_i,X(j,sample_index)+1)+value1/sample_size;
    end
 end

