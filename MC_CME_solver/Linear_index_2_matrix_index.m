function MI=Linear_index_2_matrix_index(truncated_size,Index)

Index=Index-1;
for i=1:size(truncated_size,1)
    MI(i)=mod(Index,truncated_size(i));
    Index=Index-MI(i);
    Index=Index/truncated_size(i);
end

MI=MI+1;