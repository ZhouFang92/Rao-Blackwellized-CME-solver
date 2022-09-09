function index=Matrix_index_2_linear_index(truncated_size,MI)

index=0;
unit=1;
MI=double(MI);
for i=1:max(size(MI))
    index=index+(MI(i)-1)*unit;
    unit=unit*truncated_size(i);
end

index=index+1;