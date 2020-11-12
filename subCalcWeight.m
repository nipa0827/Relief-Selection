function [avg1,K,zz] = subCalcWeight(arr1, totalData,miss, z,y)
k=7;
K=k;
Rvalue1 = zeros(size(arr1,1),1);

miass = [];
count_index = [];
hit_miss = [];
counter = 1;
[IDX,D] = knnsearch(arr1, arr1,'K', k,'Distance','cityblock');
[IDX1,D1] = knnsearch(totalData, arr1,'K', k,'Distance','cityblock');
[IDX_miss,D_miss]=knnsearch(miss, arr1,'K', k,'Distance','cityblock');
len_idx = size(arr1,1);
col_idx = size(IDX,2);
difference = max(arr1)-min(arr1);
weight = 0;
Rvalue1 = zeros(size(arr1,1),1);
hit(1:len_idx,1:col_idx) = -1;
Max = max(D1(:));
for i = 1:len_idx
    temp_weight = 0;
    for j = 1:col_idx
        %weight = weight - D(i,j)+D_miss(i,j); %%weight using distance
        temp_weight = temp_weight -D(i,j)+D_miss(i,j);
        if z(i)==y(IDX1(i,j))
            Rvalue1(i) = Rvalue1(i) + 1;
            if D(i,j) == 0
                D(i,j) = 1;
            end
            hit(i,j) = D(i,j);
        end
    end
    weight = weight+(temp_weight/7);
end
for i = 1:len_idx
    hit_max(i) = max(hit(i,:));
end
hit_max = hit_max';

%weight = sum(weight,2)/7;
Rvalue1 = weight;
avg1 = mean(Rvalue1);
zz = Rvalue1;