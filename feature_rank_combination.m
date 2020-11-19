function selected_fea = feature_rank_combination(xa, y)
total_fea = size(xa,2);
no_of_data = size(xa,1);
idx_val = 0;

% for i=1: size(xa,2)
%    [r_val(i)] = relieff_test(xa(:,i), y,xa,7,'method','classification');
% end

[test,r_val]=relieff_our(xa,y,7,'method','classification');
[val1, idx1] = max(r_val);
selected_fea = [];
selected_fea = [selected_fea idx1];
idx_val = idx_val+val1;
remain = setdiff(1:size(xa,2),selected_fea);
for jj = 1:size(remain, 2)
    remain1 = setdiff(1:size(xa,2),selected_fea);
    r_val1 = [];
    mrmr = [];
    c=[];
    count = 1;
    for j = 1:size(remain1, 2)
        zz=[];
        zz = [selected_fea remain1(j)];
        [r_val1(j)] = relieff_complementary(xa(:,zz),y,xa,7,'method','classification');
       % [r_val1(j)] = relieff_test(xa(:,zz),y,7,'method','classification');
        index_be(count) = remain1(j);
        count = count+1;
    end
    [val, idx] = max(r_val1);
    max_inter_idx = index_be(idx);
    val2 = (idx_val+r_val(max_inter_idx))/(size(selected_fea,2)+1);
    %     if val>val1 | val > val2
    %     if val > puloma
    
%     if val >1
%         print('hi');
%     end
    if val>val1 
        pp = remain1(idx);
        selected_fea = [selected_fea pp];
        idx_val = idx_val + r_val(max_inter_idx);
        val1 = val;
    else
        break;
    end
end