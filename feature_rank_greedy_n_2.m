function [selected_fea, score] = feature_rank_greedy_n_2(xa, y)
total_fea = size(xa,2);
no_of_data = size(xa,1);
idx_val = 0;
score = [];
% for i=1: size(xa,2)
%     [r_val(i)] = relieff_test(xa(:,i), y,xa,7,'method','classification');
% end

[test,r_val]=relieff(xa,y,7,'method','classification');

[val1, idx1] = max(r_val);
score  = [score val1];
selected_fea = [];
selected_fea = [selected_fea idx1];
idx_val = idx_val+val1;
remain = setdiff(test,selected_fea);
test(1) = [];
for jj = 1:size(remain, 2)
    %remain1 = setdiff(1:size(xa,2),selected_fea);
    r_val1 = [];
    mrmr = [];
    c=[];
    count = 1;
    
    zz=[];
    zz = [selected_fea test(jj)];
    [r_val1] = relieff_test(xa(:,zz),y,xa,7,'method','classification');
    % [r_val1(j)] = relieff_test(xa(:,zz),y,7,'method','classification');
    
    
    [val, idx] = max(r_val1);
    %     max_inter_idx = index_be(idx);
    %     val2 = (idx_val+r_val(max_inter_idx))/(size(selected_fea,2)+1);
    %     if val>val1 | val > val2
    %     if val > puloma
    
    %     if val >1
    %         print('hi');
    %     end
    if val>val1
        %pp = remain1(idx);
        selected_fea = [selected_fea test(jj)];
        score  = [score val];
        %idx_val = idx_val + r_val(max_inter_idx);
        val1 = val;
    end
end