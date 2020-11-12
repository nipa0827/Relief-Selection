function [selected_fea, score] = feature_rank(xa, y,k)
total_fea = size(xa,2);
no_of_data = size(xa,1);
idx_val = 0;
score = [];
for i=1: size(xa,2)
    [r_val(i)] = relieff_test(xa(:,i), y,xa,k,'method','classification');
end

%[test,r_val]=relieff(xa,y,7,'method','classification');
% 
% corr = corrcoef(xa);
% corr(corr==1)=0;
% reduncancy = sum(corr,1)/(size(corr,2)-1);
% for i=1: size(xa,2)
%     if reduncancy(i)>0.5
%         r_val(i) = r_val(i)-reduncancy(i);
%     end
% end

[r_val, test] = sort(r_val,'descend');

[val1, idx1] = max(r_val);
score  = [score val1];
selected_fea = [];
selected_fea = [selected_fea test(1)];
idx_val = idx_val+val1;
remain1 = setdiff(1:size(xa,2),selected_fea);
for jj = 2:size(xa,2)
    zz = [selected_fea test(jj)];
    [r_val1] = relieff_test(xa(:,zz),y,xa,k,'method','classification');
    % [r_val1(j)] = relieff_test(xa(:,zz),y,7,'method','classification');
    
    
    [val, idx] = max(r_val1);
    
    if val>val1
        %pp = remain1(idx);
        selected_fea = [selected_fea test(jj)];
        score  = [score val];
        val1 = val;
        %     else
        %         break;
    end
end