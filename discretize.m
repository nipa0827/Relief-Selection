function dis = discretize(x, edges)
    len = length(x);
    dis = zeros(len, 1);
    empty_count = 0;
    count = 1;
    pp = [];
    aa = [];
    
% % without Empty cell %
%     for j=1:len
%        for i=1:length(edges)-1
%             if x(j) >= edges(i) && x(j) <= edges(i+1)
%                 dis(j) = i;
%                 break;
%             end
%         end 
%     end
% % without Empty cell %


% % with Empty cell %
       for i=1:length(edges)-1
            a = find( x >= edges(i) & x < edges(i+1) );
            
            aa = [a;pp];
                if length(unique(aa)) > empty_count
                dis(unique(aa)) = i;
%                 count = count +1;
                pp = [];
                aa = [];
                elseif length(unique(aa)) <= empty_count & i == length(edges)-1
                    dis(unique(aa)) = i-1;
                else 
               pp = aa;
%                count = count -1;
            end
        end 
% % with Empty cell %
end