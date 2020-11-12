function [discretized, edges] = quantize(x, k)

discretized = zeros(size(x,1),size(x,2));
 

for ii=1:size(x,2)
    minVal = min(x(:,ii));
    maxVal = max(x(:,ii));
    stepSize = (maxVal-minVal)/k;
    
%     if minVal == 0 
%         edges(:,ii) = 0;
%         continue;
%     end
%     if maxVal == 0 
%         edges(:,ii) = 0;
%         continue;
%     end
    if (minVal == 0 & maxVal == 0) | stepSize == 0 
        edges(:,ii) = zeros(6,1);
        continue;
    end
    edges(:,ii) = (minVal:stepSize:maxVal);
    edges(end,ii) = edges(end,ii)+ 1;
    discretized(:,ii) = discretize(x(:,ii), edges(:,ii));
end
end