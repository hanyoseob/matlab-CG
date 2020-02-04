% function d = Dxt(u)
% d               = u;
% d(:,1:end-1)    = u(:,2:end);
% d(:,end)        = u(:,1);
% d(:,2:end)      = -(d(:,1:end-1) - u(:,1:end-1));
% d(:,1)          = u(:,end) - u(:,1);
% return

function d = Dxt(u)
[rows,cols] = size(u); 
d = zeros(rows,cols, 'like', u);
d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);
d(:,cols) = u(:,cols)-u(:,1);
return