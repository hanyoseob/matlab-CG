% function d = Dyt(u)
% d               = u;
% d(1:end-1,:)    = u(2:end,:);
% d(end,:)        = u(1,:);
% d(2:end,:)      = -(d(1:end-1,:) - u(1:end-1,:));
% d(1,:)          = u(end,:) - u(1,:);
% return

function d = Dyt(u)
[rows,cols] = size(u); 
d = zeros(rows,cols, 'like', u);
d(1:rows-1,:) = u(1:rows-1,:)-u(2:rows,:);
d(rows,:) = u(rows,:)-u(1,:);
return