% function d = Dy(u)
% d               = u;
% d(1:end-1,:)    = u(2:end,:);
% d(end,:)        = u(1,:);
% d               = d - u;
% return

function d = Dy(u)
[rows,cols] = size(u); 
d = zeros(rows,cols, 'like', u);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);
d(1,:) = u(1,:)-u(rows,:);
return