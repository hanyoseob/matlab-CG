% % function d = Dx(u)
% % d               = u;
% % d(:,1:end-1)    = u(:,2:end);
% % d(:,end)        = u(:,1);
% % d               = d - u;
% % return

function d = Dx(u)
[rows,cols] = size(u); 
d = zeros(rows,cols, 'like', u);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
d(:,1) = u(:,1)-u(:,cols);
return