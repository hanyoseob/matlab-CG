%% REFERENCE
% https://en.wikipedia.org/wiki/Conjugate_gradient_method

%%
function [x, obj]  = CG(A,b,x,n,COST,bfig)

if (nargin < 6)
    bfig = false;
end

if (nargin < 5 || isempty(COST))
    COST.function	= @(x) (0);
    COST.equation	= [];
end

if (nargin < 4)
    n   = 1e2;
end


% r       = b - A*x;
r       = b - A(x);
p       = r;

rsold   = r(:)'*r(:);
obj     = zeros(n, 1);

for i = 1:n
    %    Ap   = A*p;
    Ap   = A(p);
    a    = rsold/(p(:)'*Ap(:));
    
    x    = x + a*p;
    r    = r - a*Ap;
    
    rsnew= r(:)'*r(:);
    
    if (sqrt(rsnew) < eps)
        break;
    end
    
    p    = r + (rsnew/rsold)*p;
    rsold= rsnew;
    
    obj(i)  = COST.function(x);
    
    if bfig
        figure(1); colormap gray;
        subplot(121); imagesc(abs(x));     	title([num2str(i) ' / ' num2str(n)]);
        subplot(122); semilogy(obj, '*-');  title(COST.equation);  xlabel('# of iteration');   ylabel('Objective'); 
                                            xlim([1, n]);   grid on; grid minor;
        drawnow();
    end
end

x   = gather(x);
obj = gather(obj);

end