# conjugate-gradient-method

## Reference 
Conjugate Gradient Method `(ENG)`
- https://en.wikipedia.org/wiki/Conjugate_gradient_method#Example_code_in_MATLAB_/_GNU_Octave

켤레기울기법 `(KOR)`
- https://ko.wikipedia.org/wiki/%EC%BC%A4%EB%A0%88%EA%B8%B0%EC%9A%B8%EA%B8%B0%EB%B2%95

## Conjugate Gradient Method (CG)
The conjugate gradient method is an algorithm for the numerical solution of particular systems of linear equations, namely those whose matrix is [symmetric](https://en.wikipedia.org/wiki/Symmetric_matrix) and [positive-definite](https://en.wikipedia.org/wiki/Positive-definite_matrix). The conjugate gradient method is often implemented as an [iterative algorithm](https://en.wikipedia.org/wiki/Iterative_method), applicable to sparse systems that are too large to be handled by a direct implementation or other direct methods such as the Cholesky decomposition. 

## Cost function
Suppose we want to solve the [system of linear equations](https://en.wikipedia.org/wiki/System_of_linear_equations)

        (P1) A * x   = b    : matrix ver.
        
or,

        (P2) A( x )  = b    : function ver. 
        
for the vector `x`, where the known n x n matrix `A` is symmetric (i.e., A^T = A), positive-definite (i.e., x^T A x > 0 for all non-zero vectors x in R^n), and real, and `b` is known as well. We denote the unique solution of this system by `x^*`.

## The basic iteration CG for solving problem (matrix ver.)

        function [x] = conjgrad(A, b, x)
            r = b - A * x;
            p = r;
            rsold = r' * r;

            for i = 1:length(b)
                Ap = A * p;
                alpha = rsold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                rsnew = r' * r;
                if sqrt(rsnew) < 1e-10
                    break;
                end
                p = r + (rsnew / rsold) * p;
                rsold = rsnew;
            end
        end
 

## The basic iteration CG for solving problem (function ver.)

        function [x] = conjgrad(A, b, x, N)
            r = b - A ( x );
            p = r;
            rsold = r(:)' * r(:);

            for i = 1:N
                Ap = A ( p );
                alpha = rsold / (p(:)' * Ap(:));
                x = x + alpha * p;
                r = r - alpha * Ap;
                rsnew = r(:)' * r(:);
                if sqrt(rsnew) < 1e-10
                    break;
                end
                p = r + (rsnew / rsold) * p;
                rsold = rsnew;
            end
        end

