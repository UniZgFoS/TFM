function R = function_Rfactor( A )
    % Q*A = [ R ; 0 ],
    % Q^T*Q = Q*Q^T = I,
    % R je gornje trokutasta
    assert( size(A,1) > size(A,2) );
    
    for j = 1 : size(A,2)
        
        [v,alpha] = Householder( A( j:end , j ) );
        
        A(j,j)       = alpha;
        A(j+1:end,j) = 0;
        
        for jj = j+1 : size(A,2)
            A( j:end , jj ) = A( j:end , jj ) + (transp(v)*A( j:end , jj )) * v;
        end
    end
    
    R = A( 1:size(A,2) , 1:size(A,2) );
    
    return;
end

function [v,alpha] = Householder( x )
    % Q = I + v*v^T.
    % Q^T*Q = Q*Q^T = I
    % Q*x = alpha*e1
    
    n = length(x);
    x = reshape( x , [length(x),1] );
    
    alpha = -sqrt( transp(x)*x );
    
    v = x;
    v(1) = v(1) - alpha;
    
    fac = -2 / (transp(v)*v);
    
    v = v .* sqrt(fac);
    return;
end



