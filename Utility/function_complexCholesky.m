function R = function_complexCholesky( A )
        
    assert( size(A,1) == size(A,2) );
    
    n = length(A);
    
    R = zeros(n,n);
    for j = 1 : n
        R(j,j) = sqrt( A(j,j) - sum( R(j,1:j-1).*R(j,1:j-1) ) );
        for i = j+1 : n
            R(i,j) = ( A(i,j) - sum( R(i,1:j-1).*R(j,1:j-1) ) ) / R(j,j);
        end
    end

    return;
end

