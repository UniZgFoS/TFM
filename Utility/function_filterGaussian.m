function Afiltered = function_filterGaussian( A , dt , f0 , f_sigma , f_upperCutoff , f_lowerCutoff )
    % A( t = (n-1)*dt ) = A(n), for n=1,2,...
    
    NtimePoints = length(A);
    
    X = fft(A);
    X( floor(NtimePoints/2)+1 : end ) = 0;
    
    for k = 1 : NtimePoints
        fk = (k-1)/NtimePoints * 1/dt;
        
        if( fk<f_upperCutoff && f_lowerCutoff<fk )
            X(k) = X(k) * exp( -0.5 * ((fk-f0)/f_sigma)^2 );
        else
            X(k) = 0;
        end
    end
    
    Afiltered = ifft(X);
    
    return;
end

