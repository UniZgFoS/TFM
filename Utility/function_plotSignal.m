function function_plotSignal( A , dt_MHz )
    
    A = reshape( A , [1,length(A)] );

    tk = [];
    fk = [];
    for k = 1 : length(A)
        fk = [ fk , (k-1)/length(A) * 1/dt_MHz ];
        tk = [ tk , (k+1)           *   dt_MHz ];
    end

    figure, plot( fk , abs(fft(A)) , '.-' );
    xlabel('f_k [MHz]');
    grid on; grid minor;
    xlim([0,max(fk)]);

    figure, plot( tk , real(A) , '.-' );
    xlabel('$t_n = (n-1)\delta t_s $ [$\mu s$]','Interpreter','latex','FontSize',30);
    ylabel('${\bf{A}}_{5,10}( t_n )$', 'Interpreter','latex', 'FontSize',30);
    grid on; grid minor;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',15);
    xlim([0,max(tk)]);
end

