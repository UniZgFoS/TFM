clear;
clc;
close all;
addpath('../Utility/');

load('../FMCdatabase/Bristol''s_FMC_with_64_elements_5MHz_probe.mat');

ImgParams = struct();

ImgParams.ImageLength_mm  = 45; %reconstruction zone length along x-axis
ImgParams.ImageDepth_mm   = 60; %reconstruction zone length along z-axis

ImgParams.N_Xblocks       = 10; %number of pixel blocks along x-axis
ImgParams.N_Zblocks       = 10; %number of pixel blocks along z-axis

ImgParams.N_pixPerXblock  = 20; %pixels per X block
ImgParams.N_pixPerZblock  = 20; %pixels per Z block

ImgParams.NsingvalsRetain = 50; %approximation rank




%% Calculating frequencies
fmin_MHz = FMC.Probe.LowerCutoffFrequency_MHz;
fmax_MHz = FMC.Probe.UpperCutoffFrequency_MHz;
fk = [];
for k = 1 : FMC.NtimePoints
    f_MHz = (k-1)*FMC.SamplingFreqency_MHz/FMC.NtimePoints;
    if( fmin_MHz<=f_MHz && f_MHz<=fmax_MHz )
        fk = [ fk , f_MHz ];
    end
end



%% Calculating TransducerPosition
TransducerPosition = zeros( 2 , FMC.Probe.Nelements );
for iel = 1 : FMC.Probe.Nelements
     L = ImgParams.ImageLength_mm;
     p = FMC.Probe.ElementPitch_mm;
     w = FMC.Probe.ElementWidth_mm;
     N = FMC.Probe.Nelements;

     Probelen = p*(N-1) + w;

     xpos = (L-Probelen+w)/2 + (iel-1)*p;
     zpos = 0;

     TransducerPosition(:,iel) = [ xpos ; zpos ];
end





%% Calculating U and V
U = zeros( FMC.Probe.Nelements*length(fk) , ImgParams.NsingvalsRetain                          , ImgParams.N_Xblocks*ImgParams.N_Zblocks , 'single' );
V = zeros( ImgParams.NsingvalsRetain      , ImgParams.N_pixPerXblock*ImgParams.N_pixPerZblock  , ImgParams.N_Xblocks*ImgParams.N_Zblocks , 'single' );

i_block = 0;
for ixblock = 1 : ImgParams.N_Xblocks
    for izblock = 1 : ImgParams.N_Zblocks
        
        i_block = i_block + 1;
        fprintf('Calculating low-rank approximation for pixel block: ( %3d/%3d , %3d/%3d )\n', ixblock, ImgParams.N_Xblocks, izblock, ImgParams.N_Zblocks);
        
        TimeDelays = zeros( FMC.Probe.Nelements , ImgParams.N_pixPerXblock , ImgParams.N_pixPerZblock );
        for ix = 1 : ImgParams.N_pixPerXblock
            for iz = 1 : ImgParams.N_pixPerZblock
                xpix = ( (ixblock-1)*ImgParams.N_pixPerXblock + ix-0.5 ) * ImgParams.ImageLength_mm / (ImgParams.N_pixPerXblock*ImgParams.N_Xblocks);
                zpix = ( (izblock-1)*ImgParams.N_pixPerZblock + iz-0.5 ) * ImgParams.ImageDepth_mm  / (ImgParams.N_pixPerZblock*ImgParams.N_Zblocks);

                for iel = 1 : FMC.Probe.Nelements
                    c = FMC.SpecimenUltrasonicSpeed_mmperus;
                    TimeDelays( iel , ix , iz ) = norm( TransducerPosition(:,iel) - [xpix;zpix] ) / c;
                end

            end
        end
        

        
        E = zeros( FMC.Probe.Nelements*length(fk) , ImgParams.N_pixPerXblock*ImgParams.N_pixPerZblock );
        jj = 0;
        for ix = 1 : ImgParams.N_pixPerXblock
            for iz = 1 : ImgParams.N_pixPerZblock
                jj = jj + 1;
                
                xpix = ( (ixblock-1)*ImgParams.N_pixPerXblock + ix-0.5 ) * ImgParams.ImageLength_mm / (ImgParams.N_pixPerXblock*ImgParams.N_Xblocks);
                zpix = ( (izblock-1)*ImgParams.N_pixPerZblock + iz-0.5 ) * ImgParams.ImageDepth_mm  / (ImgParams.N_pixPerZblock*ImgParams.N_Zblocks);
                
                
                ii = 0;
                for f = fk
                    for iel = 1 : FMC.Probe.Nelements
                        ii = ii + 1;
                        
                        L = ImgParams.ImageLength_mm;
                        p = FMC.Probe.ElementPitch_mm;
                        w = FMC.Probe.ElementWidth_mm;
                        N = FMC.Probe.Nelements;
                        Probelen = p*(N-1)+w;

                        xpos = (L-Probelen+w)/2 + (iel-1)*p;
                        zpos = 0;
                        
                        a_mm      = FMC.Probe.ElementWidth_mm;
                        c_mmperus = FMC.SpecimenUltrasonicSpeed_mmperus;
                        lambda_mm = c_mmperus/f;
                        sintheta  = abs(xpos-xpix) / sqrt( (xpos-xpix)^2 + (zpos-zpix)^2 );
                        
                        xx = pi * a_mm * sintheta / lambda_mm;
                        apodization = (sin(xx)/xx)^2;
                        apodization = 1.0;
                        
                        E(ii,jj) = apodization * exp( 2*pi*1j * f * TimeDelays( iel , ix , iz ) ) / sqrt(FMC.NtimePoints);
                    end
                end
                
            end
        end
        
        
        
        [u,s,v] = svd(E);
        rank = ImgParams.NsingvalsRetain;
        U(:,:,i_block) = single( u(:,1:rank) * sparse(s(1:rank,1:rank)) );
        V(:,:,i_block) = single( v(:,1:rank)' );

        
        
%         % plotting singular values
%         singvals = diag(s);
%         semilogy( singvals(1:100)/max(singvals) , 'x' );
%         ylim([0,1]);
%         line([0,100],[0.10,0.10]);
%         line([0,100],[0.01,0.01]);
%         xlabel('$i$','Interpreter','latex','FontSize',20);
%         ylabel('$\sigma_i/\sigma_{\mathrm{max}}$','Interpreter','latex','FontSize',20);
%         grid on; grid minor;
%         pause(0.1);
%         
    end
end



clc;
close all;





%% Preparation phase
kfmin = 1 + FMC.Probe.LowerCutoffFrequency_MHz / FMC.SamplingFreqency_MHz * FMC.NtimePoints;
kfmax = 1 + FMC.Probe.UpperCutoffFrequency_MHz / FMC.SamplingFreqency_MHz * FMC.NtimePoints;


Ublocks = zeros( FMC.Probe.Nelements , ImgParams.NsingvalsRetain*(ImgParams.N_Xblocks*ImgParams.N_Zblocks) , kfmax-kfmin+1 , 'single' );
for k = 1 : kfmax-kfmin+1
    for i_block = 1 : ImgParams.N_Xblocks*ImgParams.N_Zblocks
        N    = FMC.Probe.Nelements;
        rank = ImgParams.NsingvalsRetain;
        Ublocks( 1:N , 1+(i_block-1)*rank : i_block*rank , k ) = single( U( 1+(k-1)*N : k*N , : , i_block ) );
    end
end

G_fk = zeros( FMC.Probe.Nelements , FMC.Probe.Nelements , kfmax-kfmin+1 , 'single' ); 
for k = kfmin : kfmax
    f_0     = FMC.Probe.TrandsucerCentralFrequency_MHz;
    f_k     = (k-1)/FMC.NtimePoints * FMC.SamplingFreqency_MHz;
    f_sigma = FMC.Probe.GaussFilterSigmaFrequency_MHz;
    
    i = k-kfmin+1;
    G_fk(:,:,i) = exp( -0.5 * ((f_k-f_0)/f_sigma)^2 ) * ones( FMC.Probe.Nelements , FMC.Probe.Nelements );
end


A_fft   = zeros( FMC.Probe.Nelements                               , FMC.Probe.Nelements                                                 , FMC.NtimePoints                         , 'single' );
A_fk    = zeros( FMC.Probe.Nelements                               , FMC.Probe.Nelements                                                 , kfmax-kfmin+1                           , 'single' );
Image   = zeros( ImgParams.N_Xblocks*ImgParams.N_pixPerXblock      , ImgParams.N_Zblocks*ImgParams.N_pixPerZblock                        ,                                           'single' );
Pixvals = zeros( ImgParams.N_pixPerXblock*ImgParams.N_pixPerZblock , ImgParams.N_Xblocks*ImgParams.N_Zblocks                             ,                                           'single' );
A_fkU_b = zeros( FMC.Probe.Nelements                               , ImgParams.NsingvalsRetain*(ImgParams.N_Xblocks*ImgParams.N_Zblocks) , kfmax-kfmin+1                           , 'single' );
A_fkU   = zeros( FMC.Probe.Nelements*(kfmax-kfmin+1)               , ImgParams.NsingvalsRetain                                           , ImgParams.N_Xblocks*ImgParams.N_Zblocks , 'single' );
UtA_fkU = zeros( ImgParams.NsingvalsRetain                         , ImgParams.NsingvalsRetain                                           , ImgParams.N_Xblocks*ImgParams.N_Zblocks , 'single' );

fftw('planner','patient');
maxNumCompThreads('automatic');


tic;

%% Calculating S_fk
A_fft = fft( FMC.Araw , [] , 3 );
A_fk  = A_fft( : , : , kfmin:kfmax ); 
A_fk  = A_fk .* G_fk;

%% Calculating pixel values
A_fkU_b = pagemtimes( A_fk , Ublocks );
A_fkU   = reshape( permute(A_fkU_b,[1,3,2]) , size(A_fkU) );
UtA_fkU = pagemtimes( U , 'transpose' , A_fkU , 'none' );
Pixvals = abs( squeeze( sum( V .* pagemtimes(UtA_fkU,V) , 1 ) ) );

%% Collecting pixel values into an image matrix
for ixblock = 1 : ImgParams.N_Xblocks
    for izblock = 1 : ImgParams.N_Zblocks
        i_block = (ixblock-1)*ImgParams.N_Zblocks + izblock;

        Image( (ixblock-1)*ImgParams.N_pixPerXblock+1 : ixblock*ImgParams.N_pixPerXblock , ...
               (izblock-1)*ImgParams.N_pixPerZblock+1 : izblock*ImgParams.N_pixPerZblock)  ...
              = reshape( Pixvals(:,i_block) , [ImgParams.N_pixPerXblock,ImgParams.N_pixPerZblock] )';
    end
end

toc;



%% Plotting
for threshold = linspace(250,10,20)
    figure('Renderer', 'painters', 'Position', [100 100 1500 800]); imagesc( arrayfun( @(x)min(x,threshold) , Image' ) ); colormap(jet(256));
    pause();
    close all; 
end


figure;
surf(Image);
colormap(jet(256));







