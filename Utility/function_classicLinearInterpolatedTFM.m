function Image = function_classicLinearInterpolatedTFM( FMC , Npix_x , Npix_z , ImageDepth_mm , ImageLength_mm , Status )
    
    TransducerPosition = zeros( 2 , FMC.Probe.Nelements );
    for iel = 1 : FMC.Probe.Nelements
         L = ImageLength_mm;
         p = FMC.Probe.ElementPitch_mm;
         w = FMC.Probe.ElementWidth_mm;
         N = FMC.Probe.Nelements;
         
         Probelen = p*(N-1)+w;
         
         xpos = (L-Probelen+w)/2 + (iel-1)*p;
         zpos = 0;
         
         TransducerPosition(:,iel) = [ xpos ; zpos ];
    end

    TimeDelays = zeros( FMC.Probe.Nelements , Npix_x , Npix_z );
    for ix = 1 : Npix_x
        for iz = 1 : Npix_z
            xpix = (ix-0.5) * ImageLength_mm / Npix_x;
            zpix = (iz-0.5) * ImageDepth_mm  / Npix_z;
            
            for iel = 1 : FMC.Probe.Nelements
                c = FMC.SpecimenUltrasonicSpeed_mmperus;
                TimeDelays( iel , ix , iz ) = norm( TransducerPosition(:,iel) - [xpix;zpix] ) / c;
            end
            
        end
    end

    
    dt = 1 / FMC.SamplingFreqency_MHz; % [us]
    
    Image = zeros( Npix_x , Npix_z );
    for ix = 1 : Npix_x
        if( Status )
            fprintf( '%d/%d\n' , ix , Npix_x );
        end
        for iz = 1 : Npix_z
            
            for iel1 = 1 : FMC.Probe.Nelements
                t1 = TimeDelays(iel1,ix,iz);
                
                for iel2 = 1 : FMC.Probe.Nelements
                    t2 = TimeDelays(iel2,ix,iz);
                    
                    it = (t1+t2)/dt;
                    it1 = int32(it);
                    it2 = it1 + 1;
                    
                    z1 = 0;
                    z2 = 0;
                    if( it2+1 <= FMC.NtimePoints )
                        z1 = FMC.Afiltered( iel1 , iel2 , it1+1 );
                        z2 = FMC.Afiltered( iel1 , iel2 , it2+1 );
                    end
                    
                    z = z1*double(it2-it) + z2*double(it-it1);
                    
                    Image(ix,iz) = Image(ix,iz) + z;
                end
                
                 
            end
            
        end
    end
    
    Image = abs(Image);
    
    
    return;
end

