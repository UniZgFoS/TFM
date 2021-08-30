clear;
clc;
close all;
addpath('../Utility/');

load data.mat;


Probe = struct();
Probe.Nelements                      = 64;
Probe.ElementPitch_mm                = 0.63;
Probe.ElementWidth_mm                = 0.53;
Probe.ElementLength_mm               = 12;
Probe.TrandsucerCentralFrequency_MHz = 5.0;
Probe.GaussFilterSigmaFrequency_MHz  = 1;
Probe.UpperCutoffFrequency_MHz       = 7.0; 
Probe.LowerCutoffFrequency_MHz       = 3.0;

FMC = struct();
FMC.title                            = 'Bristol''s_FMC_with_64_elements_5MHz_probe';
FMC.Probe                            = Probe;
FMC.NtimePoints                      = 1000;
FMC.SamplingFreqency_MHz             = 50;
FMC.SpecimenUltrasonicSpeed_mmperus  = 6.3;



Araw      = zeros( FMC.Probe.Nelements , FMC.Probe.Nelements , FMC.NtimePoints , 'single' );
Afiltered = zeros( FMC.Probe.Nelements , FMC.Probe.Nelements , FMC.NtimePoints , 'single' );

iAscan = 0;
for rx = 1 : FMC.Probe.Nelements
    for tx = rx : FMC.Probe.Nelements        
        
        iAscan = iAscan + 1;
        Ascan  = exp_data.time_data( : , iAscan );
        
        dt_us             = 1 / FMC.SamplingFreqency_MHz;
        f0_MHz            = FMC.Probe.TrandsucerCentralFrequency_MHz;
        fsigma_MHz        = FMC.Probe.GaussFilterSigmaFrequency_MHz;
        f_upperCutoff_MHz = FMC.Probe.UpperCutoffFrequency_MHz;
        f_lowerCutoff_MHz = FMC.Probe.LowerCutoffFrequency_MHz;
        
        Ascan_filtered = function_filterGaussian( Ascan , dt_us , f0_MHz , fsigma_MHz , f_upperCutoff_MHz , f_lowerCutoff_MHz );
        
        
        Araw(tx,rx,:)      = single( Ascan          );
        Araw(rx,tx,:)      = single( Ascan          );
    
        Afiltered(tx,rx,:) = single( Ascan_filtered );
        Afiltered(rx,tx,:) = single( Ascan_filtered );
           
    end
end
Araw      = Araw      / max(max(max(abs(Araw))));
Afiltered = Afiltered / max(max(max(abs(Afiltered))));

FMC.Araw      = Araw;
FMC.Afiltered = Afiltered;

save('FMC','FMC');
movefile( "./FMC.mat" , "../FMCdatabase/"+FMC.title+".mat" );


function_plotSignal( Araw(5,10,:)      , 1 / FMC.SamplingFreqency_MHz );
function_plotSignal( Afiltered(5,10,:) , 1 / FMC.SamplingFreqency_MHz );






