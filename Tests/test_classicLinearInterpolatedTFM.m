clear;
clc;
close all;

load('../FMCdatabase/Bristol''s_FMC_with_64_elements_5MHz_probe.mat');


Npix_x         = 250;
Npix_z         = 250;
ImageDepth_mm  = 60;
ImageLength_mm = 45;

lambda_mm      = FMC.SpecimenUltrasonicSpeed_mmperus / FMC.Probe.TrandsucerCentralFrequency_MHz;
pixel_xsize_mm = ImageLength_mm / Npix_x;
pixel_zsize_mm = ImageDepth_mm  / Npix_z;

fprintf( 'lambda      = %7.3f mm\n'                 , lambda_mm                                     );
fprintf( 'pixel_xsize = %7.3f mm = 1 / %f lambda\n' , pixel_xsize_mm , 1/(pixel_xsize_mm/lambda_mm) );
fprintf( 'pixel_zsize = %7.3f mm = 1 / %f lambda\n' , pixel_zsize_mm , 1/(pixel_zsize_mm/lambda_mm) );


tic;
Image = function_classicLinearInterpolatedTFM( FMC , Npix_x , Npix_z , ImageDepth_mm , ImageLength_mm , true );
toc;


figure;
imagesc(Image');
colormap(jet(256));

figure;
surf(Image);
colormap(jet(256));













