// -------------------------------------------------------------------------------------- //
// 2024.08.30 Pre-Processing of Tree Ring Structure
// -------------------------------------------------------------------------------------- //
// @description: 
// This macro is to automate the pre-processing steps to filter, top-hat, DoG the CT Images.
// Noted that, 
// for all the detailed parameter settings, please adjust them according to the image.
// It would be recommanded to use CLIJ2 GUI interface for more interactive adjustments.
// -------------------------------------------------------------------------------------- //

// # Environment Settings...
// Initialize CLIJ2
run("CLIJ2 Macro Extensions", "cl_device=[NVIDIA GeForce GTX 1050 Ti]");
Ext.CLIJ2_clear();

// Open an image from folder
// Active this line when you do not want to input image to imageJ first.
// open();

// Get the title of the active image
image = getTitle();

// Get the directory of the input image
inputDir = getDirectory("image");
print(inputDir);

// # Pre-Processing
// median
Ext.CLIJ2_push(image);
im_median = "01_median";
radius_x = 2.0;
radius_y = 2.0;
Ext.CLIJ2_median2DSphere(image, im_median, radius_x, radius_y);
Ext.CLIJ2_pull(im_median);

// top hat
Ext.CLIJ2_push(im_median);
im_tophat = "02_top_hat_10";
radius_x = 10.0;
radius_y = 10.0;
radius_z = 10.0;
Ext.CLIJ2_topHatSphere(im_median, im_tophat, radius_x, radius_y, radius_z);
Ext.CLIJ2_pull(im_tophat);

// difference of gaussian
Ext.CLIJ2_push(im_tophat);
im_DoG = "03_DoG_210_TH10";
sigma1x = 2.0;
sigma1y = 2.0;
sigma2x = 10.0;
sigma2y = 10.0;
Ext.CLIJ2_differenceOfGaussian2D(im_tophat, im_DoG, sigma1x, sigma1y, sigma2x, sigma2y);
Ext.CLIJ2_pull(im_DoG);

// # Save the final image 
// in the same directory as the input image
selectWindow("01_median");
outputFileName = inputDir + im_median + ".tif";
saveAs("Tiff", outputFileName);

selectWindow("02_top_hat_10");
outputFileName = inputDir + im_tophat + ".tif";
saveAs("Tiff", outputFileName);

selectWindow("03_DoG_210_TH10");
outputFileName = inputDir + im_DoG + ".tif";
saveAs("Tiff", outputFileName);
