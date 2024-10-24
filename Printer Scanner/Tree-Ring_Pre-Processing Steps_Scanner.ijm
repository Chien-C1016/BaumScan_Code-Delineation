// -------------------------------------------------------------------------------------- //
// 2024.08.30 Pre-Processing of Tree Ring Structure
// -------------------------------------------------------------------------------------- //
// @description: 
// This macro is to automate the pre-processing steps to filter, DoG the CT Images.
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

// invert + add image and scalar
im_invert = "05_invert";
Ext.CLIJ2_push(im_median);
Ext.CLIJ2_invert(im_median, im_invert);
im_invert_scalar = "05_invert_scalar";
scalar = 255.0;
Ext.CLIJ2_addImageAndScalar(im_invert, im_invert_scalar, scalar)
Ext.CLIJ2_pull(im_invert_scalar)

// difference of gaussian
Ext.CLIJ2_push(im_invert_scalar);
im_DoG = "03_DoG_TH12";
sigma1x = 2.0;
sigma1y = 2.0;
sigma2x = 10.0;
sigma2y = 10.0;
Ext.CLIJ2_differenceOfGaussian2D(im_invert_scalar, im_DoG, sigma1x, sigma1y, sigma2x, sigma2y);
Ext.CLIJ2_pull(im_DoG);

// threshold otsu
Ext.CLIJ2_push(im_DoG);
im_otsu = "06_Thres_Otsu";
Ext.CLIJ2_thresholdOtsu(im_DoG, im_otsu);
Ext.CLIJ2_pull(im_otsu);

// # Save the final image 
// in the same directory as the input image
selectWindow("01_median");
outputFileName = inputDir + im_median + ".tif";
saveAs("Tiff", outputFileName);

selectWindow("03_DoG_TH12");
outputFileName = inputDir + im_DoG + ".tif";
saveAs("Tiff", outputFileName);

selectWindow("05_invert_scalar");
outputFileName = inputDir + im_invert + ".tif";
saveAs("Tiff", outputFileName);

selectWindow("06_Thres_Otsu");
outputFileName = inputDir + im_otsu + ".tif";
saveAs("Tiff", outputFileName);
