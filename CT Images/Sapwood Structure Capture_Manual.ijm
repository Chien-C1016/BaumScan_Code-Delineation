// -------------------------------------------------------------------------------------- //
// 2024.12.10 Sapwood Structure Capture (Manual)
// -------------------------------------------------------------------------------------- //
// @description: 
// This macro is to semi-automate the processing steps to:
// blur, manual supported threshold, ROI indexing
// the sapwood region. Noted that, in manual-suporting step, 
// adjust the threshold in the sliding bar, click 'Apply', 
// agree with any pop-ups from the threshold, and Click OK to continue the code.
// -------------------------------------------------------------------------------------- //
// # Environment Settings...
// Initialize CLIJ2
run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();

// Open an image
// open("path/to/your/image.tif");

// Get the title of the active image
inputImage = getTitle();

// Get the directory of the input image
inputDir = getDirectory("image");

// Push the image to GPU memory
Ext.CLIJ2_push(inputImage);

// # Padding ...
// Get the dimensions of the input image
width  = getWidth();
height = getHeight();

// Define the padding size
paddingSize = 400;

// Calculate new dimensions
newWidth  = width  + 2 * paddingSize;
newHeight = height + 2 * paddingSize;

// Create a new image with the new dimensions, filled with zeros
paddedImage = "paddedImage"
Ext.CLIJ2_create2D(paddedImage, newWidth, newHeight, 32); // 32-bit float image

// Copy the original image to the center of the new image
Ext.CLIJ2_paste2D(inputImage, paddedImage, paddingSize, paddingSize);

// Pull the processed image back to CPU memory
// Ext.CLIJ2_pull(paddedImage);

// Show the result
// run("Tile");

// # Perform Gaussian blur ...
gaussian_blur_r20 = "gaussian_blur_r20";
sigma_x = 15.0;
sigma_y = 15.0;
Ext.CLIJ2_gaussianBlur2D(paddedImage, gaussian_blur_r20, sigma_x, sigma_y);

// Pull the blurred image back to CPU memory
Ext.CLIJ2_pull(gaussian_blur_r20);

// Show the result
// run("Tile");

// Trigger the threshold adjustment interface
// *** Wait for the user to adjust and apply the threshold
run("Duplicate...", "title=threshold_IMG");
threshold_IMG = "threshold_IMG";
selectWindow("threshold_IMG");
run("Threshold...");
waitForUser("Please adjust the threshold and click 'Apply'. Click OK here to continue.");

if (isOpen("Threshold")) {
    selectWindow("Threshold");
    run("Close");
}

// Remove "gaussian_blur_r20" viewing window
selectWindow("gaussian_blur_r20");
close();

// Remove the padding using ImageJ's crop function
selectWindow("threshold_IMG");
makeRectangle(paddingSize, paddingSize, width, height);
run("Crop");

// Push the cropped image back to GPU memory
Ext.CLIJ2_push(threshold_IMG);

// Show the final result
// run("Tile");

// # Get Sapwood boundaries (inner and outer) 
// Closing: fill the tiny gaps in Sapwood
closingImage = "closingImage";
number_of_iterations = 2;
Ext.CLIJ2_closingBox(threshold_IMG, closingImage, number_of_iterations);

// connected Components Labeling Diamond
// Pull Labels to ROIManager
labelmap = "labelmap";
Ext.CLIJ2_connectedComponentsLabelingDiamond(closingImage, labelmap);
Ext.CLIJ2_pullLabelsToROIManager(labelmap);

// Get the ROIs and their areas
// If the sapwood polygon > 2, find the one with the max area
roiManager("Deselect");
maxAreaIndex = 0;
if(roiManager("count") > 1){
	
	// Measure all ROIs	
	roiManager("Measure");
	
	// Find the largest ROI
	maxArea = 0;
	maxAreaIndex = -1;
	for (i = 0; i < roiManager("Count"); i++) {
		area = getResult("Area", i);
    	if (area > maxArea) {
        	maxArea = getResult("Area", i);
        	maxAreaIndex = i;
    	}
	}
	
	// Select the largest ROI
    roiManager("Select", maxAreaIndex);

    // Remove all other ROIs
    for (i = roiManager("count") - 1; i >= 0; i--) {
        if (i != maxAreaIndex) {
            roiManager("Select", i);
            roiManager("Delete");
        }
    }
	
}

// Close Measurement Results
close("Results");

// Ensure the Default Size of the Image
selectWindow("threshold_IMG");

// Now select the largest area ROI
// Only the largest one left, position = 0 in ImageJ
maxAreaIndex = 0;
roiManager("Select", maxAreaIndex);
roiManager("Rename", "Sapwood area");

// Generate Sapwood Mask
run("Create Mask");
rename("Sapwood");

// Split the ROI of Sapwood:
// Sapwood polygon is like a donat, 
// for the calculation in R, we need inner and outer boundaries.
roiManager("Select", maxAreaIndex);
roiManager("Split");

// The ROI Manager should now contain two ROIs (inner and outer boundary)
// the first one is the inner and the second is the outer boundary
// ** There might be full sapwood condition in the future... be careful
if (roiManager("Count") > 1) { // Check if there are new ROIs after splitting
    roiManager("Select", maxAreaIndex + 1); // Select the first ROI resulting from the split
    roiManager("Rename", "Inner Boundary");
    roiManager("Select", maxAreaIndex + 2); // Select the second ROI resulting from the split
    roiManager("Rename", "Outer Boundary");
}

// Remove "threshold_intermodes" viewing window
selectWindow("threshold_IMG");
close();
// print("Done processing.");

// # Save the final image 
// in the same directory as the input image
selectWindow("Sapwood");
outputFileName = inputDir + "Sapwood_" + inputImage + ".tif";
saveAs("Tiff", outputFileName);

// # Save ROI files
roiManager("select", newArray(0, 1, 2));
outputROIName = inputDir + "RoiSet_" + inputImage + ".zip";
roiManager("Save", outputROIName);

