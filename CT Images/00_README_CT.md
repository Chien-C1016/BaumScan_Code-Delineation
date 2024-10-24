# Protocol for Processing CT Image

**This protocol is designed specifically for CT images.**  
Please follow the instructions below and ensure that the relevant data is prepared in the correct formats.  
The process is ***semi-automatic*** and requires manual input to define the outermost tree ring structure using ImageJ. For detailed steps, please refer to the protocol provided below.  

# Platforms & Data Requirement
Please cite accordingly (please refers to the **README.md** at beginning of the repository).  

- ImageJ (Fiji)
- R
- Stem Disc Image (.tiff file)  

# Protocol  
## Step01
**Aim:** Delineate the outermost tree ring boundary and create the mask.  

- **Import Image:** Drag the targeted Stem Image into ImageJ.  
- **Create Mask:** Generate the file **"04_OuterRing_Mask.tif"**  
(Please refer to the **"Protocol for Tree Ring Marking via ImageJ.pdf"** for detailed instructions.)

## Step02
**Aim:** Pre-Processing Image by ImageJ.  

- Select the Image opened in the ImageJ, then Navigate to **Plugins / Macros / Run...**
- Select the **Macro file: "Tree-Ring_Pre-Processing Steps.ijm"**  
The relevant data inputs for R Scripts will be generated & saved in the same folder as the image opened.  

**Note:**  
- Please edit the graphic card details in the Macro for CLIJ2. If you are unsure of the relevant statements, you can check them by navigating to **Plugins / Macros / Record...**,
  then selecting any function from **Plugins / ImageJ on GPU (CLIJ2)**.
- You can edit specific parameters by navigating to **Plugins / Macros / Edit...**, then selecting the **Macro file: "Tree-Ring_Pre-Processing Steps.ijm"**.  
It is recommended to use the interactive parameter adjustments provided by **CLIJ2-Assistant** from CLIJ2 when processing your own image.  

## Step03
**Aim:** Run R scripts for computing tree ring structures and generate "ring_auto.csv".  

- Open RStudio and R  
- Run Script **"Detection Code.R"**  

**Note**  
- This step may take some time to complete.  
- Please **adjust the R script to match your personal working directory**.  
The working directory should be the same as the folder containing the input stem image and all processing scripts/macros.  

## Step04
**Aim:** 
#' Open ImageJ with the targeted stem image.
#' Run Macro to input tree ring structures as ROI files into ROI Manager.
#' Navigate to Plugins / Macros / Run...
#' Select the Macro file: "R_2_ImageJ.ijm".  

#'[Step06]
#' Use the interactive platform from ImageJ,
#' Select each tree ring structure from ROI Manager and 
#' check their structures accordingly. For any wrong delineation, 
#' use “Selectin Brush Tool” to adjust the R-detected results:
#' Select the ROI with mistakes / use “Selectin Brush Tool” for adjustments / 
#' ROI Manager / Update
#' Save the Updated ROI objects by:
#' ROI Manager / Select all ROIs of tree rings / More >> / Save...
#'[Step07] (Optional)
#' In case the updated tree ring structure data is needed back into R,
#' Navigate to Plugins / Macros / Run...
#' Select the Macro file: "ExtractROICoordinate.ijm"
#' Save the pop-up table as "Results.csv",
#' In R-Script "Detection Code.R", section [2]-(5),
#' There is supportive lines to translate the csv file from ImageJ to R.
#' Remember to also import the targeted Image into R as "im" as matrix format.
# ---------------------------------------------------------------------------- #
