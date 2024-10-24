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
**Aim:** Compute tree ring structures and output as "ring_auto.csv".  

- Open RStudio and R  
- Run Script **"Detection Code.R"**  

**Note**  
- This step may take some time to complete.  
- Please **adjust the R script to match your personal working directory**.  
The working directory should be the same as the folder containing the input stem image and all processing scripts/macros.  

## Step04
**Aim:** Re-Construct Tree Ring Structures in ImageJ.  

- Open ImageJ with the targeted stem image.  
- Navigate to **Plugins / Macros / Run...**
- Select the **Macro file: "R_2_ImageJ.ijm"**.
- Select the **"ring_auto.csv"** file (Results from R) and construct tree ring structures as ROIs in ROI Manager. 

## Step05
**Aim**: Revise R-detected structures interactively using the “Selectin Brush Tool” .  

In general, in **ImageJ:**  
- Continue from Step04,  
- Select each tree ring structure (**ROI object**) from the **ROI Manager** and review their accuracy.  

For any wrong delineations,  
- Select the ROI with mistakes.  
- Use the **“Selectin Brush Tool”** to adjust the R-detected results.  
- After making corrections, navigate to **ROI Manager / Update** to save the changes.  

When all adjustments are complete for all tree rings (ROIs), **save the updated ROI objects by:**  
- Navigate to **ROI Manager /** Select all ROIs in ROI Manager **/ More >> / Save...**
- The selected ROIs will be saved as a .zip file.
- To reload the saved ROIs, simply drag the .zip file into ImageJ with the stem image opened.  

## Step06 (Optional)
**Aim:** Input updated tree ring structures back into R.  

In **ImageJ:**  
- Ensure the stem image is open with tree ring structures (ROIs in ROI Manager).  
- Navigate to **Plugins / Macros / Run...**  
- Select the **Macro file: "ExtractROICoordinate.ijm"**  
- **Save** the pop-up table as **"Results.csv"**  

In **R & RStudio:**  
- Open R-Script: **"Detection Code.R"**, and navigate to **section [2]-(5)**.  
- Run the code in this section.  
The lines are to translate the **CSV** file from ImageJ to R.  
Ensure the targeted stem image is also imported into R as "im" in matrix format.  
The imported structure can be a data frame or a list that includes both a **data frame** and a **CImage object**.


