# BaumScan_Code-Delineation

**In Brief:** Delineation of tree ring and sapwood boundaries from various conditions of stem images.  

To capture the internal structures of tree stems for nondestructive analysis, we have developed delineation codes and protocols. Since images obtained from computed tomography (CT) scans share similarities with those from standard printer scanners in terms of patterns, colors, and data formats, we tested the applicability of our delineation procedures on common image data from stem discs.  

The detection codes are demonstrated using images from both "CT" and "Printer Scanner" sources. While the core code remains the same, some parameters are adjusted to account for differences in image resolution and condition (fresh versus air-dried samples). For specific details, please refer to the protocols in each folder.  

The following sections provide an overview of the software, plugins, packages, and image types used in this stem feature delineation project.  
  
.  
  
# 1. Requirements
## ImageJ  
### (1) Fiji
If Fiji is the choice of ImageJ flavor, please cite as:  
- Schindelin, J., Arganda-Carreras, I., Frise, E. et al. Fiji: an open-source platform for biological-image analysis. Nat Methods 9, 676–682 (2012). https://doi.org/10.1038/nmeth.2019  

otherwise, please refer to [ImageJ Docs](https://imagej.net/contribute/citing) for making correct and respectful citations.  

### (2) Plugins: CLIJ2
It is crucial to activate the ImageJ update sites from **"clij" & "clij2"** as they provide very useful and interactive environments to process the stem images.
**Please navigate in ImageJ: Help / Update... / ImageJ Updater / Manage update sites**, make sure the boxes of "clij" & "clij2" are ticked, then click "Apply changes".  

To check CLIJ2, please refers to: <https://clij.github.io/>, and reference:  
- Robert Haase, Loic Alain Royer, Peter Steinbach, Deborah Schmidt, Alexandr Dibrov, Uwe Schmidt, Martin Weigert, Nicola Maghelli, Pavel Tomancak, Florian Jug, Eugene W Myers. CLIJ: GPU-accelerated image processing for everyone. Nat Methods 17, 5-6 (2020) doi:10.1038/s41592-019-0650-1     
- Daniela Vorkel, Robert Haase. GPU-accelerating ImageJ Macro image processing workflows using CLIJ. arXiv preprint  
- Robert Haase, Akanksha Jain, Stéphane Rigaud, Daniela Vorkel, Pradeep Rajasekhar, Theresa Suckert, Talley J. Lambert, Juan Nunez-Iglesias, Daniel P. Poole, Pavel Tomancak, Eugene W. Myers. Interactive design of GPU-accelerated Image Data Flow Graphs and cross-platform deployment using multi-lingual code generation. bioRxiv preprint  

## R  
### (1) R & RStudio
The scripts are developed by R (ver. 4.3.1) & RStudio (ver. 2023.09.1+494 "Desert Sunflower"):  
- R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation
  for Statistical Computing, Vienna, Austria.
- RStudio Team (2020). RStudio: Integrated Development for R. RStudio, PBC, Boston, MA URL http://www.rstudio.com/

### (2) Packages: Imager, TSP, dplyr, purrr, OpenImageR
Key Packages including: ***Imager, TSP, dplyr, purrr, OpenImageR***, the relevant referencing as:
- Barthelme S (2023). _imager: Image Processing Library Based on 'CImg'_. R package version
  0.45.2, <https://CRAN.R-project.org/package=imager>.
- Hahsler M, Hornik K (2023). _TSP: Traveling Salesperson Problem (TSP)_. R package version
  1.2-4, <https://CRAN.R-project.org/package=TSP>.
- Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data
  Manipulation_. R package version 1.1.2, <https://CRAN.R-project.org/package=dplyr>.
- Wickham H, Henry L (2023). _purrr: Functional Programming Tools_. R package version 1.0.1,
  <https://CRAN.R-project.org/package=purrr>.
- Mouselimis L (2023). _OpenImageR: An Image Processing Toolkit_. R package version 1.3.0,
  <https://CRAN.R-project.org/package=OpenImageR>.   
  
## Example Imagries
Please note that the codes are designed for **grayscale images**.  
Hence, if images from normal printer scanners are with RGB format, make sure to convert them into greyscale. 
However, it is also recommanded to use the predefined option in office printer scanners to generate grayscale images other than color ones.  

### (1) Image from CT (Larch under Living Condition)
The image in this scenario comes from a larch stem disc, preserved and scanned using a stationary CT scanner with custom-tailored settings immediately after it was harvested from the forest. 
These settings were specifically designed and tested on the disc to support the development of the mobile CT (mCT) system, which will be implemented in the forest during future project phases.  

- **Image Name:** La-01.tif  
- **Location: Folder: "CT Images"**.  
- Resolution: 150 microns.  
- Delineation target: Tree Ring Boundaries + Sapwood Boundary  

### (2) Image from Printer Scanner (Spruce under Air-Dry-Condition)
The image used in this scenario comes from a spruce stem disc that was sanded and prepared for standard dendrochronological analysis. 
It was scanned using a commercial office printer scanner, Canon imageRUNNER ADVANCE C3561, set to a resolution of 600 dpi in grayscale.  

- **Image Name:** Fi_04.tif
- **Location: Folder: "Printer Scanner"**.  
- Resolution: 600 dpi (~42.33 microns).  
- Delineation target: Tree Ring Boundaries.  
  
.   
  
# 2. Protocol for Tree Ring Delineations
### The protocol is in general the same for both CT & Printer Scanned images.  
Please follow the instructions below and ensure that the relevant data is prepared in the correct formats.  
The process is ***semi-automatic*** and requires manual input to define the outermost tree ring structure using ImageJ. 
For detailed steps, please refer to the steps provided below.  

### Step 1: Delineate the outermost tree ring boundary and create the mask.  
- **Import Image:** Drag the targeted Stem Image into ImageJ.  
- **Create Mask:** Generate the file **"04_OuterRing_Mask.tif"**  
(Please refer to the **"Protocol for Tree Ring Marking via ImageJ.pdf"** for detailed instructions.)

### Step 2: Pre-Processing Image by ImageJ.  
- Select the Image opened in the ImageJ, then Navigate to **Plugins / Macros / Run...**
- Select the **Macro file: "Tree-Ring_Pre-Processing Steps.ijm"**  
The relevant data inputs for R Scripts will be generated & saved in the same folder as the image opened.  

**Note:**  
- Please edit the graphic card details in the Macro for CLIJ2. If you are unsure of the relevant statements, you can check them by navigating to **Plugins / Macros / Record...**,
  then selecting any function from **Plugins / ImageJ on GPU (CLIJ2)**.
- You can edit specific parameters by navigating to **Plugins / Macros / Edit...**, then selecting the **Macro file: "Tree-Ring_Pre-Processing Steps.ijm"**.  
It is recommended to use the interactive parameter adjustments provided by **CLIJ2-Assistant** from CLIJ2 when processing your own image.  
- Procedure differences between **CT** & **Printer Scanner** images can be found from their **Macro file: "Tree-Ring_Pre-Processing Steps.ijm"** in their folders.  

### Step 3: Compute tree ring structures and output as "ring_auto.csv".  
- Open RStudio and R  
- Run Script **"Detection Code.R"**  

**Note**  
- This step may take some time to complete.  
- Please **adjust the R script to match your personal working directory**.  
The working directory should be the same as the folder containing the input stem image and all processing scripts/macros.  
- Different parameter adjustments between **CT** & **Printer Scanner** images can be found from their **"Detection Code.R"** in their folders.

### Step 4: Re-Construct Tree Ring Structures in ImageJ.  
- Open ImageJ with the targeted stem image.  
- Navigate to **Plugins / Macros / Run...**
- Select the **Macro file: "R_2_ImageJ.ijm"**.
- Select the **"ring_auto.csv"** file (Results from R) and construct tree ring structures as ROIs in ROI Manager. 

### Step 5: Revise R-detected structures interactively using the “Selectin Brush Tool” .  
In **ImageJ:**  
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

### Step 6: Input updated tree ring structures back into R.  
In **ImageJ:**  
- Ensure the stem image is open with tree ring structures (ROIs in ROI Manager).
- Navigate to **ROI Manager /** Select all ROIs in ROI Manager   
- Navigate to **Plugins / Macros / Run...**  
- Select the **Macro file: "ExtractROICoordinate.ijm"**  
- **Save** the pop-up table as **"Results_Auto.csv"**  

In **R & RStudio:**  
- Open R-Script: **"Detection Code.R"**, and navigate to **section [2]-(4)_Import Adjusted structures**.  
- Run the code in this section.  
The lines are to translate the **CSV** file from ImageJ to R.  
Ensure the targeted stem image is also imported into R as "im" in matrix format.  
The exported structure in R will be either a data frame or a list that includes both a **data frame** and a **CImage object**.    

**Note:**  
If you are performing a **full manual delineation via ImageJ**,  
- **In ImageJ,** follow the instructions above to save the structure, and name the **.zip** file as **"Results.csv"**.  
- **In R & RStudio:** open **"Detection Code.R"**, and navigate to **section [2]-(5)_Import Manual Ring Positions (Optional)**.  
  
.  
  
# 3. Protocol for Sapwood Delineations
### The protocol is ONLY designed for CT images.  
Please follow the instructions below and ensure that the relevant data is prepared in the correct formats.  
The process is ***semi-automatic*** and requires manual input to define the outermost tree ring structure using ImageJ. 
For detailed steps, please refer to the steps provided below.  

### Step 1: Delineate the outermost tree ring boundary and create the mask.  
- **Import Image:** Drag the targeted Stem Image into ImageJ.    
- **Create Mask:** Generate the file **"04_OuterRing_Mask.tif"**  
(Please refer to the **"Protocol for Tree Ring Marking via ImageJ.pdf"** for detailed instructions.)

**Note:**  
- If the outermost tree ring boundary was already defined during the tree ring delineation, simply skip **Create Mask:** step. Instead,
- Import the created Outer-Ring Mask.
- This Outer-Ring Mask will be considered as the **"Outer Boundary"** of the sapwood.  
  (Further development might be possible to delineate the Outer Boundary from the Macro results...)  

### Step 2: Capture Sawpood Structure by ImageJ.  
- Select the Image opened in the ImageJ, then Navigate to **Plugins / Macros / Run...**
- Select the **Macro file: "Sapwood Structure Capture.ijm"**  
The relevant sapwood structure will be created & saved in the same folder as the image opened.  

**Note:**  
- Two types of results will be saved: **Sapwood Image in TIFF** & **Sapwood RoiSet.zip**.
- In **Sapwood RoiSet.zip**, please ***ONLY** make use of **"Inner Boundary"**.  
(The Outer Boundary is saved for further development opportunity to automate the delineation of the outermost Tree Ring.)  

### Step 3: Revise Sapwood Structures Interactively using the “Selectin Brush Tool” .  
- Select **ROI object:** "Inner Boundary" from the **ROI Manager** and review their accuracy.  

For any wrong delineations,  
- Select the ROI with mistakes.  
- Use the **“Selectin Brush Tool”** to adjust the results.  
- After making corrections, navigate to **ROI Manager / Update** to save the changes.  

When all adjustments are complete, **save the updated ROI object by:**  
- Navigate to **ROI Manager /** Select **ROI object:** "Inner Boundary" **/ More >> / Save...**
- The selected ROI will be saved as a .roi file.
- To reload the saved ROI Object, simply drag the .roi file into ImageJ with the stem image opened.  

### Step 4: Input Updated Sapwood Inner Boundary Structure into R.  
In **ImageJ:**  
- Ensure the stem image is open with "updated" sapwood structures (ROIs in ROI Manager).
- Navigate to **ROI Manager /** Select **ROI object:** "Inner Boundary" **
- Navigate to **Plugins / Macros / Run...**  
- Select the **Macro file: "ExtractROICoordinate.ijm"**  
- **Save** the pop-up table as **"Results_Auto.csv"**  

In **R & RStudio:**  
- Open R-Script: **"Detection Code.R"**, and navigate to **section [2]-(6)_Import Sapwood Structures**.  
- Run the code in this section.  
The lines are to translate the **CSV** file from ImageJ to R.  
Ensure the targeted stem image is also imported into R as "im" in matrix format.  
The exported structure in R will be either a data frame or a list that includes both a **data frame** and a **CImage object**.  

**Note:**  
If you are performing a **full manual delineation via ImageJ**,  
- **In ImageJ,** follow the instructions above to save the structure, and name the **.zip** file as **"Results.csv"**.  
- **In R & RStudio:** open **"Detection Code.R"**, and navigate to **section [2]-(7)_Import Manual Sapwood Positions (Optional)**.  
  
.  
  







