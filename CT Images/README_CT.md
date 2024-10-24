# General Protocol for Processing CT Image
The protocol aims for "CT images".
Please follow the following instructions and preparing relevant data with correct data formats.
# ---------------------------------------------------------------------------- #
#'@Platform
#' (1) ImageJ (Fiji)
#' (2) R
#'@note
#' Please cite plugins and libraries used in above platforms accordingly.
#' It is really important to respect great efforts in the scientific field.
# ---------------------------------------------------------------------------- #
#'@Requirements for computation
#' (1) Stem Disc Images.tif
# ---------------------------------------------------------------------------- #
#'@Protocol 
#'[Step01]
#' Import Image: Drag the targeted Stem Image into ImageJ.
#'[Step02]
#' Create 04_OuterRing_Mask.tif.
#' Please follow: Protocol for Tree Ring Marking via ImageJ.
#'[Step03]
#' Pre-Processing Image by ImageJ.
#' Use the target Image opened in the ImageJ,
#' Navigate to Plugins / Macros / Run...
#' Select the Macro file: "Tree-Ring_Pre-Processing Steps.ijm"
#' The relevant data inputs for R Scripts will be generated & saved in the folder.
#' ps. Specific parameters can be edited by:
#' Navigate to Plugins / Macros / Edit...
#'[Step04]
#' Run R scripts for computing tree ring structures and generate "ring_auto.csv".
#' (1) This step takes some time;
#' (2) Please Adjust the R-script according to the personal work directory.
#'     The work directory should be the same as the folder with targeted stem image &
#'     all processing scripts/ Macros.
#'[Step05]
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
