# BaumScan_Code-Delineation

**In Brief:** Delineation of tree ring and sapwood boundaries from various conditions of stem images.  

To capture the internal structures of tree stems for nondestructive analysis, we have developed delineation codes and protocols. Since images obtained from computed tomography (CT) scans share similarities with those from standard printer scanners in terms of patterns, colors, and data formats, we tested the applicability of our delineation procedures on common image data from stem discs.  

The detection codes are demonstrated using images from both "CT" and "Printer Scanner" sources. While the core code remains the same, some parameters are adjusted to account for differences in image resolution and condition (fresh versus air-dried samples). For specific details, please refer to the protocols in each folder.  

The following sections provide an overview of the software, plugins, packages, and image types used in this stem feature delineation project.  

# 1. Platforms & Required Packages / Plugins
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

# 2. Example Imagries
Please note that the codes are designed for **grayscale images**.  
Hence, if images from normal printer scanners are with RGB format, make sure to convert them into greyscale. However, it is also recommanded to use the predefined option in office printer scanners to generate grayscale images other than color ones.  

## Scenario 1: Image from CT (Larch under living-condition)
The image in this scenario comes from a larch stem disc, preserved and scanned using a stationary CT scanner with custom-tailored settings immediately after it was harvested from the forest. 
These settings were specifically designed and tested on the disc to support the development of the mobile CT (mCT) system, which will be implemented in the forest during future project phases.  

- Resolution: 150 microns.
- Delineation target: Tree Ring Boundaries + Sapwood Boundary  

## Scenario 2: Image from Printer Scanner (Spruce under Air-Dry-Condition)
The image used in this scenario comes from a spruce stem disc that was sanded and prepared for standard dendrochronological analysis. 
It was scanned using a commercial office printer scanner, Canon imageRUNNER ADVANCE C3561, set to a resolution of 600 dpi in grayscale.  

- Resolution: 600 dpi (~42.33 microns).  
- Delineation target: Tree Ring Boundaries.  

## Requires images for computations
### (1) For Tree Rings & Pith:
- Stem Disc Image (.tif).  
- Outer Ring Mask (.tif, masked in ImageJ).  

### (2) For Sapwood:
- Stem Disc Images with Sapwood pattern (.tif).  
- Outer Ring Mask (.tif, masked in ImageJ).  
