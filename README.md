# BaumScan_Code-Delineation
Delineation of tree ring and sapwood boundaries from various conditions of stem images

# Description
The detection codes are demonstrated separately into "CT images" and "Printer Scanner". Their primary codes are the same but with some adjustments based on the image's resolution. Please follow the specific protocol in each folder for more specified info.

## Platforms & Required Packages / Plugins
### ImageJ  
#### (1) Fiji
If Fiji is the choice of ImageJ flavor(s), please cite as:  
- Schindelin, J., Arganda-Carreras, I., Frise, E. et al. Fiji: an open-source platform for biological-image analysis. Nat Methods 9, 676–682 (2012). https://doi.org/10.1038/nmeth.2019

otherwise, please refers to [ImageJ Docs](https://imagej.net/contribute/citing) for making correct and respectful citations.

#### (2) CLIJ2
It is crucial to activate the ImageJ update sites from **"clij" & "clij2"** as they provide very useful and interactive environments to process the stem images.
**Please navigate in ImageJ: Help / Update... / ImageJ Updater / Manage update sites**, make sure the boxes of "clij" & "clij2" are ticked, then click "Apply changes".  

To check CLIJ2, please refers to: <https://clij.github.io/>, and reference:  
- Robert Haase, Loic Alain Royer, Peter Steinbach, Deborah Schmidt, Alexandr Dibrov, Uwe Schmidt, Martin Weigert, Nicola Maghelli, Pavel Tomancak, Florian Jug, Eugene W Myers. CLIJ: GPU-accelerated image processing for everyone. Nat Methods 17, 5-6 (2020) doi:10.1038/s41592-019-0650-1  
- Daniela Vorkel, Robert Haase. GPU-accelerating ImageJ Macro image processing workflows using CLIJ. arXiv preprint  
- Robert Haase, Akanksha Jain, Stéphane Rigaud, Daniela Vorkel, Pradeep Rajasekhar, Theresa Suckert, Talley J. Lambert, Juan Nunez-Iglesias, Daniel P. Poole, Pavel Tomancak, Eugene W. Myers. Interactive design of GPU-accelerated Image Data Flow Graphs and cross-platform deployment using multi-lingual code generation. bioRxiv preprint  

### R  
### Note
Please cite plugins and libraries used in above platforms accordingly. It is really important to respect great efforts in the scientific field.

# Example Image Description 
## Resolution
(1) CT Image: 150 microns  
(2) Printer Scanner: 600 dpi (~42.33 microns)  
## note
As the code is defined for the greyscale computation, if the Printer images are color images, make sure to convert it. For several office scanners, there is an option to generate grayscale image.

## Requirements for computation
### [For TreeRing & Pith]
#' (1) Stem Disc Images.tif
#' (2) Outer Ring Mask.tif (marked in ImageJ)
### [For Sapwood]
#' (1) Stem Disc Images with Sapwood pattern.tif
#' (2) Outer Ring Mask.tif (marked in ImageJ)
