# Environment Setting ----
# Package
require(dplyr)
require(purrr)
require(imager)

# Work directory
# Select the Image & Code directory
setwd(setwd("D:/Schreiben/Results/Fi_01"))

# Self-defined Scripts
source("Functions_Pith.R")
source("Functions_Ring-Clust.R")
source("Functions_Ring-Edge.R")
source("Functions_Segment-Grouping.R")
source("Function_Segment.Morph&Search.R")
source("Functions_Data Conversion.R")

# Self-defined Functions ----
## GetGradientval.v2 
#' Get unit gradient vector from the gradient output:
#'@description:
#  Gradient form computation can be 2 elements: 
#  "magnitude(intensity)" + "direction(Angular)"
#  Generally, the output can be used as direction for further spline usage
#'@note:
#  Such gradient can be compute in Cimage and Raster domain respectively. 
#  Choose the RIGHT one!!!
#'@import 
#  A gradient result from gradient computation (A CImage)
#'@export 
#  A dataframe structured for further spline approach 
#  (magnitude, unit dx, unit dy)
#'@source 
GetGradientVal <- function(gr, Edge, con.Raster = F){
  
  # ASSUMED the "xy" gradient had been computed and saved as a imlist
  # NAMED the image in the imlist FIRST !!!
  gr <- gr %>% setNames(c("x","y")) 
  # IMPORTANT when use get_gradient() OR USE imgradient() #
  # Codes modified from source code of cannyEdge and nonmax from imager
  
  # Magnitude and Gradient (dx,dy) preparation:
  #'@note: 3 CImages -> 3 dataframes
  .mag <- with(gr, sqrt(x^2 + y^2))
  
  # Thanks to the help of "list" property from imager
  if(con.Raster == F){
    .grs <- list( x = gr$x/.mag, y = gr$y/.mag)
  } else {
    # Reverse Y axes when consider the "AXES" change from CImage to Raster
    .grs <- list( x = gr$x/.mag, y = gr$y/.mag*(-1)) 
  }
  
  #'[EdgeMask Position preparation:]
  # Edge is a logi pixlsets indicates 
  # which cell values to acquire -> range [0,1] CImage
  # If the Edge is NOT a Pixset class, 
  # use Matrix2CImage() to convert it into CImage
  .EdgeMask <- Edge %>% imager::as.cimg()
  .df.Edge  <- .EdgeMask %>% as.data.frame()
  .edge.p   <- which(.df.Edge$value > 0) 
  
  # 3 Targeted df (magnitude, unit dx, unit dy):
  .df.mag <- (.mag   * .EdgeMask) %>% as.data.frame()
  .df.grX <- (.grs$x * .EdgeMask) %>% as.data.frame()
  .df.grY <- (.grs$y * .EdgeMask) %>% as.data.frame()
  
  # Acquire the values from the Edge positions
  .df <- data.frame(x   = .df.mag$x[.edge.p],
                    y   = .df.mag$y[.edge.p],
                    mag = .df.mag$value[.edge.p],
                    dx  = .df.grX$value[.edge.p],
                    dy  = .df.grY$value[.edge.p])
  
  return(.df)
}



# ---------------------------------------------------------------------------- #
# Data Input ----
# Image Import
im  = OpenImageR::readImage("01_median.tif") # Median_3D Processed
cim = im %>% imager::as.cimg()

mask.outer =
  OpenImageR::readImage("04_OuterRing_Mask.tif")%>%
  imager::as.cimg() %>%
  imager::imeval(~(.>0))
mask.cim <- mask.outer %>% imager::as.cimg()

# Edge Import
im.imj  = OpenImageR::readImage("03_DoG_210_TH10.tif")
imj.cim = (im.imj %>% imager::as.cimg()) * mask.cim

# ---------------------------------------------------------------------------- #
# Preparation ----
#'[Image:]
im.edge <- imj.cim %>% imager::imeval(~(.>0))
df.edge <- im.edge %>% imager::get.locations(., function(v) v>0)

gr    <- imager::get_gradient(cim, axes = "xy", scheme = 2)
df.gr <- GetGradientVal(gr, Edge = im.edge, con.Raster = F)



# ---------------------------------------------------------------------------- #
# [1] Pith Computation ----
#'[Pith: DBSCAN]
dat <- df.edge[,c("x", "y")]

eps_scale <- dat %>% dbscan.eps()
frnn_dat  <- dbscan::frNN(dat, eps = eps_scale)
db_dat    <- dbscan::dbscan(frnn_dat, minPts = 3)

df.db <- data.frame(x = dat$x,
                    y = dat$y,
                    clst = as.factor(db_dat$cluster))

df <- df.db %>% pith.CenterCheck()

# For checking and visualization (im.make = TRUE):
df.pith <- pith.GetPith(df, df.db = df.db, im.make = TRUE, im = im)

# Pith Computation
pith <- pith.ExtractPith(df.pith)

#'[Plot_Check:Not Run]
# Setting
.plotcheck = FALSE
if(.plotcheck == TRUE){
  
  plot(cim)
  imager::highlight(df.pith$im.pith>0)
  imager::circles(x = pith$x, y = pith$y, radius = 5, fg = "red")
  
}



# ---------------------------------------------------------------------------- #
# [2] Tree Ring Delineation ----
#'[Vector_Clean: Reduce unnecessary edges from noise]
#'[Vector Clean]
im.edge.vc <- edge.vectorclean(xc = pith$x,
                               yc = pith$y,
                               df.gradient = df.gr,
                               cos.sigma = 0.5,
                               mag.thres = 0.1,
                               as.mask = TRUE,
                               im.ref = cim)

#'[Get edge info from c-image]
df.edge <- 
  im.edge.vc %>% 
  imager::get.locations(., function(v) v>0) %>% 
  dplyr::select(., c("x", "y"))

df.edge %>% head()

#'[Perform_dbscan:]
eps_scale    <- 2 #df.edge %>% dbscan.eps()
frnn_df.edge <- dbscan::frNN(df.edge, eps = eps_scale)
db_df.edge   <- dbscan::dbscan(frnn_df.edge, minPts = 6)

#'[Plot_Check:Not Run]
# Setting
.plotcheck = FALSE
if(.plotcheck == TRUE){
  
  db_df.edge$cluster %>% unique() %>% length()
  plot(df.edge[,c("x", "y")], 
       col = as.factor(db_df.edge$cluster), 
       cex = 0.3)
  
}

#'[Filter the noise clusters identified by dbscan]
df.edge <-
  df.edge %>% 
  dplyr::mutate(clst = db_df.edge$cluster) %>% 
  dplyr::filter(clst != 0)

#'[Add MANUAL Outer-Ring data]
mask.cim
mask.df <- mask.cim %>% mask.as.Ring()

# Add manual outer tree ring structure
df.edge <- mask.join(df.edge = df.edge, mask.df = mask.df)

# Remove heavily clustered segments around "Pith"
# by radius to pith of 10 pixel
.radius2pith <- 10
df.edge <- mask.join(df.edge = df.edge, 
                     mask.df = data.frame(x = pith$x,
                                          y = pith$y,
                                          clst = 999999),
                     group.dist = .radius2pith)
df.edge <- df.edge %>% dplyr::filter(clst != 999999)

#'[Plot_Check:Not Run]
# Setting
.plotcheck = FALSE
if(.plotcheck == TRUE){
  
  plot(df.edge$x,
       df.edge$y,
       col = as.factor(df.edge$clst),
       cex = 0.3)
  
}

#'[Ring_Segments_Smoothing:Re-sample Tree Ring Segments]
# Re-sample Tree Ring Segments
theta.range = 5
ReSample.Size = 0.0001
ring.ReSample <-
  df.edge %>% 
  edge.FeatureCompute(xc = pith$x, 
                      yc = pith$y) %>% 
  edge.filter(., min.npts = 1) %>% 
  segment.AngularCompute(theta.range = theta.range) %>% 
  ring.ResampleSegments(pith = pith[c("x", "y")],
                        shrink_by = 0.01,
                        theta.range = theta.range,
                        AngularPositionUse = 3,
                        ReSample.Size = ReSample.Size)

#'[Plot_Check:Not Run]
# Setting
.plotcheck = FALSE
if(.plotcheck == TRUE){
  
  # Polarized System
  plot(ring.ReSample$theta,
       ring.ReSample$dist,
       col = as.factor(ring.ReSample$clst),
       cex = 0.3)
  
}



## (1)_Segment Grouping ----
## Preparation
#'[Preparation_01]
# ReSample.Size = 0.0001
# .digits = log10(ReSample.Size) %>% abs()
# .M.Step = ReSample.Size * 0.5
# .pi = round(pi, digits = .digits) %>% as.character() %>% as.numeric()

#'[Preparation_02]
#'@description 
#' (1) segment.size: How large should a big segment be.
#' (2) search.range: How much pixel range should the reference search begin.
#'                   If there is not enough segments to be considered, 
#'                   the threshold will be automatically enlarged.
segment.size <- 15
# search.range <- 25  

#'[Preparation_03]
#' Trusted Tree Ring Segments with Complete Structure
clst.trust <- 
  ring.ReSample %>% segment.ClustTrust(theta.range = theta.range)

#'[Preparation_04]
#' sapwood condition control
#'@note
#' This defines how many "true" segments should be included to
#' construct the reference tree ring.
sapw.control <- append(0.3, rep(0.6, (length(clst.trust)-1)))

#'[Preparation_05]
#' Reverse Control
#'@note
#' Originally, the code is designed to search from polar coordinates 
#' close to pith up-wards (search in between "trust Clusters).
#' And followed by a reverse call to process from the last inner 
#' trust tree ring towards pith.
#' However, as now there is a manual outer ring reference,
#' the segment grouping can start directly from the most outer part
#' inwards. Therefore, 
#' a "reverse" call is activated at beginning.
#'@signal
#' (1) reverse.base
#' (2) rev.control
#' (3) ring.ReSample.rev

rev.control <- TRUE
if(rev.control){
  
  message("# Search from Outer Ring to Pith #")
  reverse <- TRUE
  
  ring.2 <- ring.ReSample %>% dplyr::filter(., clst == clst.trust[1])
  ring.2 <- ring.2[,c("clst", "theta", "dist")]
  ring.2$clst <- "ref"
  
  .bp <- ring.2$dist %>% max() %>% log10() %>% floor()
  .bc <- (ring.2$dist %>% max() / 10^.bp) %>% ceiling()
  reverse.base <- .bc * 10^.bp
  
  ring.ReSample.rev <- ring.ReSample
  ring.ReSample.rev$dist <- reverse.base - ring.ReSample.rev$dist
  
}

#'[Segment_Grouping]
ring.rev.list <- segment.Morph.Search(ring.segments = ring.ReSample.rev,
                                      rev.control   = TRUE,
                                      theta.range   = theta.range,
                                      segment.size  = segment.size,
                                      ReSample.Size = ReSample.Size,
                                      sapw.control  = sapw.control)



## (2)_Structure Adjustments ----
#'[Structure Adjustments]
ring.auto <-
  ring.rev.list %>%
  purrr::pluck("ring.b")  %>% 
  ring.BaseAdjust(reverse.base = reverse.base) %>% 
  ring.WriteID() %>% 
  dplyr::bind_rows() %>% 
  ring.PositionBuild(polar.ring = ., pith = pith[c("x", "y")]) %>% 
  dplyr::select(., c("x", "y", "clst")) %>% 
  dplyr::distinct()

#'[Plot_Check:Not Run]
.plotcheck = FALSE
if(.plotcheck == TRUE){
  im.result <- 
    ring.MatrixBuild(im.ref  = cim, ring.df = ring.auto) %>% 
    imager::as.cimg()
  plot(cim)
  imager::highlight(im.result > 0)
}



## (3)_Output for Validation in "ImageJ" ----
# Export to ImageJ
#'@note 
#' The difference in between R and ImageJ
#' Transpose Matrix Read
#'@description 
#' The result from R could be better improved manually in ImageJ
#' Once it is validated and improved, 
#' the results can be imported back to R.
ring.ijm <- R2IJ(ring.auto, im.matrix = cim)
write.csv(ring.ijm, file = "ring_auto.csv", row.names = F)



## (4)_Import Adjusted structures ----
#'@note
#' Adjusted data frame named as "rings.R" within the list
#' generated by Ij2R.
#'[Import Adjusted structures]
R.auto.data <- read.csv("Results_Auto.csv")
ring.auto <-
  R.auto.data %>% 
  dplyr::filter(Label != "Bark") %>% 
  Ij2R(., im.matrix = im)

# Revise the "im.result" in R env.
im.result <-
  ring.auto$rings.R %>% 
  dplyr::bind_rows() %>% 
  ring.MatrixBuild(im.ref  = cim, 
                   ring.df = .) %>% 
  imager::as.cimg()



## (5)_Import Manual Ring Positions (Optional) ----
#'[Import Manual Ring Positions]
manual.data <- read.csv("Results.csv")
ring.manual <- 
  manual.data %>% 
  dplyr::filter(Label != "Bark") %>% 
  Ij2R(., im.matrix = im)

#'[Plot_Check:Not Run]
.plotcheck = FALSE
if(.plotcheck == TRUE){
  # Over-lay
  im.manual <-
    ring.manual$rings.R %>% 
    dplyr::bind_rows() %>% 
    ring.MatrixBuild(im.ref  = cim, 
                     ring.df = .) %>% 
    imager::as.cimg()
  plot(cim)
  imager::highlight(im.manual > 0)
}



## (6)_Import Sapwood Structures ----
sw.data <- read.csv("Results_SWAuto.csv")
sw.auto <-
  sw.data %>% 
  # In case all SW ROIs were extracted by "ExtractROICoordinate.ijm"
  dplyr::filter(Label == "Inner Boundary") %>%
  Ij2R(., im.matrix = .y)



## (7)_Import Manual Sapwood Positions (Optional) ----
sw.manual.data <- read.csv("Results_SW.csv")
sw.manual <-
  sw.manual.data %>% 
  Ij2R(., im.matrix = .y)


