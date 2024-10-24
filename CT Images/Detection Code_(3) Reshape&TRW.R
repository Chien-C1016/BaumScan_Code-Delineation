# Environment Setting ----
# Package
require(dplyr)
require(purrr)
require(imager)

# Work directory
setwd("E:/PhD_WFH/La_03_R/R_TRW")

# Self-defined Function ----
source("AngularCompute.R")
source("ImageJ_2_R.R")
source("R_2_ImageJ.R")
source("XY Conversion & CImg Transformation.R")
source("Ring Distance Computation.R")

ring.BaseAdjust <- function(ring.b.list, reverse.base){
  
  .ring.b.list <-
    lapply(ring.b.list, 
           function(v){
             v$dist <- reverse.base - v$dist
             return(v)
           })
  return(.ring.b.list)
  
}

# Data ----
#' Image
im  = OpenImageR::readImage("01_median.tif") # Median_3D Processed
cim = im %>% imager::as.cimg()

#' Computed Pith Pos
pith <- readRDS("Pith.RData")

#' Manual Outer-ring Mask
mask.outer =
  OpenImageR::readImage("04_OuterRing_Mask.tif")%>%
  imager::as.cimg() %>%
  imager::imeval(~(.>0))
mask.df <-
  mask.outer %>% 
  imager::as.cimg() %>%
  mask.as.Ring() %>% 
  edge.FeatureCompute(xc = pith$x, 
                      yc = pith$y)

Ring.ref <-
  mask.df %>% 
  AngularCompute(xc = pith$x,
                 yc = pith$y,
                 df.mt = .,
                 theta.range = 5) %>% 
  purrr::pluck("df.clst") %>% 
  ring.ResampleSegments(df.clst = .,
                        pith = pith[c("x", "y")],
                        shrink_by = 0.01, # 0.005
                        theta.range = 5,
                        .plotcheck = FALSE) #TRUE
# plot(Ring.ref$theta, Ring.ref$dist, cex = 0.3)



# Computed Ring Positions ----
reverse.base <- readRDS("reverse.base.RData")

ring.b.auto <- 
  readRDS("ring.b.list.RData") %>% 
  ring.BaseAdjust(reverse.base = reverse.base) %>% 
  ring.WriteID()
ring.auto <- ring.b.auto %>% bind_rows()

# plot(ring.auto$theta, 
#      ring.auto$dist, 
#      col = as.factor(ring.auto$clst),
#      cex = 0.3)

ring.df <- 
  ring.auto %>% 
  ring.PositionBuild(polar.ring = ., pith = pith[c("x", "y")]) %>% 
  dplyr::select(., c("x", "y", "clst")) %>% 
  dplyr::distinct()

# plot(ring.df$x,
#      ring.df$y,
#      col = as.factor(ring.df$clst),
#      cex = 0.3)

# Over-lay
im.result <- ring.MatrixBuild(im.ref  = cim, 
                              ring.df = ring.df) %>% imager::as.cimg()
plot(cim)
imager::highlight(im.result > 0)

# Validate in Interactive platform "ImageJ"
# Export to ImageJ
#'@note 
#' The difference in between R and ImageJ
#' Transpose Matrix Read
#'@description 
#' The result from R could be better improved manually in ImageJ
#' Once it is validated and improved, 
#' the results can be imported back to R.
ring.ijm <- R2IJ(ring.df, im.matrix = cim)
write.csv(ring.ijm, file = "ring_auto.csv", row.names = F)

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


# Manual Ring Positions ----
manual.data <- read.csv("Results.csv")
ring.manual <- 
  manual.data %>% 
  dplyr::filter(Label != "Bark") %>% 
  Ij2R(., im.matrix = im)

# plot(ring.df$x, ring.df$y, col = "red", cex = 0.3)
# points(manual.data$X, manual.data$Y, cex = 0.3)

# Over-lay
im.manual <-
  ring.manual$rings.R %>% 
  dplyr::bind_rows() %>% 
  ring.MatrixBuild(im.ref  = cim, 
                   ring.df = .) %>% 
  imager::as.cimg()
plot(cim)
imager::highlight(im.manual > 0)



# Tree Ring Accuracy ----
#' Over-lay (size: 2000*2000)
plot(cim)
imager::highlight(im.result > 0)
imager::highlight(im.manual > 0, col = "orange")

#' Vertical distance between each other
i = 1
distcompute.list <- list()
# n.i = ring.df$clst %>% unique() %>% length()
n.i = ring.auto$rings.R %>% length()
while (TRUE) {
  
  ta.ring   <- ring.manual$rings.R[[i]]
  # ne.ring   <- ring.df %>% dplyr::filter(clst %in% i)
  ne.ring   <- ring.auto  $rings.R[[i]]
  temp.dist <- RingDistance(Target.ring = ta.ring, 
                            Next.ring   = ne.ring, 
                            center = c(pith$x, pith$y))
  temp.dist$clst <- i
  distcompute.list[[i]] <- temp.dist
  
  # i-loop
  i = i+1
  if(i > n.i){break}
  
}

# Check
i = 1
distcompute.list[[i]]$Ring.dist %>% hist()
distcompute.list[[i]]$point.cos %>% hist()
i = i + 1

# 2-dimentional Box-plot
source("R_Plot Results.R")
distcompute.list %>% TRW_2dimentionBoxplot()
distcompute.list %>% TRW_2dimentionBoxplot(relative.x = TRUE)
saveRDS(distcompute.list, file = "distcompute.RData")

# Check
# test <- readRDS("distcompute.RData")
# identical(distcompute.list, test)

res_accuracy <- 
  distcompute.list %>% 
  ring.RelativeDist(Target.ID = 1, .keep = TRUE)
write.csv(res_accuracy, file = "La03_res_accuracy.csv", row.names = F)

# Check
plot(res_accuracy$relative.position, res_accuracy$Ring.dist, cex = 0.3)
