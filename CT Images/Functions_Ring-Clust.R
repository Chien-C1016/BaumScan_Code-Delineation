# ---------------------------------------------------------------------------- #
# Functions with Tree ring cluster segments #
# ---------------------------------------------------------------------------- #
# General ----
#'@description 
#' (1) Working on tree ring cluster segment objects.
#' (2) provide conversion from Mask to Ring;
#' (3) Polar system to Cartesian coordinate system (or vice versa)*
#' (4) Cartesian coordinate system to Image (vice versa is handled by imager)

# ---------------------------------------------------------------------------- #
# Functions ----
## Build ring Positions 
#'@description 
#' Convert a ring segment to Cartesian coordinate system, by
#' pith location or without but with reference tree ring data
#'@source 
ring.PositionBuild <- function(pith = NULL,
                               polar.ring,
                               reference.df = NULL){
  
  # Package
  require(purrr)
  require(dplyr)
  
  #'[Pith:Condition Check]
  if(purrr::is_empty(pith) & purrr::is_empty(reference.df)){
    
    stop("Please provide either Pith or reference ring data")
    
  }
  
  if(purrr::is_empty(pith) & 
     (purrr::is_empty(reference.df$x) | purrr::is_empty(reference.df$dist))){
    
    stop("Please provide the correct reference data frame with: ",
         "\n(1) pixel position in x-y format",
         "\n(2) polar position in theta-dist format",
         "\n(3) Targeted tree ring ID to be processed")
    
  }
  
  if(purrr::is_empty(pith)){
    
    message("Estimation of Pith is derived & used from reference ring data")
    
    # Estimate the Pith Position from reference ring data
    temp.refring <- 
      reference.df %>% 
      dplyr::mutate(.,
             .keep = c("none"),
             dx = dist * cos(theta),
             dy = dist * sin(theta),
             x  = (x - dx) %>% round(),
             y  = (y - dy) %>% round())
    
    pith <- c(temp.refring$x %>% mean(),
              temp.refring$y %>% mean())
    
  }

  pith
  # Fix pith data structure
  if(is.data.frame(pith)){
    
    pith <- 
      pith %>% dplyr::select("x", "y") %>% unlist() %>% unname()
    
  }
  if(is.list(pith)){
    
    pith <- pith %>% unlist() %>% unname()
    
  }
  
  #'[Cartesian_coordinate_system:]
  ring.pos <- 
    polar.ring %>% 
    dplyr::mutate(x = (pith[1] + dist * cos(theta)) %>% round(),
                  y = (pith[2] + dist * sin(theta)) %>% round())
  
  # Select and unify the columns
  ring.pos <- ring.pos[,append(c("x", "y"), names(polar.ring))]
  
  #'[Output:]
  return(ring.pos)
  
}



## Re-Sample Ring Segments
#'@description 
#  The internal point distance is not equal.
#  It is actually getting larger towards pith.
#'@method 
#  (1) Angular Window Computation
#  Spline in between certain Angular Positions
#  The Window Move 1 Angular Position for each step
#'@import 
#' (1) AngularPositionUse 
#' the size of the Search Window
#' If the Segment has Angular Position less than the window,
#' the overall spline in such Segment is applied.
#' (2) ReSample.Size
#' How fine scale does the re-sampling aiming at
#' (3) shrink_by
#' shrink re-sampled ring segments to avoid strange behavior from segment edges
#' (4) theta.range
#' How large theta range (radius) for angular grouping
#' (5) pith
#' Pith positions for "ring.PositionBuild"
#' (6) .plotcheck
#' An option for always drawing the smoothing between 0 and 2pi
#' (7) .plotexport
#' Export PNG figures for each tree ring segments
#'@source
ring.ResampleSegments <- function(df.clst,
                                  theta.range = 5,
                                  pith        = vector(),
                                  shrink_by   = vector(),
                                  .plotcheck  = FALSE,
                                  .plotexport = FALSE,
                                  AngularPositionUse = 3,
                                  ReSample.Size      = 0.0001){
  
  # Package
  require(dplyr)
  require(purrr)
  
  # Ring.ID
  Ring.ID <- df.clst$clst %>% unique()
  
  # Corresponding Round digits to ReSample.Size
  .digits = log10(ReSample.Size) %>% abs()
  
  #'[Adjustments:range of Angular Blocks]
  if(as.integer(360/theta.range) != 360/theta.range){
    
    while (as.integer(360/theta.range) != 360/theta.range) {
      theta.range = theta.range + 1
    }
    
    # Report Adjustments:
    message(paste0("## Adjustment: theta.range = ", theta.range))
    
  }
  
  # Angular Blocks
  n.theta <- 360/theta.range
  
  #'[Angular_Compute:Angular Blocks]
  theta.range
  n.theta
  Ang.pos <- replicate(length(df.clst$dist), n.theta)
  i = 1
  while (TRUE) {
    
    .theta <- 2*pi - i*theta.range/180*pi
    temp.p <- which(df.clst$theta < .theta)
    Ang.pos[temp.p] <- (n.theta-i)
    
    # Loop
    i = i+1
    if(i == n.theta){break}
    
  }
  df.clst$Ang.pos <- Ang.pos
  
  #'[Ring.Summary: Summary of Tree Ring Segments]
  Ring.sum <-
    df.clst %>% 
    dplyr::group_by(clst) %>% 
    dplyr::summarise(npts = n(),
                     Ang.min   = min(Ang.pos),
                     Ang.max   = max(Ang.pos),
                     Ang.range = Ang.max - Ang.min + 1)
  
  
  # -------------------------------------------------------------------------- #
  #'[Loop_Control:]
  # rs controls which Tree Ring Segment to be re-sampled
  # i  controls the position to Spline
  # k  controls the feed of Angular Position
  
  # Control the digits of pi
  .pi <- round(pi, digits = .digits)
  # -------------------------------------------------------------------------- #
  
  # rs-loop
  rs  = 1
  nrs = length(Ring.ID)
  ring.ReSample <- list()
  while (TRUE) {
    
    # "rs-th" Ring Segment df
    df.rs <- 
      df.clst %>% 
      dplyr::filter(clst == Ring.ID[rs]) %>% 
      dplyr::mutate(theta = theta %>% round(., digits = .digits))
    
    # rs Ring.Summary
    rs.rs <- Ring.sum %>% dplyr::filter(clst == Ring.ID[rs])
    
    # rs-Ring Segment Check
    ## For segments that cross between 0 and 2pi, padding
    ## Signal rw.check is indicated for removing padding.
    Ring.ID.rs <- df.rs$Ang.pos %>% unique()
    rw.check = FALSE
    if(length(Ring.ID.rs) != rs.rs$Ang.range &
       rs.rs$Ang.min == 1 &
       rs.rs$Ang.max == n.theta){
      
      message("\n", "RingID ", rs, " is Cut by theta = 0\n")
      
      ## Where the Ring Segment breaks
      k = 1L
      while (TRUE) {
        if((rs.rs$Ang.max - k) %in% Ring.ID.rs){
          k = k + 1
        }else{
          break
        }
      }
      Ang.min <- rs.rs$Ang.max - k + 1
      
      ## Re-Write df.rs (Padding)
      .rwp <- which(df.rs$Ang.pos < Ang.min)
      df.rs$theta[.rwp]   <- df.rs$theta[.rwp] + 2 * .pi
      df.rs$Ang.pos[.rwp] <- df.rs$Ang.pos[.rwp] + n.theta
      
      ## Re-Write rs.rs
      rs.rs <-
        df.rs %>% 
        group_by(clst) %>% 
        summarise(npts = n(),
                  Ang.min   = min(Ang.pos),
                  Ang.max   = max(Ang.pos),
                  Ang.range = Ang.max - Ang.min + 1)
      
      ## Signal
      rw.check = TRUE
      
    }
    
    # How many points per Angular range
    rs.rs <-
      rs.rs %>% 
      dplyr::group_by(clst) %>% 
      dplyr::mutate(Ang.npts = npts / Ang.range)
    
    # Segment Edges of rs (ensure format by "ReSample.Size"&.digits)
    theta.l <- df.rs$theta %>% min() %>% round(., digits = .digits)
    theta.r <- df.rs$theta %>% max() %>% round(., digits = .digits)
    
    # "shrink_by" to shrink the segments
    # Only when segment size is large enough or shrink_by is indicated
    if(purrr::is_empty(shrink_by) |
       (2*shrink_by / (theta.r - theta.l)) > 0.1){
      
      seq.rs  <- seq(theta.l, theta.r, by = ReSample.Size)
      
    }else{
      
      if(shrink_by > 1){shrink_by <- round(shrink_by)*ReSample.Size}
      seq.rs  <- seq((theta.l + shrink_by), 
                     (theta.r - shrink_by), 
                     by = ReSample.Size)
      
    }
    
    # rs segment control (points and Angular blocks)
    npts.rs <- length(seq.rs)
    
    # Re-mark the Ang.pos for reference
    Ang.prs <- rep(rs.rs$Ang.max, npts.rs)
    if(rs.rs$Ang.range > 1){
      
      k  = 1
      nk = rs.rs$Ang.range
      while (k <= nk) {
        
        thres.angpos <- (rs.rs$Ang.max - k) * theta.range / 180 * .pi
        remark.p     <- which(seq.rs <= (thres.angpos + ReSample.Size/10))
        Ang.prs[remark.p] <- (rs.rs$Ang.max - k)
        
        # loop
        k = k + 1
        
      }
    }
    
    # Smooth.spline ----
    if(rs.rs$Ang.range <= AngularPositionUse |
       rs.rs$Ang.npts < 20){
      ## Small Ring Segments ----
      
      # Fix
      #'@note
      #' For some segments, there are really small deviations (~0.0001)
      #' These do not contribute "unique" theta range to a considered segment
      len.unique_theta <-
        length(unique(round(df.rs$theta, digits = (.digits-1))))
      
      #'@note 
      #' If the segments are too small (npts <= 3)
      #' the Re-Sampling will use lm
      if(len.unique_theta == 1){
        
        print(paste("Only 1 point exists, skip"))
        temp.df.rs <- data.frame(clst    = df.rs$clst,
                                 theta   = df.rs$theta,
                                 dist    = df.rs$dist,
                                 Ang.pos = df.rs$Ang.pos)
        
      }else if(len.unique_theta <= 4){
        
        # Linear lm
        lm.rs <- lm(dist ~ theta, data = df.rs)
        
        # Re-sample the Ring Segment
        dist.rs    <- predict(lm.rs, data.frame(theta = seq.rs)) %>% unname()
        temp.df.rs <- data.frame(clst    = Ring.ID[rs],
                                 theta   = seq.rs,
                                 dist    = dist.rs,
                                 Ang.pos = Ang.prs)
        
      }else if(rs.rs$Ang.npts < 20 |
               dist(range(df.rs$theta)) < 0.2){
        
        # Smooth.spline
        smsp.rs <- smooth.spline(x = df.rs$theta,
                                 y = df.rs$dist,
                                 df = 3)
        
        # Re-sample the Ring Segment
        dist.rs    <- predict(smsp.rs, x = seq.rs)$y
        temp.df.rs <- data.frame(clst    = Ring.ID[rs],
                                 theta   = seq.rs,
                                 dist    = dist.rs,
                                 Ang.pos = Ang.prs)
        
      }else{
        
        # Smooth.spline
        smsp.rs <- smooth.spline(x = df.rs$theta,
                                 y = df.rs$dist,
                                 df = 5)
        
        # Re-sample the Ring Segment
        dist.rs    <- predict(smsp.rs, x = seq.rs)$y
        temp.df.rs <- data.frame(clst    = Ring.ID[rs],
                                 theta   = seq.rs,
                                 dist    = dist.rs,
                                 Ang.pos = Ang.prs)
        
      }
      
    }else{
      ## Large Ring Segments ----
      
      # Only smooth.spline within the search Window(AngularPositionUse)
      seq.window <- c(1L:AngularPositionUse)
      
      # i-loop (Angular Blocks to be Splined & Merged)
      i = 1
      ring.seg.api <- list()
      while (TRUE) {
        
        #'[Smoothing:]
        # Angular Positions for Smooth.spline
        .smsp.api <- seq.int(from = (rs.rs$Ang.min + (i-1)),
                             along.with = seq.window)
        df.api    <- df.rs %>% dplyr::filter(., Ang.pos %in% .smsp.api)
        
        # Smooth.spline
        #'@note Size-depended Degree of Freedom Usage
        #'@method npts >= 100 : df == 5
        if(nrow(df.api) < 100){ # use 100 for none-dist_resampled data
          
          smsp.api  <- smooth.spline(x  = df.api$theta,
                                     y  = df.api$dist,
                                     df = 3)
          
        }else{
          
          smsp.api  <- smooth.spline(x  = df.api$theta,
                                     y  = df.api$dist,
                                     df = 5)
          
        }
        
        # Re-sample the Ring Segment
        dist.api   <- predict(smsp.api, x = seq.rs)$y
        temp.df.rs <- data.frame(clst    = Ring.ID[rs],
                                 theta   = seq.rs,
                                 dist    = dist.api,
                                 Ang.pos = Ang.prs)
        
        #'[Angular_Window_Computation:]
        #'[ap.i Window Section]
        # sec01 (Overlapping Area)
        df.i.sc01 <- temp.df.rs %>% 
          dplyr::filter(., Ang.pos == .smsp.api[1])
        
        # sec02 
        df.i.sc02 <- temp.df.rs %>% 
          dplyr::filter(., Ang.pos == .smsp.api[2])
        
        # secNN (Rest Area)
        df.i.scNN <- temp.df.rs %>% 
          dplyr::filter(., Ang.pos %in% .smsp.api[-c(1,2)])
        
        #'[Process01]
        ring.seg.api[[i+1]] <- df.i.sc02
        
        #'[Process02 - Merging 2 Overlapping Section]
        # (1) ring.seg.api[[i]] (Previous df.i.sc02)
        # (2) Current df.i.sc01
        #'@method 
        # (1) Distance based weighting
        # (2) keep the edge.l from ring.seg.api[[i]]
        # (3) keep the edge.r from sec01
        if(i == 1){
          
          ring.seg.api[[i]] <- df.i.sc01
          
        }else{
          
          seq.api <- df.i.sc01$theta
          D = max(seq.api) - min(seq.api)
          
          # Merging 2 Spline Segments
          M.dist <- 
            ring.seg.api[[i]]$dist * 
            ((max(seq.api)-ring.seg.api[[i]]$theta) / D) + 
            df.i.sc01$dist * 
            ((df.i.sc01$theta - min(seq.api)) / D)
          
          ring.seg.api[[i]]$dist <- M.dist
          
        }
        
        #'[Process03 - Ending]
        if((i+(AngularPositionUse-1)) == rs.rs$Ang.range){
          
          ring.seg.api[[i+2]] <- df.i.scNN
          break
          
        }
        
        # i-Loop
        i = i + 1
        
      }
      
      temp.df.rs <- ring.seg.api %>% dplyr::bind_rows()
      
      # Padding Spline Correction
      if(rs.rs$Ang.range == n.theta){
        
        #' [Previous Smooth.spline Product]
        #' keep:
        #' (1) min theta position in Angular Position 72
        #' (2) max theta position in Angular Position 01 (73)
        .df.01 <- temp.df.rs %>% dplyr::filter(., Ang.pos == 1)
        .df.72 <- temp.df.rs %>% dplyr::filter(., Ang.pos == n.theta)
        
        # Padding
        .df.01$theta   <- .df.01$theta + 2 * .pi
        .df.01$Ang.pos <- .df.01$Ang.pos + n.theta
        
        #' [Original df.rs]
        # Section Selection
        .df.bin <- 
          df.rs %>% dplyr::filter(., Ang.pos %in% c(1,2))
        .df.fin <- 
          df.rs %>% dplyr::filter(., Ang.pos %in% c((n.theta-1),n.theta))
        
        # Padding to >2pi region
        .df.bin$theta   <- .df.bin$theta + 2 * .pi
        .df.bin$Ang.pos <- .df.bin$Ang.pos + n.theta
        .df.temp <- rbind(.df.fin, .df.bin)
        
        # rs theta sequence with "ReSample.Size"
        theta.l <- .df.temp$theta %>% min() %>% round(., digits = .digits)
        theta.r <- .df.temp$theta %>% max() %>% round(., digits = .digits)
        seq.rs  <- seq(theta.l, theta.r, by = ReSample.Size)
        npts.rs <- length(seq.rs)
        Ang.prs <- rep(max(.df.temp$Ang.pos), npts.rs)
        
        # Re-mark the Ang.pos
        k  = 1
        nk = 4
        while (k < nk) {
          
          thres.angpos <- (max(.df.temp$Ang.pos) - k) * theta.range / 180 * .pi
          remark.p     <- which(seq.rs <= (thres.angpos + ReSample.Size/10))
          Ang.prs[remark.p] <- (max(.df.temp$Ang.pos) - k)
          
          # loop
          k = k + 1
          
        }
        
        # Smooth.spline
        if(length(unique(.df.temp$theta)) <= 4){
          
          # Linear lm
          lm.rs   <- lm(dist ~ theta, data = df.rs)
          
          # Re-sample the Ring Segment
          dist.pad <- predict(lm.rs, 
                              data.frame(theta = .df.temp$theta)) %>% unname()
          
        }else if(nrow(.df.temp) < 20 |
                 dist(range(.df.temp$theta)) < 0.2){
          
          smsp.pad <- smooth.spline(x  = .df.temp$theta,
                                    y  = .df.temp$dist,
                                    df = 3)
          dist.pad <- predict(smsp.pad, x = seq.rs)$y
          
        }else{
          
          smsp.pad <- smooth.spline(x  = .df.temp$theta,
                                    y  = .df.temp$dist,
                                    df = 5)
          dist.pad <- predict(smsp.pad, x = seq.rs)$y
          
        }
        
        
        # Re-sample the Ring Segment
        temp.df.pad <- data.frame(clst    = Ring.ID[rs],
                                  theta   = seq.rs,
                                  dist    = dist.pad,
                                  Ang.pos = Ang.prs)
        
        # Subset from Original
        .rs01p <- which(temp.df.pad$Ang.pos == (1+n.theta))
        .rs72p <- which(temp.df.pad$Ang.pos == n.theta)
        temp.rs01 <- temp.df.pad[.rs01p,]
        temp.rs72 <- temp.df.pad[.rs72p,]
        
        #' [Make Sure the theta range are the same]
        #'@note 
        #' Most likely, due to theta + 2*pi (where theta = 0)
        #' There will be additional positions produced.
        #'@note 
        #' Use as.character() to avoid the problem
        
        # as.character()
        temp.rs01$theta <- temp.rs01$theta %>% as.character()
        temp.rs72$theta <- temp.rs72$theta %>% as.character()
        .df.01$theta <- .df.01$theta %>% as.character()
        .df.72$theta <- .df.72$theta %>% as.character()
        
        # Temp.Save
        .temp.rs01 <- temp.rs01
        .temp.rs72 <- temp.rs72
        
        # filter
        temp.rs01 <- 
          temp.rs01 %>% dplyr::filter(., theta %in% .df.01$theta)
        temp.rs72 <- 
          temp.rs72 %>% dplyr::filter(., theta %in% .df.72$theta)
        
        # as.numeric()
        temp.rs01$theta <- temp.rs01$theta %>% as.numeric()
        temp.rs72$theta <- temp.rs72$theta %>% as.numeric()
        .df.01$theta <- .df.01$theta %>% as.numeric()
        .df.72$theta <- .df.72$theta %>% as.numeric()
        
        # Plot Check
        if(.plotcheck == TRUE){
          
          plot(x = .df.temp$theta,y = .df.temp$dist)
          lines(smsp.pad, col = "red", lwd = 2)
          points(.df.72$theta, .df.72$dist, col = "blue")
          points(.df.01$theta, .df.01$dist, col = "green")
          
        }
        
        #' [Merging 2 Spline Segments for 01 and 72]
        #' @note only Spline df has 0-2pi theta range
        #' (1) dist.72
        seq.api <- temp.rs72$theta
        D = max(seq.api) - min(seq.api)
        M.dist72 <- 
          .df.72$dist * 
          ((max(seq.api)-.df.72$theta) / D) + 
          temp.rs72$dist * 
          ((temp.rs72$theta - min(seq.api)) / D)
        
        .p72 <- which(.temp.rs72$theta %in% temp.rs72$theta)
        .temp.rs72$dist[.p72] <- M.dist72
        
        #' (2) dist.01 (73)
        seq.api <- temp.rs01$theta
        D = max(seq.api) - min(seq.api)
        M.dist01 <- 
          temp.rs01$dist * 
          ((max(seq.api)-temp.rs01$theta) / D) + 
          .df.01$dist * 
          ((.df.01$theta - min(seq.api)) / D)
        
        .p01 <- which(.temp.rs01$theta %in% temp.rs01$theta)
        .temp.rs01$dist[.p01] <- M.dist01
        
        #' [Insert the section 01 and 72 (n.theta) back]
        # Convert Chr. & padding back
        .temp.rs01$theta <- .temp.rs01$theta %>% as.numeric()
        .temp.rs72$theta <- .temp.rs72$theta %>% as.numeric()
        
        .temp.rs01$theta   <- .temp.rs01$theta   - 2*.pi
        .temp.rs01$Ang.pos <- .temp.rs01$Ang.pos - n.theta
        
        # Remove old ones
        .p01       <- which(temp.df.rs$Ang.pos == 1)
        temp.df.rs <- temp.df.rs[-.p01,]
        .p72       <- which(temp.df.rs$Ang.pos == n.theta)
        temp.df.rs <- temp.df.rs[-.p72,]
        
        # Add the new spline section
        temp.add   <- rbind(.temp.rs01, .temp.rs72)
        # temp.add   <- RingPositionBuild(polar.ring   = temp.add,
        #                                 reference.df = df.clst)
        temp.df.rs <- rbind(temp.df.rs, temp.add)
        
      }
      
    }
    
    # Plot Export Check
    if(.plotexport == TRUE){
      
      temp.plt      <- paste("RingID 0", rs, sep = "")
      temp.plt.name <- paste(temp.plt, ".png", sep = "")
      png(filename = temp.plt.name, width = 3000, height = 1500)
      
      plot(df.rs$theta, df.rs$dist)
      points(temp.df.rs$theta, temp.df.rs$dist, col = "red")
      
      dev.off()
      
    }
    
    # rw.check
    if(rw.check == TRUE){
      
      .rwp <- which(temp.df.rs$Ang.pos > n.theta)
      temp.df.rs$Ang.pos[.rwp] <- 
        temp.df.rs$Ang.pos[.rwp] - n.theta
      temp.df.rs$theta  [.rwp] <- 
        temp.df.rs$theta[.rwp] - 2 * .pi
      
      rw.check = FALSE
      
    }
    
    # rs-Loop
    ring.ReSample[[rs]] <- temp.df.rs
    rs = rs + 1
    print(rs)
    if(rs > nrs){break}
    
  }
  
  ring.ReSample <- dplyr::bind_rows(ring.ReSample)
  if(purrr::is_empty(pith)){
    
    ring.ReSample <- ring.PositionBuild(polar.ring   = ring.ReSample,
                                        reference.df = df.clst)
    
  }else{
    
    ring.ReSample <- ring.PositionBuild(pith = pith,
                                        polar.ring   = ring.ReSample,
                                        reference.df = df.clst)
    
  }
  
  # Ensure the class (int) and the digits (numeric)
  ring.ReSample$theta    <- ring.ReSample$theta   %>% round(., digits = .digits)
  ring.ReSample$dist     <- ring.ReSample$dist    %>% round(., digits = .digits)
  ring.ReSample$x        <- ring.ReSample$x       %>% as.integer()
  ring.ReSample$y        <- ring.ReSample$y       %>% as.integer()
  ring.ReSample$clst     <- ring.ReSample$clst    %>% as.integer()
  ring.ReSample$Ang.pos  <- ring.ReSample$Ang.pos %>% as.integer()
  
  # Output Result
  return(ring.ReSample)
  
}



ring.MatrixBuild <- function(pith = NULL,
                             im.ref  = NULL,
                             ring.df = data.frame(x = vector(),
                                                  y = vector(),
                                                  theta = vector(),
                                                  dist  = vector())){
  
  # -------------------------------------------------------------------------- #
  #' [Input Check]
  if(is.matrix(im.ref)!=TRUE){
    
    if(purrr::is_empty(im.ref)){
      
      stop("Reference Image is Missing")
      
    }else if(imager::is.cimg(im.ref)){
      
      message("CImage is Detected, Convert to Matrix Form")
      im.ref <- im.ref %>% as.matrix()
      
    }else{
      
      stop("Unknown Image Type")
      
    }
    
  }
  
  if(purrr::is_empty(ring.df$x) &
     purrr::is_empty(ring.df$y)){
    
    if(purrr::is_empty(ring.df$theta) &
       purrr::is_empty(ring.df$dist)){
      
      stop("Dataframe 'ring.df' is needed")
      
    }else if(purrr::is_empty(pith)){
      
      stop("Pith is needed when X-Y not exist")
      
    }else{
      
      #' [XY-Coordinates]
      ring.df <-
        ring.df %>% 
        dplyr::mutate(x = pith[1] + round(dist * cos(theta)),
                      y = pith[2] + round(dist * sin(theta)))
      
    }
    
  }
  
  if(purrr::is_empty(ring.df$clst)){ring.df$clst <- 1L}
  
  # -------------------------------------------------------------------------- #
  
  
  
  
  #' [Conversion]
  im    <- im.ref
  im[,] <- 0
  
  .Ring.ID  <- ring.df$clst %>% unique()
  n.Ring.ID <- .Ring.ID %>% length()
  Ring.ID   <- .Ring.ID %>% as.integer() #c(1:n.Ring.ID)
  
  k = 1
  while (k <= n.Ring.ID) {
    
    temp.df <- ring.df %>% dplyr::filter(clst %in% .Ring.ID[k])
    temp.df <- temp.df[, c("x", "y")] %>% as.matrix()
    im[temp.df] <- Ring.ID[k]
    
    # k-loop
    k = k + 1
    
  }
  
  im.matrix <- im
  # im.cimg   <- im %>% imager::as.cimg()
  return(im.matrix)
  
}



ring.WriteID <- function(ring.list = NULL,
                         direction = "inwards"){
  
  message("Writing Ring ID ", direction)
  # -------------------------------------------------------------------------- #
  #' [Input Check]
  if(is.list(ring.list) != TRUE){
    
    stop("A List Ring Data Input Expected")
    
  }
  # -------------------------------------------------------------------------- #
  
  #' [Re-write the Ring.ID]
  k  = 1
  nk = length(ring.list)
  ID.start  = 1L
  Ring.ID <- seq(ID.start, by = 1L, length.out = nk)
  while (k <= nk) {
    
    ring.list[[k]]$clst <- Ring.ID[k] %>% as.character()
    
    # k-loop
    k = k + 1
    
  }
  
  ring.list.ID <- 
    ring.list %>% 
    dplyr::bind_rows() %>% 
    dplyr::group_by(clst) %>% 
    dplyr::summarise(m.dist = mean(dist)) %>% 
    dplyr::arrange(desc(m.dist)) %>% 
    dplyr::mutate(ID   = Ring.ID,
                  clst = clst %>% as.integer()) %>% 
    dplyr::arrange(clst)
  
  #' [Re-write the Ring.ID]
  k  = 1
  nk = length(ring.list)
  Ring.ID <- ring.list.ID$ID
  while (k <= nk) {
    
    ring.list[[k]]$clst <- Ring.ID[k] %>% as.character()
    
    # k-loop
    k = k + 1
    
  }
  
  return(ring.list)
  
}



ring.RelativeDist <- function(ring.TRW,
                              Target.ID,
                              .keep = FALSE){
  
  require(dplyr)
  
  # Check
  if(is.data.frame(ring.TRW) == FALSE){
    
    if(is.list(ring.TRW)){
      ring.TRW <- ring.TRW %>% dplyr::bind_rows()
    }else{
      stop("Error: 'A data.frame' or 'List' is expected")
    }
    
  }
  if((Target.ID %in% ring.TRW$clst) != TRUE){
    stop("Error: Target.ID is wrong")
  }
  
  out.df   <- ring.TRW %>% dplyr::filter(clst == Target.ID)
  Ring.ref <-
    out.df %>% 
    ring.ResampleSegments(df.clst = .,
                          shrink_by = 0.01, # 0.005
                          .plotcheck = FALSE) %>% 
    dplyr::mutate(theta = round(theta, digits = 4))
  # Fix resampleSegment product range (0.0001 - 2*pi)
  .ref <- 
    Ring.ref %>% 
    filter(theta == 2*round(pi, digits = 4)) %>% 
    mutate(theta = 0.0000)
  Ring.ref <- rbind(Ring.ref, .ref)
  temp.df <- 
    ring.TRW %>% 
    dplyr::select("theta") %>% 
    dplyr::mutate(theta = round(theta, digits = 4)) %>% 
    dplyr::left_join(., Ring.ref[, c("theta", "dist")])
  .ring.TRW <-
    ring.TRW %>% 
    dplyr::mutate(dist.ref = temp.df$dist,
                  relative.position = 1 - abs(dist.ref - dist)/dist.ref)
  
  # Output Result
  if(.keep == FALSE){
    ring.TRW$dist <- .ring.TRW$relative.position
  }else{
    ring.TRW <- .ring.TRW
  }
  return(ring.TRW)
  
}

ring.RelativeDist.Mask <- function(ring.TRW,
                                   Mask.Matrix,
                                   pith,
                                   .keep = FALSE){
  
  require(dplyr)
  
  # Self-defined Function
  source("AngularCompute.R")
  source("RingEdge_Functions.R")
  # mask.as.Ring() in Detection Code.R
  
  out.df   <- 
    Mask.Matrix %>% 
    imager::as.cimg() %>%
    mask.as.Ring() %>% 
    edge.FeatureCompute(xc = pith$x, 
                        yc = pith$y) %>% 
    AngularCompute(xc = pith$x,
                   yc = pith$y,
                   df.mt = .,
                   theta.range = 5) %>% 
    purrr::pluck("df.clst")
  
  Ring.ref <-
    out.df %>% 
    ring.ResampleSegments(df.clst = .,
                          shrink_by = 0.01, # 0.005
                          .plotcheck = FALSE,
                          pith = pith) %>% 
    dplyr::mutate(theta = round(theta, digits = 4))
  # Fix resampleSegment product range (0.0001 - 2*pi)
  .ref <- 
    Ring.ref %>% 
    filter(theta == 2*round(pi, digits = 4)) %>% 
    mutate(theta = 0.0000)
  Ring.ref <- rbind(Ring.ref, .ref)
  temp.df <- 
    ring.TRW %>% 
    dplyr::select("theta") %>% 
    dplyr::mutate(theta = round(theta, digits = 4)) %>% 
    dplyr::left_join(., Ring.ref[, c("theta", "dist")])
  .ring.TRW <-
    ring.TRW %>% 
    dplyr::mutate(dist.ref = temp.df$dist,
                  relative.position = 1 - abs(dist.ref - dist)/dist.ref)
  
  # Output Result
  if(.keep == FALSE){
    ring.TRW$dist <- .ring.TRW$relative.position
  }else{
    ring.TRW <- .ring.TRW
  }
  return(ring.TRW)
  
}



mask.as.Ring <- function(mask.cim){
  
  # Data Check
  if(imager::is.pixset(mask.cim) |
     is.matrix(mask.cim)){mask.cim <- imager::as.cimg(mask.cim)}
  
  temp.contour <- 
    mask.cim %>% 
    imager::contours(.,nlevels = 1) %>% 
    as.data.frame() %>% 
    dplyr::mutate(.keep = c("none"),
                  clst = 9999,
                  x = x,
                  y = y)
  
  return(temp.contour)
  
}



mask.join <- function(df.edge, 
                      mask.df,
                      group.dist = 2){
  
  require(dplyr)
  require(purrr)
  
  ta.clst <- mask.df
  df.clst <- df.edge
  
  j = 1
  clst.j  <- vector()
  dist.j  <- vector()
  temp.df <- df.clst[, c("x", "y")] %>% as.matrix()
  while (TRUE) {
    
    temp.ta <- matrix(data = ta.clst[j, c("x", "y")] %>% as.matrix(),
                      nrow = nrow(temp.df),
                      ncol = 2,
                      byrow = TRUE)
    
    temp.dist     <- sqrt(rowSums((temp.df - temp.ta)^2))
    
    .pos <- which(temp.dist < group.dist)
    
    # j-loop updates
    print(j)
    clst.j <- append(clst.j, df.clst$clst[.pos])
    dist.j <- append(dist.j, temp.dist[.pos])
    j = j + 1
    
    if(j > nrow(ta.clst)){break}
    
  }
  
  # Reference Table
  sum.i <- rowsum(dist.j, group = clst.j)
  df.i  <- clst.j %>% table() %>% as.data.frame()
  
  group.clst <- 
    data.frame(clst = as.numeric(levels(df.i$.))[df.i$.],
               dist = sum.i / df.i$Freq)
  
  # New Group Cluster ID
  .pos <- which(df.edge$clst %in% clst.j)
  df.edge$clst[.pos] <- mask.df$clst %>% unique()
  df.edge <- rbind(df.edge, mask.df)
  
  return(df.edge)
  
}
