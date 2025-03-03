





vector.PatternShift <- function(vec){
  
  # (0) Calculation of changes and signs
  sum.vec <- 
    data.frame(vec = sort(vec)) %>% 
    dplyr::transmute(change_vec = c(0, diff(vec)),
                     sign_vec   = sign(change_vec))
  
  # Get the most frequent non-zero sign for x and y
  # (1) Process for assigned vector
  table_vec <- table(sum.vec$sign_vec[sum.vec$sign_vec != 0])
  if(length(table_vec) == 0) {
    vec.sign <- NA  # No changes for such vector
    .det = FALSE
  } else {
    vec.sign <- as.integer(names(table_vec)[which.max(table_vec)])
  }
  
  # (2) Compute the cumulative sum
  cumsum_vec <- cumsum(sum.vec$change_vec)
  
  # Check
  # plot(cumsum_vec)
  
  # Loop to find local-minima
  i = 1
  .step <- FALSE
  local_pick <- vector()
  local_pos  <- vector()
  local_step <- vector()
  temp_pick  <- cumsum_vec[1]
  while (is.na(vec.sign)!=TRUE) {
    
    # if(i == 18){break}
    if(cumsum_vec[i+1]*vec.sign >= temp_pick*vec.sign){
      
      # Step
      if(.step == TRUE){
        
        local_step <- append(local_step, (i - temp_pos))
        .step <- FALSE
        
      }
      
      temp_pick <- cumsum_vec[i+1]
      temp_pos  <- i+1
      
    }else{
      
      local_pos  <- append(local_pos, temp_pos)
      local_pick <- append(local_pick, temp_pick)
      
      # Signal
      .step <- TRUE
      
    }
    
    # i-loop
    i = i + 1
    if(i > (length(cumsum_vec)-2)){
      
      local_pos  <- local_pos %>% unique()
      local_pick <- local_pick %>% unique()
      
      # Step
      if(.step == TRUE){
        
        local_step <- 
          append(local_step, (length(cumsum_vec) - temp_pos))
        .step <- FALSE
        
      }
      
      # .det
      .det <- sum(local_step)
      break
      
    }
    
  }
  
  return(.det)
  
}

RD.xy.det <- function(target.ring) {
  
  require(dplyr)
  require(purrr)
  
  det.x <- target.ring$x %>% vector.PatternShift()
  det.y <- target.ring$y %>% vector.PatternShift()
  
  
  # Determine rotation need based on det.x and det.y
  if(det.x == 0) {
    .det <- FALSE
  } else if(det.y == 0) {
    .det <- TRUE
  }else if(det.x < 0.1 * nrow(target.ring) | 
           det.y < 0.1 * nrow(target.ring) ) {
    
    if(det.x < det.y){
      .det <- FALSE
    }else{
      .det <- TRUE
    }
    
  }else{
    plot(target.ring$x, target.ring$y)
    stop("det.x & det.y Over-size, stopping execution.")
  }
  
  return(.det)
  
}

RingDistance <- function(Target.ring, 
                         Next.ring, 
                         center    = NULL,
                         im.mask   = NULL,
                         view.ang  = 5,
                         spline.df = NULL){
  
  # -------------------------------------------------------------------------- #
  # [Import]
  # Target.ring can be a df system 
  # / or cimg as an option for manual measurement.
  # If we want to display (highlight) the ring
  # use matrix to convert.
  # [See show rings]
  #'@note The manual rings processed by "IJ2R" is in the same system.
  
  # -------------------------------------------------------------------------- #
  #'[XYcoordinates: Ring Pixel Positions]
  if(imager::is.cimg(Target.ring) == TRUE){
    Target.ring <- Target.ring %>% imager::get.locations(., function(v) v>0)
  }
  if(imager::is.cimg(Next.ring) == TRUE){
    Next.ring <- Next.ring %>% imager::get.locations(., function(v) v>0)
  }
  
  # -------------------------------------------------------------------------- #
  #'[Center: Pith]
  # Default center for splitting rings for spline
  # The mean is used,
  # It is possible to have Pith position as input.
  if(purrr::is_empty(center)){
    
    # If the center is a vector (Pith.x, Pith.y)
    # This will not Run
    tr <- Target.ring
    nr <- Next.ring
    center <- c(min(mean(tr$x), mean(nr$X)),
                min(mean(tr$y), mean(nr$Y))) %>% round()
    
  }
  
  # -------------------------------------------------------------------------- #
  #'[Ring Fill for Rings from Fuji]
  # This is done in "Ij2R". No need to combine them.
  
  # -------------------------------------------------------------------------- #
  #'[ImageMasking: The positions to be kept]
  if(purrr::is_empty(im.mask) != TRUE){
    
    # Passive Adjustment
    # To remove the given positions that are breaking-Gap/ Branches
    message("# Acquiring Unwanted Positions from im.mask #")
    
    # Input im.mask
    if(is.matrix(im.mask) | imager::is.cimg(im.mask)){
      
      im.mask <- im.mask %>% imager::as.cimg()
      df.mask <- im.mask %>% imager::get.locations(., function(v) v > 0)
      
    }else{
      
      stop("Only Matrix or Cimage classes are supported")
      
    }
    
    print("# Generate mask position ID")
    ID.mask <- paste(df.mask$x, "_", df.mask$y, sep = "")
    
  }
  
  # -------------------------------------------------------------------------- #
  #'[Data_Preparation:]
  # Check the data structure
  if(is.data.frame(Target.ring)!=TRUE){
    Target.ring <- Target.ring %>% unname() %>% as.data.frame()
  }
  if(is.data.frame(Next.ring)!=TRUE){
    Next.ring <- Next.ring %>% unname() %>% as.data.frame()
  }
  
  # Check the column names
  # especially when they are from ImageJ
  if(all(c("x", "y") %in% names(Target.ring)) != TRUE |
     all(c("x", "y") %in% names(Next.ring)) != TRUE){
    stop("Please Adjust the ring Column name")
  }
  Target.ring <- Target.ring[,c("x", "y")]
  Next.ring   <- Next.ring  [,c("x", "y")]
  
  # Cleaned ring data positions
  tr <- Target.ring %>% edge.FeatureCompute(xc = center[1], yc = center[2])
  nr <- Next.ring   %>% edge.FeatureCompute(xc = center[1], yc = center[2])
  
  # Check
  head(tr)
  head(nr)
  
  # Data preparation
  tr <- tr[order(tr$theta),]
  nr <- nr[order(nr$theta),]
  
  # -------------------------------------------------------------------------- #
  #'[Ring_distance_computation:]
  #'@description 
  #  The ring distance (or differences between Automated and manual ones) is
  #  computed vertically from the "Target" to "Next".
  #'@note 
  #  To take advantage from matrix computation,
  #  The Vertical distance will be acquired from:
  #  2 ring-based spline points for all given points in "Target.ring"
  
  message("# Compute Vertical Distance")
  
  # Acquire data from such gap to perform spline
  # Range of data acquisition:
  degree.range = view.ang
  theta.range  = (degree.range/180)*pi
  cat(paste(paste("# Degree range =", degree.range),
            paste("# Radius range =", theta.range %>% round(., digits = 4)),
            "\n",
            sep = "\n"))
  
  # i-loop controls the sequence of position as input
  i = 1
  n.i = dim(tr)[1]
  Ring.dist <- vector()
  point.cos <- vector()
  while (TRUE) {
    
    # Report the progress
    print(paste("Processing: position ", i, "/", n.i, sep = ""))
    
    # Check
    # if(i == 5292){stop()}
    
    # data acquired in between the angular range
    theta.lb <- tr$theta[i] - theta.range
    theta.rb <- tr$theta[i] + theta.range
    
    # Target.ring
    tr.i <- tr  [tr  $theta > theta.lb,]
    tr.i <- tr.i[tr.i$theta < theta.rb,]
    
    # Next.ring
    nr.i <- nr  [nr  $theta > theta.lb,]
    nr.i <- nr.i[nr.i$theta < theta.rb,]
    
    # Adjustment for theta < theta.range & close to 2pi
    if(theta.lb < 0){
      
      temp.lb <- 2*pi + theta.lb
      temp.tr <- tr[tr$theta > temp.lb,]
      temp.nr <- nr[nr$theta > temp.lb,]
      
      # Add back
      tr.i <- rbind(tr.i, temp.tr)
      nr.i <- rbind(nr.i, temp.nr)
      
    }
    
    if(theta.rb > 2*pi){
      
      temp.rb <- theta.rb - 2*pi
      temp.tr <- tr[tr$theta < temp.rb,]
      temp.nr <- nr[nr$theta < temp.rb,]
      
      # Add back
      tr.i <- rbind(tr.i, temp.tr)
      nr.i <- rbind(nr.i, temp.nr)
      
    }
    
    # Target ring position [i] for distance computation
    tr.pi <- c(tr$x[i], tr$y[i])
    
    # Check the need of exchange the axis
    xy.rotate = FALSE
    # if(length(unique(tr.i$x))*sqrt(3) < length(unique(tr.i$y)) |
    #    length(unique(nr.i$x))*sqrt(3) < length(unique(nr.i$y))){
    if((IQR(tr.i$x) < 1 | IQR(nr.i$x) < 1) &
       (length(unique(tr.i$y)) != 1 & length(unique(nr.i$y)) != 1) |
       RD.xy.det(tr.i) | RD.xy.det(nr.i)){
      
      xy.rotate = TRUE
      temp.trx <- tr.i$x
      temp.try <- tr.i$y
      tr.i$x <- temp.try
      tr.i$y <- temp.trx
      
      temp.nrx <- nr.i$x
      temp.nry <- nr.i$y
      nr.i$x <- temp.nry
      nr.i$y <- temp.nrx
      
      tr.pi <- c(tr.pi[2], tr.pi[1])
      
    }else{
      
      tr.pi <- tr.pi
      
    }
    
    # [Important]
    #'@description 
    # (1) sp.tr: Compute the vertical vector
    # (2) sp.nr: Compute the vertical distance
    .lm.ctrl <- FALSE
    if(length(unique(tr.i$x)) < 4 | length(unique(nr.i$x)) < 4){
      
      # unique value is not able to support smooth.spline fitting
      # use lm instead
      .lm.ctrl <- TRUE
      sp.tr <- lm(y ~ x, data = tr.i[,c("x", "y")])
      sp.nr <- lm(y ~ x, data = nr.i[,c("x", "y")])
      
    }else if(purrr::is_empty(spline.df)){
      
      # smooth.spline df auto (This might cause over fitting)
      # sp.tr <- smooth.spline(x = tr.i$x, y = tr.i$y)
      # sp.nr <- smooth.spline(x = nr.i$x, y = nr.i$y)
      
      if(length(unique(tr.i$x)) < 100 | length(unique(nr.i$x)) < 100){
        
        sp.tr <- smooth.spline(x = tr.i$x, y = tr.i$y, df = 3)
        sp.nr <- smooth.spline(x = nr.i$x, y = nr.i$y, df = 3)
        
      }else{
        
        sp.tr <- smooth.spline(x = tr.i$x, y = tr.i$y, df = 5)
        sp.nr <- smooth.spline(x = nr.i$x, y = nr.i$y, df = 5)
        
      }
      
    }else{
      
      sp.tr <- smooth.spline(x = tr.i$x, y = tr.i$y, df = spline.df)
      sp.nr <- smooth.spline(x = nr.i$x, y = nr.i$y, df = spline.df)
      
    }
    
    # Distance computation loop
    #'@note
    # j-loop is used to control the 10-step search 
    # from 0.1 -> 0.01 -> ..., 
    # until tn.cos reaches the angle of 0.9995 (1 degree)
    j = 1
    .step   = 10^(-j)
    .search = 10*.step
    while (TRUE) {
      
      #'@note
      # Use j = 1 to build (1) tr.vvec (2) nr.pi
      # j-loop to approach the possible "vertical" situation
      if(j == 1){
        
        #'[VerticalVector: tr.vvec]
        #'@note 
        # xy.rotate == FALSE
        # -> input tr.pi "x" value for prediction: tr.pi[1]
        # xy.rotate == TRUE
        # -> input tr.pi "y" value for prediction: tr.pi[2]
        # [Important]
        # Vectors refers to the original tr.pi
        # When j-loop starts, it moves and re-predict the values in nr.fit
        # The tr.fit do not do anything to the distance computation anymore.
        xt.seq  <- seq((tr.pi[1] - .search), (tr.pi[1] + .search), by = .step)
        
        # Prediction control
        if(.lm.ctrl == TRUE){
          tr.fit <- data.frame(x = xt.seq,
                               y = predict(sp.tr, 
                                           newdata = data.frame(x = xt.seq)))
        }else{
          tr.fit <- predict(sp.tr, x = xt.seq) %>% as.data.frame()
        }
        
        
        # Vector computation
        d.tr <- vector()
        d.nr <- vector()
        n = 1
        n.p = dim(tr.fit)[1]
        while (TRUE) {
          
          temp.dtr <- tr.fit$y[n+1] - tr.fit$y[n]
          d.tr     <- append(d.tr, temp.dtr)
          
          # temp.dnr <- nr.fit$y[n+1] - nr.fit$y[n]
          # d.nr     <- append(d.nr, temp.dnr)
          
          # loop
          n = n+1
          if(n == n.p){break}
          
        }
        
        tr.vec <- c(0.1, mean(d.tr))
        # nr.vec <- c(0.1, mean(d.nr))
        
        #'@note 
        # Due to the fact that we just need "absolute" difference
        # between two rings,
        # Use "abs()" value instead of 
        # distinguishing the vector position on the ring.
        # Therefore, 
        # The vertical vectors above can be defined as:
        tr.vvec <- c(mean(d.tr), -0.1)
        # nr.vvec <- c(mean(d.nr), -0.1)
        
        #'[Most_vertical_Position: nr.pi]
        #'@description 
        # The position where to search the intersects on next.ring (j-loop)
        nr.i <- nr.i[, c("x", "y")]
        tr.pim  <- matrix(tr.pi, nrow = nrow(nr.i), ncol = 2, byrow = TRUE)
        tn.dist <- sqrt( rowSums( (nr.i - tr.pim)^2 ) ) %>% as.numeric()
        
        tr.vecm <- matrix(tr.vvec, nrow = nrow(nr.i), ncol = 2, byrow = TRUE)
        tn.vecm <- nr.i - tr.pim
        
        tn.cos  <- 
          as.numeric(rowSums(tr.vecm * tn.vecm)) / 
          ((tn.dist)* sqrt(sum(tr.vvec^2)))
        
        det.cos <- max(abs(tn.cos))
        det.tnp <- which.max(abs(tn.cos))
        nr.pi   <- nr.i[det.tnp,] %>% as.numeric()
        
      }
      
      #'[Same_Position_Check: tr.pi == nr.pi]
      if(purrr::is_empty(which(is.na(tn.cos))) != TRUE){
        
        # Two Points Overlaps
        #'@description 
        # The calculation tolerate this situation,
        # For the purpose of:
        # (1) Error computation
        # (2) Some missing tree rings might happen.
        # Therefore, the loop should break and go to next i
        det.cos  <- 1
        det.tnp  <- which(is.na(tn.cos))
        nr.pi    <- nr.i[det.tnp,] %>% as.numeric()
        res.dist <- tn.dist[det.tnp]
        break
        
      }
      
      # Start Approaching the most Vertical Position:
      xn.seq <- seq((nr.pi[1] - .search), (nr.pi[1] + .search), by = .step)
      if(.lm.ctrl == TRUE){
        nr.fit <- data.frame(x = xn.seq,
                             y = predict(sp.nr, 
                                         newdata = data.frame(x = xn.seq)))
      }else{
        nr.fit <- predict(sp.nr, x = xn.seq) %>% as.data.frame()
      }
      
      
      # [Important]
      # As the distance (TRW) is based on the "target.ring",
      # only "tr.vvec" will be used.
      
      #'[Vertical_Distance: tn.dist]
      tr.pim  <- matrix(tr.pi, nrow = dim(nr.fit)[1], ncol = 2, byrow = TRUE)
      tn.dist <- sqrt( rowSums( (nr.fit - tr.pim)^2 ) )
      
      #'[Supervis_Selected_Pair: det.cos & det.tnp]
      tr.vecm <- matrix(tr.vvec, nrow = dim(nr.fit)[1], ncol = 2, byrow = TRUE)
      tn.vecm <- nr.fit - tr.pim
      tn.cos  <- rowSums(tr.vecm * tn.vecm) / ((tn.dist)* sqrt(sum(tr.vvec^2)))
      
      #'[Fix 2-points overlapped]
      if(any(is.na(tn.cos)) | any(tn.dist == 0)){
        
        det.cos <- 1
        det.tnp <- 0
        
      }else{
        
        det.cos <- max(abs(tn.cos))
        det.tnp <- which.max(abs(tn.cos))
        
      }
      
      if(j > 3){ 
        
        res.dist <- tn.dist[det.tnp]
        break 
        
      }else if(det.cos < 0.9995){
        
        j = j + 1
        .step   = 10^(-j)
        .search = 10*.step
        # Update tr.pi for the search range change in nr
        nr.pi <- nr.fit[det.tnp,] %>% as.numeric()
        
      }else{ 
        
        res.dist <- tn.dist[det.tnp]
        break
        
      }
      
    }
    
    # Check j-loop results
    det.cos
    det.tnp
    res.dist
    nr.pi
    
    # Debug [Not Run]
    # if(res.dist > 50){break}
    
    # Output Results
    Ring.dist <- append(Ring.dist, res.dist)
    point.cos <- append(point.cos, det.cos)
    
    # i-loop control
    i = i+1
    if(i > n.i){break}
    
  }
  
  # Check
  # hist(Ring.dist)
  # hist(point.cos)
  # mean(Ring.dist)
  # mean(point.cos)
  # sd(Ring.dist)
  
  res.ringdist <- tr %>% dplyr::mutate(Ring.dist = Ring.dist,
                                       point.cos = point.cos)
  
  # Return Output
  return(res.ringdist)
  
}