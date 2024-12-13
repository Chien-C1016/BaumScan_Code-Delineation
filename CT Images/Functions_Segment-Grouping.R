# ---------------------------------------------------------------------------- #
# Functions of Segment Grouping #
# ---------------------------------------------------------------------------- #
# General ----
#'@note 
#' Please refers to "Detection Code_(2) Segment Grouping.R".
#' The functions defined here are the supporters for the scripts.
#'@description 
#' Based on characteristics of each selected tree ring segments for 
#' morphological grouping, the reference ring, both edges of each segments 
#' were used to 
#' (1) filter,
#' (2) select,
#' (3) padding,
#' (4) overlap check,
#' (5) perform pre-defined spline function,
#' (6) fill the gap within the morphed product,
#' (7) define the gap size and segment length.
#' Additionally, to select tree ring segments based on their distance to
#' the reference tree ring (to reduce unnecessary time for morph check),
#' (1) a slice function is defined;
#' Lastly,
#' (1)
#' (2)

# ---------------------------------------------------------------------------- #

# Functions ----
## Distance Grouping Among all tree ring segments
#'@description 
#' As all tree ring segments were over-segmented by dbscan
#' due to the previous vector-clean function,
#' This offers a short cut to control the tolerance grouping by 
#' (1) distance to pith (dist), and
#' (2) theta (radius, input in degree for intuitive input consideration)
#' After grouping over-segmented segments, it is often found that both 
#' edges behave strangely. Therefore, 
#' (1) shrink_by is designed for reducing such impacts.
#' For fill-up the gaps within the grouped segments,
#' (1) Morph.method can be decided between "linear" or "spline" 
#'@source
dist.group <- function(ring.data, 
                       tol.degree = 1, 
                       dist.thres = 1.5,
                       Morph.step = 0.0001,
                       shrink_by  = 0,
                       Morph.method = "linear"){
  
  #'[Preparation:]
  # Corresponding Round digits to Morph.step
  .digits <- log10(Morph.step) %>% abs()
  .pi     <- round(pi, digits = .digits)
  .pi2    <- round(2*.pi, digits = .digits)
  
  # Morph.method Check
  if(length(Morph.method) != 1){
    stop("Only ONE input for Morph.method is Allowed")
  }
  if((Morph.method %in% c("linear", "spline")) != TRUE){
    stop("For Morph.method, please choose either linear or spline")
  }
  
  # Radius Conversion & pi reduced by Morph.step
  .pi = round(pi, digits = .digits)
  theta.thres = (tol.degree/360*2*.pi) %>% round(., digits = .digits)
  
  #'[Group_Ring_Segments:]
  #' while-loop
  #'@note 
  #' while-loop is controlled by group.ID.
  #' It is updated through the while loop.
  while (TRUE) {
    
    #'[Update the temp.sum summary]
    temp.sum <-
      ring.data %>% 
      dplyr::group_by(clst) %>% 
      dplyr::summarise(max.theta = max(theta),
                       min.theta = min(theta),
                       npts      = n(),
                       dist.L    = dist[which.min(theta)],
                       dist.R    = dist[which.max(theta)])
    
    #'[Get targeted tree ring segments (group.ID) for dist.group]
    #'@note
    #' i is the i-th tree ring segment
    i = 1
    ni = nrow(temp.sum)
    group.ID <- vector()
    while (i <= ni) {
      
      delta.theta = temp.sum$min.theta[i] - temp.sum$max.theta
      delta.dist  = (temp.sum$dist.L[i] - temp.sum$dist.R) %>% abs()
      
      # theta & distance threshold
      temp.group <- which(abs(delta.theta) < theta.thres & 
                            delta.dist < dist.thres)
      
      # Remove Self-Considering if ANY
      if(purrr::is_empty(temp.group) != TRUE){
        
        temp.group <- dplyr::setdiff(temp.group, i)
        
      }
      
      if(purrr::is_empty(temp.group) != TRUE &
         length(temp.group) > 0){
        
        # Save the group Ring Segment ID
        group.ID <- append(group.ID, temp.sum$clst[i])
        
        # Merge the surrounding segments as ID = i
        #'(1) temp.summary
        temp.group.clst <- temp.sum$clst[temp.group]
        temp.sum$clst[temp.group] <- temp.sum$clst[i]
        #'(2) ring.data
        temp.datafix <- which(ring.data$clst %in% temp.group.clst)
        ring.data$clst[temp.datafix] <- temp.sum$clst[i]
        
      }
      
      #i-loop
      #'@note 
      #' temp.sum is already updated accordingly if group.ID exists.
      print(paste("Process Segment no. ", i, sep = ""))
      i = i + 1
      
    }
    
    # How much ring segments exist after merging? (Only to know)
    temp.sum$clst %>% unique() %>% length()
    
    #'[while_Loop:control]
    if(purrr::is_empty(group.ID)){break}else{print(group.ID)}
    
    
    #'[Grouped_Segments_Refine:] 
    #'@description 
    #' To fill-up the gaps in between each ring segments 
    #' that were grouped into the same segment ID.
    Refine.ID <- 
      ring.data %>% 
      dplyr::filter(., clst %in% group.ID) %>% 
      dplyr::select(clst) %>% 
      unique() %>% 
      unlist() %>% 
      as.integer()
    
    # i-loop
    #' i controls which segment ID to be considered to be "smoothed",
    #' by then, the little gaps in between will be filled up.
    i = 1
    ni = length(Refine.ID)
    list.i <- list()
    while (i <= ni) {
      
      # Check
      # if(Refine.ID[i] == 1322){stop()}
      
      temp.i <- 
        ring.data %>% 
        dplyr::filter(., clst %in% Refine.ID[i]) %>% 
        dplyr::arrange(theta)
      
      # plot(temp.i$theta, temp.i$dist)
      
      # Inner (within cluster) gaps more than Morph.step
      delta.theta <- temp.i$theta %>% diff()
      
      # ---------------------------------------------------------------------- #
      #'[Ensure Gap Size]
      #' Ensure the Gap Positions and *Size*
      #' There is no padding control here, the one needs to be padded 
      #' will not find close-by neighbors.
      #' However, it will be grouped by other near-by segments...
      gap.pos <- which(delta.theta > theta.thres)
      
      if(purrr::is_empty(gap.pos) != TRUE){
        
        #'@note
        #' Huge Gap exists in the ring segment structure.
        #' Hence, Padding is needed.
        # Signal
        segment.padding <- TRUE
        
        # Check
        if(length(gap.pos) != 1){
          
          plot(temp.i$theta, temp.i$dist)
          stop("Special scenario, Check segment structures")
          
        }
        
        # Clip the segment structure
        temp.i.r <- temp.i[c((gap.pos+1):nrow(temp.i)),]
        temp.i.l <- temp.i[c(1:gap.pos),]
        
        # Padding...
        temp.i.l$theta <- (temp.i.l$theta + .pi2) %>% round(digits = .digits)
        temp.i <- rbind(temp.i.r, temp.i.l)
        
        # Update gap.pos
        delta.theta <- temp.i$theta %>% diff()
        gap.pos     <- which(delta.theta > Morph.step*1.1)
        
      }else{
        
        # Signal
        segment.padding <- FALSE
        
        # Update gap.pos
        gap.pos <- which(delta.theta > Morph.step*1.1)
        
      }
      # ---------------------------------------------------------------------- #
      
      #'[Morph the Gaps in between]
      temp.k <- temp.i[vector(),]
      if(Morph.method == "linear"){
        
        if(purrr::is_empty(gap.pos)!=TRUE){
          
          k = 1
          nk = length(gap.pos)
          list.k <- list()
          while (k <= nk) {
            
            temp.theta <- seq(from = temp.i$theta[gap.pos[k]], 
                              to   = temp.i$theta[gap.pos[k]+1],
                              by   = Morph.step)
            temp.dist  <- seq(from = temp.i$dist[gap.pos[k]],
                              to   = temp.i$dist[gap.pos[k]+1],
                              length.out = length(temp.theta))
            temp.k <- data.frame(clst  = temp.i$clst[1],
                                 theta = temp.theta,
                                 dist  = temp.dist)
            temp.k <- temp.k[-c(1,nrow(temp.k)),]
            
            # k-loop
            # ring.data <- rbind(ring.data, temp.k) ## Slow
            list.k[[k]] <- temp.k
            k = k + 1
            
          }
          
          temp.k <- dplyr::bind_rows(list.k)
          
        }
        
      }else if(Morph.method == "spline"){
        
        
        
        temp.x    <- seq(min(temp.i$theta),
                         max(temp.i$theta),
                         by = Morph.step)
        temp.range  <- temp.x %>% range() %>% diff()
        # if(temp.range < 1.5){temp.df  = 10}else{temp.df = 40}
        # temp.smsp <- smooth.spline(x  = temp.i$theta,
        #                            y  = temp.i$dist,
        #                            df = temp.df)
        
        # Shrink fun()
        .shrink_limit <- 0.1
        if((2*shrink_by / temp.range) < .shrink_limit){
          
          if(shrink_by > 1){shrink_by <- round(shrink_by)*Morph.step}
          temp.x  <- seq((min(temp.i$theta) + shrink_by), 
                         (max(temp.i$theta) - shrink_by), 
                         by = Morph.step)
          
        }else{
          
          .shrink_by <- (temp.range * .shrink_limit) %>% round(digits = 4)
          temp.x  <- seq((min(temp.i$theta) + .shrink_by), 
                         (max(temp.i$theta) - .shrink_by), 
                         by = Morph.step)
          
        }
        
        if(temp.range < 2*(5/180*pi)){
          
          temp.smsp <- smooth.spline(x  = temp.i$theta,
                                     y  = temp.i$dist,
                                     df = 3)
          
        }else if(temp.range < 0.7){
          
          temp.smsp <- smooth.spline(x  = temp.i$theta,
                                     y  = temp.i$dist,
                                     df = 5)
          
        }else{
          
          temp.smsp <- smooth.spline(x  = temp.i$theta,
                                     y  = temp.i$dist,
                                     df = 15)
          
        }
        
        temp.sm   <- predict(temp.smsp, x = temp.x) %>% as.data.frame()
        temp.k    <- data.frame(clst  = temp.i$clst[1],
                                theta = temp.sm$x,
                                dist  = temp.sm$y)
        
      }
      
      # Padding Adjustments
      if(segment.padding == TRUE){
        
        # Signal-off
        segment.padding <- FALSE
        
        # Padding Back
        .pos <- which(temp.k$theta > .pi2)
        temp.k$theta[.pos] <- temp.k$theta[.pos] - .pi2
        temp.k$theta       <- temp.k$theta %>% round(digits = .digits)
        
      }
      
      # i-loop
      print(paste("Process Segment no. ", i, sep = ""))
      list.i[[i]] <- temp.k
      i = i + 1
      
    }
    
    # Remove the old ring Segments in Refine.ID
    # If method is spline
    if(Morph.method == "spline"){
      
      keep.ID   <- dplyr::setdiff(ring.data$clst, Refine.ID)
      ring.data <- ring.data %>% dplyr::filter(., clst %in% keep.ID)
      
    }
    
    temp.list.i <- dplyr::bind_rows(list.i)
    ring.data   <- rbind(ring.data, temp.list.i)
    
  }
  
  return(ring.data)
  
}



## Join 2 ring segments by "theta" value
#'@description 
#' This is to add reference ring to the "ring.ck" for dist comparison.
#' As the theta value will be unified for all ring segments,
#' It is possible to get differences between any segment and "ring.segment".
#' For example,
#' the dist value differences between 
#' ring.ck & reference ring (ring.segment)
#'@import 
#' (1) ring.ck: the target ring segment(s) to be compared.
#' (2) ring.segment: the reference ring segments.
#'@source
ring.join <- function(ring.ck, ring.segment, .digits = 4L){
  
  temp.ck <- ring.ck
  temp.ck$theta <- 
    temp.ck$theta %>% 
    round(digits = .digits) %>% 
    format(scientific = FALSE) %>% 
    as.character() %>% 
    gsub(" ", "", .)
  
  k = 1
  ring.ID = ring.segment$clst %>% unique()
  nk = ring.ID %>% length()
  while (k <= nk) {
    
    ring.segment.k <- 
      ring.segment %>% dplyr::filter(., clst == ring.ID[k])
    dist.name <- 
      ring.segment.k$clst %>% 
      unique() %>% 
      paste("dist.", ., sep = "")
    temp.rs <- ring.segment.k[,c("theta", "dist")]
    names(temp.rs) <- c("theta", dist.name)
    temp.rs$theta  <- 
      temp.rs$theta %>% 
      round(digits = .digits) %>% 
      format(scientific = FALSE) %>% 
      as.character() %>% 
      gsub(" ", "", .)
    
    temp.ck <-
      temp.ck %>% 
      dplyr::left_join(., temp.rs, by = "theta")
    
    k = k + 1
    
  }
  
  
  temp.ck$theta <- temp.ck$theta %>% as.numeric()
  return(temp.ck)
  
}



## Ring Substitute
#'@description 
#' This function is only used to get tree ring segments that are closed
#' to the morphed tree ring into the morphed ring product.
#'@note
#' This is ONLY used in the "Detection Code 02", updating "ring.morph"
#'@source
ring.substitute <- function(ring.1, ring.segment, .digits = 4L){
  
  require(dplyr)
  
  temp.r1 <- ring.1
  temp.r1$theta  <- 
    temp.r1$theta %>% 
    round(digits = .digits) %>% 
    format(scientific = FALSE) %>% 
    as.character() %>% 
    gsub(" ", "", .)
  
  temp.rs <- ring.segment[, c("theta", "dist")]
  temp.rs$theta  <- 
    temp.rs$theta %>% 
    round(digits = .digits) %>% 
    format(scientific = FALSE) %>% 
    as.character() %>% 
    gsub(" ", "", .)
  names(temp.rs) <- c("theta", "rs.dist")
  
  temp.r1 <- 
    temp.r1 %>% 
    dplyr::left_join(., temp.rs[c("theta", "rs.dist")]) %>% 
    dplyr::filter(., is.na(rs.dist) == TRUE)
  temp.r1$rs.dist <- NULL
  
  # Ensure the theta data structure
  temp.r1$theta <- 
    temp.r1$theta %>%
    as.numeric() %>% 
    round(., digits = .digits) %>% 
    as.character() %>% as.numeric()
  temp.r1 <- rbind(temp.r1, ring.segment[,c("clst","theta", "dist")])
  
  return(temp.r1)
  
}



## Ring Structure Spline smoothing
#'@description 
#' This is to automate the smoothing spline process,
#' This helps to reduce the jagged boundaries of the tree ring structure.
#' Hence, help to generalized the "reference ring".
#'@import 
#' (1) ring.base: Target Ring Structure to be smoothed.
#' (2) window.size: Smoothing by slide-Window method.
#' (3) Morph.step: How fine to re-sample the "ring.base".
#' (4) .plotcheck: Whether to plot the splined result.
#'@source
ring.spline <- function(ring.base, 
                        window.size = 12, 
                        Morph.step  = 0.0001,
                        .plotcheck  = FALSE){
  
  # Corresponding Round digits to Morph.step
  .digits = log10(Morph.step) %>% abs()
  
  # Radius conversion
  .pi <- round(pi, digits = .digits)
  window.size <- (window.size/180*.pi) %>% round(., digits = .digits)
  
  # Always take the center part of the spline product
  #'@note keep 1st of 1st and 3rd of the end, (<= + <)
  
  k = 1
  start = Morph.step # min(ring.base$theta)
  ring.sp.list <- list()
  .k.control   <- FALSE
  while (TRUE) {
    
    temp.lb <- start   + window.size * (k-1)
    temp.rb <- temp.lb + window.size * 3 - Morph.step
    if((max(ring.base$theta) - temp.rb) < window.size){
      temp.rb    <- max(ring.base$theta)
      .k.control <- TRUE
    }
    
    temp.df <- 
      data.frame(theta = seq(temp.lb, temp.rb, by = Morph.step)) %>% 
      ring.join(., ring.base)
    
    # Spline
    temp.sm <- smooth.spline(x  = temp.df$theta,
                             y  = temp.df$dist.ref,
                             df = 15)
    
    # Start Keeping
    # ex. 123 | 456 | 789
    nrow(temp.df)
    
    #'[Spline-Predicted Sections]
    temp.sec1 <-  
      predict(temp.sm,
              x = seq(temp.lb, 
                      temp.lb + window.size*1 - Morph.step, 
                      by = Morph.step)) %>% 
      as.data.frame()
    
    names(temp.sec1) <- c("theta", "dist")
    temp.sec1$theta  <- temp.sec1$theta %>% round(digits = .digits)
    
    temp.sec3 <-  
      predict(temp.sm,
              x = seq(temp.lb + window.size * 2, 
                      temp.rb, 
                      by = Morph.step)) %>% 
      as.data.frame()
    
    names(temp.sec3) <- c("theta", "dist")
    temp.sec3$theta  <- temp.sec3$theta %>% round(digits = .digits)
    
    temp.sec2 <-  
      predict(temp.sm,
              x = seq(temp.lb + window.size * 1, 
                      temp.lb + window.size * 2 - Morph.step, 
                      by = Morph.step)) %>% 
      as.data.frame()
    
    names(temp.sec2) <- c("theta", "dist")
    temp.sec2$theta  <- temp.sec2$theta %>% round(digits = .digits)
    
    #'[Ensure data structure]
    temp.sec1$theta <- 
      temp.sec1$theta %>% 
      round(., digits = .digits) %>% 
      as.character() %>% 
      as.numeric()
    
    temp.sec2$theta <- 
      temp.sec2$theta %>% 
      round(., digits = .digits) %>% 
      as.character() %>% 
      as.numeric()
    
    temp.sec3$theta <- 
      temp.sec3$theta %>% 
      round(., digits = .digits) %>% 
      as.character() %>% 
      as.numeric()
    
    
    #'[Save & Merge the Spline Sections]
    ring.sp.list[[k+1]] <-temp.sec2
    
    if(k == 1){
      
      ring.sp.list[[1]]   <-temp.sec1
      
    }else{
      
      prev.sec2   <- ring.sp.list[[k]]
      r.merge <- 
        prev.sec2 * seq(1, 0,length.out = nrow(prev.sec2)) +
        temp.sec1 * seq(0, 1,length.out = nrow(prev.sec2))
      
      # Ensure data structure
      r.merge$theta <- 
        r.merge$theta %>% 
        round(., digits = .digits) %>% 
        as.character() %>% 
        as.numeric()
      
      ring.sp.list[[k]] <- r.merge
      
    }
    
    if(temp.rb == max(ring.base$theta)){
      
      ring.sp.list[[k+2]] <-temp.sec3
      
    }
    
    
    
    # k-loop control
    k = k + 1
    if(.k.control == TRUE){break}
    
  }
  
  #'[Fix 0-2pi Position]
  #'@note pad by 4 k sections
  
  #'[Preparation]
  temp.dfL <- dplyr::bind_rows(ring.sp.list[c(1:2)])
  temp.dfR <- dplyr::bind_rows(ring.sp.list[c(k:(k+1))])
  
  temp.dfL <- data.frame(theta = temp.dfL$theta) %>% ring.join(., ring.base)
  temp.dfR <- data.frame(theta = temp.dfR$theta) %>% ring.join(., ring.base)
  
  temp.pad <- temp.dfL
  temp.pad$theta <- temp.pad$theta + 2*.pi
  
  temp.df  <- rbind(temp.dfR, temp.pad)
  temp.sm  <- smooth.spline(x  = temp.df$theta,
                            y  = temp.df$dist.ref,
                            df = 15)
  temp.pad <- 
    predict(temp.sm,
            x = seq(min(ring.sp.list[[k+1]]$theta),
                    max(ring.sp.list[[1]]$theta) + 2*.pi,
                    by = Morph.step)) %>% 
    as.data.frame()
  names(temp.pad) <- c("theta", "dist")
  temp.pad$theta  <- temp.pad$theta %>% round(digits = .digits)
  
  .convp <- which(temp.pad$theta > (2*.pi+Morph.step/10))
  temp.pad$theta[.convp] <- temp.pad$theta[.convp] - 2*.pi
  
  # Ensure data structure
  temp.pad$theta <- 
    temp.pad$theta %>% 
    round(., digits = .digits) %>% 
    as.character() %>% 
    as.numeric()
  
  ring.sp.list[[1]]   <- temp.pad[ .convp,]
  ring.sp.list[[k+1]] <- temp.pad[-.convp,]
  
  #'[Compile & Output]
  ring.sp <- dplyr::bind_rows(ring.sp.list)
  
  # Ensure data structure
  ring.sp$theta <-
    ring.sp$theta %>% 
    round(., digits = .digits) %>% 
    as.character() %>% 
    as.numeric()
  
  ring.sp$clst <- ring.base$clst %>% unique()
  ring.sp <- ring.sp[, c("clst", "theta", "dist")]
  
  # Plot Check
  if(.plotcheck == TRUE){
    
    plot  (ring.base$theta, ring.base$dist, cex = 0.3)
    points(ring.sp$theta,   ring.sp$dist,   cex = 0.3, col = "red3")
    
  }
  
  return(ring.sp)
  
}



## Adjust distance-2-pith-reversed tree ring objects
#'@description 
#' To compute Morph & Search based on the manually delineated 
#' most-outer tree ring, the function is created to simply code structures &
#' allow pipe writing style.
#'@source
ring.BaseAdjust <- function(ring.b.list, reverse.base){
  
  .ring.b.list <-
    lapply(ring.b.list, 
           function(v){
             v$dist <- reverse.base - v$dist
             return(v)
           })
  return(.ring.b.list)
  
}



## Morph/ build tree ring segment at specific section (Theta range)
#'@description 
#' Based on the reference tree ring structure, 
#' this function provide morphological structure for 2 given edges 
#' (with theta & dist). 
#'@source
Section.Morph <- function(Edge.01, Edge.02,
                          dist.01, dist.02,
                          ref.ring,
                          Morph.step = 0.0001){
  
  #'[Preparation:]
  # Corresponding Round digits to Morph.step
  .digits = log10(Morph.step) %>% abs()
  
  #'[Gap_Size:]
  if((Edge.02 - Edge.01) < Morph.step*1.01){
    
    ring.morph <- ref.ring[vector(),]
    
  }else{
    
    .step <- Morph.step*0.1
    temp.ring2 <- ref.ring   %>% dplyr::filter(theta > (Edge.01 + .step))
    temp.ring2 <- temp.ring2 %>% dplyr::filter(theta < (Edge.02 - .step))
    temp.ring2 <- temp.ring2[order(temp.ring2$theta),]
    
    dist.2.L <- temp.ring2$dist[which.min(temp.ring2$theta)]
    dist.2.R <- temp.ring2$dist[which.max(temp.ring2$theta)]
    
    d.shift.L <- dist.01 - dist.2.L
    d.shift.R <- dist.02 - dist.2.R
    
    temp.23.L <- temp.ring2$dist + d.shift.L
    temp.23.R <- temp.ring2$dist + d.shift.R
    
    seq.x <- seq(1, 0, length.out = nrow(temp.ring2))
    
    morph.dist <- 
      temp.23.L * seq.x + 
      temp.23.R * (1 - seq.x)
    
    ring.morph <- data.frame(clst  = "Morph",
                             theta = temp.ring2$theta,
                             dist  = morph.dist)
    
    # Ensure the ring.2 theta structure
    ring.morph$theta <- 
      ring.morph$theta %>% 
      round(., digits = .digits)
    
  }
  
  return(ring.morph)
  
}



## Slice tree ring segments based on "Full-tree-ring" cluster
#'@description 
#' This is specifically used to slice the tree ring segments based on
#' certain tree ring structure.
#' The slice methods are:
#' (1) inner
#' (2) outer
#' With these options, the function can slice tree ring segments in between
#' 2 known "Full-tree-ring" cluster.
#'@source
segment.slice <- function(df.clst, 
                          cluster   = NULL,
                          slice.dir = "inner"){
  
  if(purrr::is_empty(cluster)){
    stop("Please Provide a defined Tree Ring ID to Slice the data")
  }
  
  ring.clst <- df.clst %>% dplyr::filter(., clst %in% cluster)
  ring.clst$clst <- "ref"
  
  df.sum <-
    df.clst %>% 
    ring.join(., ring.segment = ring.clst) %>% 
    dplyr::mutate(delta = dist - dist.ref) %>% 
    dplyr::group_by(clst) %>% 
    dplyr::summarise(npts  = n(),
                     delta = mean(delta)) %>% 
    dplyr::filter(delta > 0)
  
  if(slice.dir == "inner"){
    
    res.clst <- dplyr::setdiff(df.clst$clst, df.sum$clst)
    res.clst <- dplyr::setdiff(res.clst, cluster)
    
  }else if(slice.dir == "outer"){
    
    res.clst <- df.sum$clst
    
  }else{
    
    stop("Undefined Slice Method, please ONLY use 'inner' or 'outer'")
    
  }
  
  df.res   <- df.clst %>% dplyr::filter(clst %in% res.clst)
  return(df.res)
  
}



## Clean (Filter) tree ring segment summary
#'@description 
#' The function is used for filtering the segment summary table
#' for choosing the segments that approximately fit the morphed ring segment. 
#'@import 
#' (1) segment.size: Give tolerance for segments with certain segment size.
#' (2) dist.thres: the average difference between segments and morphed ring.
#'@source
segment.clean <- function(ring.summary, 
                          segment.size = NULL, 
                          dist.thres   = 5){
  
  ring.summary <- ring.summary %>% dplyr::arrange(min.theta)
  
  if(purrr::is_empty(segment.size)){
    
    ring.summary <- 
      ring.summary %>% 
      dplyr::mutate(., dist.1sd = m.delta + sd.delta) %>% 
      dplyr::filter(., dist.1sd < dist.thres) %>% 
      dplyr::mutate(., dist.1sd = NULL) %>% 
      dplyr::arrange(min.theta)
    
  }else{
    
    temp.keep    <- ring.summary %>% dplyr::filter(., segment.l > segment.size)
    ring.summary <- 
      ring.summary %>% 
      dplyr::filter(., segment.l < segment.size) %>% 
      dplyr::mutate(., dist.1sd = m.delta + sd.delta) %>% 
      dplyr::filter(., dist.1sd < dist.thres) %>% 
      dplyr::mutate(., dist.1sd = NULL) %>% 
      rbind(., temp.keep) %>% 
      dplyr::arrange(min.theta)
    
  }
  
  return(ring.summary)
  
}



## Filter tree ring segments within theta range
#'@description 
#' The function is used for filtering the segment data,
#' keeping only the ones within a targeted theta range (Edge01, Edge02)
#'@import 
#' (1) Edge.01, Edge.02: Targeted theta range (Edge01, Edge02)
#' (2) ring.segments: all tree ring segments.
#'@source
segment.filter <- function(Edge.01, Edge.02,
                           ring.segments,
                           Morph.step = 0.0001){
  
  #'[Preparation:]
  .step <- Morph.step*0.1
  
  #'[Segment_Filtering:]
  ring.section <- 
    ring.segments %>% 
    dplyr::filter(theta > (Edge.01 + .step)) %>% 
    dplyr::filter(theta < (Edge.02 - .step))
  
  #'[Filter from left (edge close to 0)]
  rest.clst.L <- 
    ring.segments %>% 
    dplyr::filter(theta < (Edge.01 + .step)) %>% 
    dplyr::select(clst) %>% 
    unique() %>% 
    unlist() %>% 
    unname()
  
  #'[Filter from right (edge close to 2pi)]
  rest.clst.R <- 
    ring.segments %>% 
    dplyr::filter(theta > (Edge.02 - .step)) %>% 
    dplyr::select(clst) %>% 
    unique() %>% 
    unlist() %>% 
    unname()
  
  #'[Target Clusters = All - (Clusters_L + Clusters_R)]
  rest.clst   <- append(rest.clst.L, rest.clst.R)
  select.clst <- ring.section$clst %>% unique()
  select.clst <- dplyr::setdiff(select.clst, rest.clst)
  
  ring.segments <- 
    ring.section %>% 
    dplyr::filter(., clst %in% select.clst)
  
  ring.segments <- 
    ring.segments[, c("clst", "theta", "dist")]
  return(ring.segments)
  
}



## Select tree ring segments with certain requirements
#'@description 
#' The function is used for filtering the segment data summary,
#' furthermore, supporting the neighboring search for 
#' (1) small neighboring segments, 
#' (2) Middle segments,
#' (3) small segments between found segments and its neighboring edge
#'@import 
#' (1) ring.summary
#' (2) small.segment.search
#' (3) neighbor.thres (*): close-neighboring threshold,
#' (4) rev: to search from Edge02 (Edge close to 2pi)
#' (5) dist.check: TRUE or FALSE to activate,
#' (6) dist.thres: 5 pixels.
#'@source
segment.select <- function(ring.summary,
                           # search.range,
                           small.segment.search = 2 * pi / 10,
                           neighbor.thres = small.segment.search / 10,
                           rev = FALSE,
                           dist.check = FALSE,
                           dist.thres = 5){
  
  if(purrr::is_empty(ring.summary$clst)!=TRUE){
    
    if(nrow(ring.summary) >= 2){
      
      #'@description Method 2
      #' The Modified code use only "neighbor.thres"
      #' to identified how many close-neighboring tree ring segments 
      #' should be taken for smooth.spline computation. 
      
      #'@description (rev option)
      #' In some cases, the trust segments are from behind (larger theta).
      #' Hence, to provide possibility to handle the segment group
      #' from such direction, rev option is developed.
      
      if(rev == TRUE){
        
        ring.summary <- ring.summary[order(-ring.summary$max.theta),]
        
        # rev adjustments:
        temp.min  <- ring.summary$min.theta[-nrow(ring.summary)]
        temp.max  <- ring.summary$max.theta[-1]
        temp.diff <- temp.min - temp.max
        
      }else{
        
        ring.summary <- ring.summary[order( ring.summary$max.theta),]
        
        temp.min  <- ring.summary$min.theta[-1]
        temp.max  <- ring.summary$max.theta[-nrow(ring.summary)]
        temp.diff <- temp.min - temp.max
        
      }
      
      # determinate
      temp.det  <- which(temp.diff < neighbor.thres)
      
      k  = 1
      if(purrr::is_empty(temp.det)!=TRUE){
        
        k  = 1
        nk = length(temp.det)
        while (temp.det[k] == k) {
          
          k = k + 1
          
          if(k > nk){
            break
          }
          
        }
        
      }
      
      temp.search.pos <- k
      ring.summary <- ring.summary[c(1:temp.search.pos),]
      
    }
  }
  
  if(dist.check == TRUE){
    
    if(rev == TRUE){
      
      ring.dist <- 
        abs(ring.summary$dist.r - dplyr::lag(ring.summary$dist.l)) %>% 
        round()
      
    }else{
      
      ring.dist <- 
        abs(ring.summary$dist.l - dplyr::lag(ring.summary$dist.r)) %>% 
        round()
      
    }
    
    # if(purrr::is_empty(dist.thres)){dist.thres <- 5}
    .pos <- which(ring.dist >= dist.thres)
    
    if(purrr::is_empty(.pos)!=TRUE){
      ring.summary <- ring.summary[c(1:(.pos[1]-1)),]
    }
    
  }
  
  return(ring.summary)
  
}



## Padding tree ring segment from Edge02
#'@description 
#' The function is used for padding segments that are not across 0 & 2pi.
#' It is used for padding tree ring segments & reference tree ring.
#'@import 
#' (1) ring.segments: targeted tree ring segments
#' (2) Padding.from: From where to be padded.
#'@source
segment.padding <- function(ring.segments, 
                            Padding.from,
                            Morph.step = 0.0001){
  
  #'[Preparation:]
  # Corresponding Round digits to Morph.step
  .digits = log10(Morph.step) %>% abs()
  
  # Keep 10 points from Padding position (Padding.from)
  Padding.pos = Padding.from + 10 * Morph.step
  
  # Generalize the pi
  .pi = round(pi, digits = .digits) %>% as.character() %>% as.numeric()
  
  #'[Padding:]
  temp.ck.L <- ring.segments %>% dplyr::filter(.,theta <= Padding.pos)
  temp.ck.R <- ring.segments %>% dplyr::filter(.,theta >  Padding.pos)
  
  temp.ck.L$theta <- temp.ck.L$theta + 2*.pi
  
  temp.ring <- rbind(temp.ck.L, temp.ck.R)
  temp.ring <- temp.ring[order(temp.ring$theta),]
  
  return(temp.ring)
  
}



## Overlap-check for tree ring segments listed in ring.summary
#'@description 
#' The function is used for Overlap-check, 
#' designing for tree ring segments listed in ring.summary
#'@import 
#' (1) segment.size: Influence the code functioning, please check @note below.
#' (2) overlap.check: only check the 1st one or All
#'@source
segment.overlap <- function(ring.summary, 
                            segment.size = NULL,
                            overlap.check = c("1st")){
  
  temp.sum <- ring.summary
  
  if(purrr::is_empty(temp.sum$clst)!=TRUE){
    
    if(overlap.check == "1st"){
      
      #'@note
      #' When segment.size exist,
      #' the overlap.check only checks big.segments
      if(purrr::is_empty(segment.size)!=TRUE){
        temp.sum <- 
          ring.summary %>% 
          dplyr::filter(., segment.l > segment.size) %>% 
          dplyr::arrange(min.theta)
      }else{
        temp.sum <- 
          ring.summary %>% 
          dplyr::arrange(min.theta)
      }
      
      temp.keep <- temp.sum[1,]
      temp.sum <-
        temp.sum[-1,] %>% 
        dplyr::mutate(overlap = (min.theta >= temp.keep$min.theta & 
                                   min.theta <= temp.keep$max.theta))
      temp.overlap <- 
        temp.sum %>% 
        dplyr::filter(., overlap == TRUE)
      temp.overlap$overlap <- NULL
      temp.overlap <- rbind(temp.keep, temp.overlap)
      
      if(nrow(temp.overlap) != 1){
        
        temp.remove  <- temp.overlap[-which.min(temp.overlap$m.delta),]
        keep.clst    <- dplyr::setdiff(ring.summary$clst, temp.remove$clst)
        ring.summary <- ring.summary %>% dplyr::filter(., clst %in% keep.clst)
        
      }
      
    }else if(overlap.check == "All"){
      
      #'@note
      #' When segment.size exist,
      #' the overlap.check only ensures 
      #' no small segments overlaps big.segments
      if(purrr::is_empty(segment.size)!=TRUE){
        
        temp.big <- 
          ring.summary %>% 
          dplyr::filter(., segment.l >= segment.size) %>% 
          dplyr::arrange(min.theta)
        
        # Check temp.big first
        temp.sum <- temp.big
        
        k = 1
        while (TRUE) {
          
          temp.keep <- temp.sum[k,]
          temp.sum.k  <-
            temp.sum[-k,] %>% 
            dplyr::mutate(overlap = (min.theta >= temp.keep$min.theta & 
                                       min.theta <= temp.keep$max.theta))
          temp.overlap <- 
            temp.sum.k %>% 
            dplyr::filter(., overlap == TRUE)
          temp.overlap$overlap <- NULL
          temp.overlap <- rbind(temp.keep, temp.overlap)
          
          if(nrow(temp.overlap) != 1){
            
            temp.remove <- temp.overlap[-which.max(temp.overlap$segment.l),]
            keep.clst   <- dplyr::setdiff(temp.big$clst, temp.remove$clst)
            temp.big <- 
              temp.big %>% 
              dplyr::filter(., clst %in% keep.clst)
            
            # Update temp.sum
            temp.sum <- temp.sum %>% 
              dplyr::filter(., clst %in% keep.clst)
            k = 1
            
          }else{
            
            # loop-control
            k = k + 1
            
          }
          
          # Loop-control
          if(k > nrow(temp.sum)){break}
          
        }
        
        temp.small <-
          ring.summary %>% 
          dplyr::filter(., segment.l < segment.size) %>% 
          dplyr::arrange(min.theta)
        
        k = 1
        while (TRUE) {
          
          temp.keep <- temp.big[k,]
          temp.sum.k  <-
            temp.small %>% 
            dplyr::mutate(overlap = (min.theta >= temp.keep$min.theta & 
                                       min.theta <= temp.keep$max.theta))
          
          temp.overlap <- 
            temp.sum.k %>% 
            dplyr::filter(., overlap == TRUE)
          
          if(purrr::is_empty(temp.overlap$clst)!=TRUE){
            
            keep.clst   <- dplyr::setdiff(temp.small$clst, temp.overlap$clst)
            temp.small  <- 
              temp.small %>% 
              dplyr::filter(., clst %in% keep.clst)
            
          }
          
          # Loop-control
          k = k + 1
          if(k > nrow(temp.sum)){break}
          
        }
        
        ring.summary <- 
          rbind(temp.big, temp.small) %>% 
          dplyr::arrange(min.theta)
        
      }else{
        
        temp.sum <- 
          ring.summary %>% 
          dplyr::arrange(min.theta)
        
        k = 1
        while (TRUE) {
          
          temp.keep <- temp.sum[k,]
          temp.sum.k  <-
            temp.sum[-k,] %>% 
            dplyr::mutate(overlap = (min.theta >= temp.keep$min.theta & 
                                       min.theta <= temp.keep$max.theta))
          temp.overlap <- 
            temp.sum.k %>% 
            dplyr::filter(., overlap == TRUE)
          temp.overlap$overlap <- NULL
          temp.overlap <- rbind(temp.keep, temp.overlap)
          
          if(nrow(temp.overlap) != 1){
            
            temp.remove  <- temp.overlap[-which.min(temp.overlap$m.delta),]
            keep.clst    <- dplyr::setdiff(ring.summary$clst, temp.remove$clst)
            ring.summary <- 
              ring.summary %>% 
              dplyr::filter(., clst %in% keep.clst)
            
            # Update temp.sum
            temp.sum <- temp.sum %>% 
              dplyr::filter(., clst %in% keep.clst)
            k = 1
            
          }else{
            
            # loop-control
            k = k + 1
            
          }
          
          # Loop-control
          if(k > nrow(temp.sum)){break}
          
        }
        
      }
      
      ring.summary
      
    }else{
      
      stop("Undefined Method: Please choose either '1st' or 'All'")
      
    }
    
  }
  
  return(ring.summary)
  
}



## Segment Structure Spline smoothing
#'@description 
#' This is to automate the smoothing spline process,
#' This helps to smooth jagged boundaries & missing sections 
#' within grouped ring segments.
#'@source
segment.spline <- function(keep.segment, 
                           smooth.spline.df = 60, 
                           Padding = NULL,
                           name.write = "ref",
                           Morph.step = 0.0001){
  
  # Ensure Data Structure
  .step   = Morph.step*0.5
  .digits = log10(Morph.step) %>% abs()
  .pi <- pi %>% round(digits = .digits)
  keep.segment$theta <- keep.segment$theta %>% round(digits = .digits)
  
  # Padding For Smooth.Spline
  if(purrr::is_empty(Padding)){Padding <- 1.0000}
  if(Padding < min(keep.segment$theta)){
    message("Note: Spline with the Start-Missing more than 1 radius")
    while (Padding < min(keep.segment$theta)) {
      Padding <- Padding + 1.0000
    }
    print(paste("Padding at:", Padding))
  }
  ..keep.segment <- 
    keep.segment %>% 
    dplyr::filter(theta < Padding) %>% 
    dplyr::mutate(theta = theta + 2 * .pi)
  keep.segment <- rbind(keep.segment, ..keep.segment)
  
  # Check
  # plot(keep.segment$theta, keep.segment$dist)
  # range(keep.segment$theta)
  
  # Smooth.spline
  keep.smsp <- smooth.spline(keep.segment$theta,
                             keep.segment$dist,
                             df = smooth.spline.df)
  keep.sp   <- predict(keep.smsp, 
                       x = seq(min(keep.segment$theta),
                               max(keep.segment$theta), 
                               by = Morph.step)) %>% as.data.frame()
  keep.sp$x <- keep.sp$x %>% round(digits = .digits)
  names(keep.sp) <- c("theta", "dist")
  # keep.sp$clst   <- name.write
  
  # Get Spline Segments at Original Position
  .keep.L <- keep.sp %>% dplyr::filter(theta < (Padding - .step))
  .keep.R <- 
    keep.sp %>% 
    dplyr::filter(theta > (2 * .pi + .step)) %>% 
    dplyr::mutate(theta = theta - 2 * .pi)
  
  # Ensure Data Structure
  .keep.L$theta <- .keep.L$theta %>% round(digits = .digits)
  .keep.R$theta <- .keep.R$theta %>% round(digits = .digits)
  
  # Check
  # .keep.L$theta %>% range()
  # .keep.R$theta %>% range()
  
  # .keep.L & .keep.R
  .seq   <- seq(min(.keep.L$theta), max(.keep.R$theta), by = Morph.step)
  .weight.O <- seq(from = 1, to = 0, along.with = .seq)
  
  # .keep.L
  #'@note
  #' .keep.L splined after padding, 
  #' Hence, section between original edge (now as max(.keep.R$theta))
  #' Should have weight as 1.
  if(max(.keep.R$theta) != (Padding - Morph.step)){
    
    .seq.L    <- seq((max(.keep.R$theta) + Morph.step), 
                     (Padding - Morph.step), 
                     by = Morph.step)
    .weight.L <- append(rev(.weight.O), rep(1, length(.seq.L)))
    
  }else{
    
    .weight.L <- rev(.weight.O)
    
  }
  
  # .keep.R
  #'@note
  #' .keep.R has the section between 0.0001 to seq.
  if(min(.keep.L$theta) != Morph.step){
    
    .seq.R <- seq(Morph.step, 
                  (min(.keep.L$theta) - Morph.step), 
                  by = Morph.step)
    .weight.R <- append(rep(1, length(.seq.R)), .weight.O)
    
  }else{
    
    .weight.R <- .weight.O
    
  }
  
  .keep.L   %>% nrow()
  .weight.L %>% length()
  .keep.R   %>% nrow()
  .weight.R %>% length()
  
  .keep.L <-
    .keep.L %>% 
    dplyr::mutate(.keep = c("none"),
                  theta = theta %>% round(digits = .digits),
                  dist.L = dist * .weight.L)
  .keep.R <-
    .keep.R %>% 
    dplyr::mutate(.keep = c("none"),
                  theta = theta %>% round(digits = .digits),
                  dist.R = dist * .weight.R)
  ..keep <-
    .keep.L %>% 
    dplyr::full_join(.keep.R, by = "theta") 
  ..keep$dist.L[is.na(..keep$dist.L)] <- 0
  ..keep$dist.R[is.na(..keep$dist.R)] <- 0
  ..keep <- 
    ..keep %>% 
    dplyr::mutate(.keep = c("none"),
                  clst  = name.write,
                  theta = theta %>% round(digits = .digits),
                  dist  = dist.L + dist.R) %>% 
    dplyr::select(names(keep.segment))
  
  # plot(..keep$theta, ..keep$dist)
  # range(..keep$theta)
  
  keep.sp <- 
    keep.sp %>% 
    dplyr::filter(theta > (Padding - .step)) %>% 
    dplyr::filter(theta < (2 * .pi + .step)) %>% 
    dplyr::mutate(clst  = name.write) %>% 
    dplyr::select(names(keep.segment)) %>% 
    rbind(..keep) %>% 
    dplyr::arrange(theta)
  
  keep.sp$theta <- keep.sp$theta %>% round(digits = .digits)
  
  # plot(keep.sp$theta, keep.sp$dist)
  # keep.sp$theta %>% range()
  
  return(keep.sp)
  
}



## Get Gap-Size from a grouped tree ring segment 
#'@description 
#' This function is used only in the 
#' grouped tree ring segment determination section 
#' (Hold Ring Segments Control).
#' This indicate where segments to be filled.
#'@note
#' See also "segment.gapfill"
#'@source
segment.gapsize <- function(ring.keep,
                            gap.thres  = 0.5,
                            Morph.step = 0.0001){
  
  # This is to compute a "list" object
  # which will indicate where segments to be filled.
  
  require(dplyr)
  
  #'[(0) Prepare data structure]
  .digits = log10(Morph.step) %>% abs()
  .pi <- pi %>% round(digits = .digits)
  
  #'[(1) Sorting]
  .ring.keep <- 
    ring.keep %>% 
    dplyr::arrange(theta)
  .ring.pad <- 
    .ring.keep[1,] %>% 
    dplyr::mutate(theta = (theta + 2*.pi))
  .ring.keep <- 
    ring.keep %>% 
    rbind(.ring.pad) %>% 
    dplyr::mutate(theta = theta %>% round(digits = .digits)) %>% 
    dplyr::arrange(theta)
  .keep.sort <- .ring.keep$theta
  
  #'[(2) Calculate differences]
  .keep.diff <- .keep.sort %>% diff()
  
  #'[(3) Identify gaps (differences > gap.thres)]
  gap_indices <- which(.keep.diff > gap.thres)
  
  # Return a data frame with positions and sizes of gaps
  if (length(gap_indices) > 0) {
    
    gap.index <- data.frame(Start   = .keep.sort[gap_indices],
                            End     = .keep.sort[gap_indices + 1],
                            Size    = .keep.diff[gap_indices],
                            d.start = .ring.keep$dist[gap_indices],
                            d.end   = .ring.keep$dist[gap_indices + 1])
    
  }else{
    
    gap.index <- data.frame(Start   = numeric(0), 
                            End     = numeric(0), 
                            Size    = numeric(0),
                            d.start = numeric(0),
                            d.end   = numeric(0))
    
  }
  
  return(gap.index)
  
}


## Fill Gaps from the index (segment.gapsize) 
#'@description 
#' This function is used together with segment.gapsize,
#' which help to secure the grouped tree ring segment for better smoothing.
#' As if there are huge gaps, 
#' the smoothing results will not cover the pattern from reference ring. 
#' This is not good as the final product will be considered 
#' as the next reference ring.
#'@note
#' See also "segment.gapsize"
#'@source
segment.gapfill <- function(gap_index, 
                            ring.ref, 
                            coverage = 1,
                            fill.clst = "fill",
                            Morph.step = 0.0001){
  
  # This function fills the gap that is indexed by the reference ring.
  # The coverage controls how much to be extracted from the center.
  
  #' [Prepare data structure]
  .digits = log10(Morph.step) %>% abs()
  .pi <- pi %>% round(digits = .digits)
  ring.pad <- 
    ring.ref %>% 
    dplyr::mutate(theta = (theta + 2*.pi) %>% round(digits = .digits))
  ring.ref <- 
    ring.ref %>% 
    rbind(ring.pad) %>% 
    dplyr::mutate(theta = theta %>% round(digits = .digits)) %>% 
    dplyr::arrange(theta)
  
  k = 1
  GapFill.list <- list()
  while (k <= nrow(gap_index)) {
    
    #' [Create Morphed ring.ref]
    # section.Morph for k gap
    gap.morph <- Section.Morph(Edge.01  = gap_index$Start  [k],
                               Edge.02  = gap_index$End    [k],
                               dist.01  = gap_index$d.start[k],
                               dist.02  = gap_index$d.end  [k],
                               ref.ring = ring.ref)
    
    #' [Extract Morphed ring.ref]
    # Gap Size -> Gap Sequence Length
    length.gap <- (gap_index$Size[k]*(10^.digits) - 1) %>% as.integer()
    
    # Calculate the length of the central sequence
    length.center <- ceiling(length.gap * coverage)
    
    # Calculate start and end indices
    start_index <- 1 + ceiling((length.gap - length.center) / 2)
    end_index   <- start_index + length.center - 1
    
    # Extract and return the central sequence
    gap.fill      <- gap.morph[start_index:end_index,]
    gap.fill$clst <- fill.clst
    
    # Constrain the fill in 0-2*pi
    .pos <- which(gap.fill$theta > 2*.pi)
    gap.fill$theta[.pos] <- 
      (gap.fill$theta[.pos] - 2*.pi) %>% round(digits = .digits)
    
    # Save
    GapFill.list[[k]] <- gap.fill
    
    # loop-control
    k = k + 1
    
  }
  
  return(GapFill.list)
  
}



## Wrapper for the 3-above functions to make detection code cleaner.
#'@description 
#' This function is used to wrap 
#' (1) segment.spline, (2) segment.gapsize, & (3) segment.gapfill
#' together for easier usage and cleaner code env.
#'@source
segment.smooth <- function(keep.segment, 
                           ring.ref,
                           gap.thres  = 0.5,
                           coverage   = 1,
                           name.write = "ref",
                           Morph.step = 0.0001,
                           smooth.spline.df = 60){
  
  # ENV
  require(dplyr)
  require(purrr)
  
  # Check Gaps in the keep.segment
  gap_index  <- 
    keep.segment %>% 
    segment.gapsize(gap.thres  = gap.thres,
                    Morph.step = Morph.step)
  
  # Spline Ring.keep differently based on gap_index
  if(purrr::is_empty(gap_index$Size) != TRUE){
    
    # message
    message("Huge Gaps Encountered during Segment.Smooth...")
    
    # Fill the gap with 70% reference from ring.2
    keep.smooth <-
      segment.gapfill(gap_index  = gap_index, 
                      ring.ref   = ring.ref, 
                      coverage   = coverage,
                      Morph.step = Morph.step) %>% 
      dplyr::bind_rows() %>% 
      rbind(., keep.segment) %>% 
      segment.spline(name.write = name.write,
                     Morph.step = Morph.step,
                     smooth.spline.df = smooth.spline.df)  %>% 
      dplyr::arrange(theta)
    
  }else{
    
    keep.smooth <- 
      keep.segment %>% 
      segment.spline(name.write = name.write,
                     Morph.step = Morph.step,
                     smooth.spline.df = smooth.spline.df) %>% 
      dplyr::arrange(theta)
    
  }
  
  # Ensure Data Structure
  .digits = log10(Morph.step) %>% abs()
  keep.smooth$theta <- keep.smooth$theta %>% round(digits = .digits)
  
  # Output
  return(keep.smooth)
  
}



## Compute each tree ring segment length in Cartesian coordinates
#'@description 
#' This function is used together with segment.gapsize,
#' which help to secure the grouped tree ring segment for better smoothing.
#' As if there are huge gaps, 
#' the smoothing results will not cover the pattern from reference ring. 
#' This is not good as the final product will be considered 
#' as the next reference ring.
#'@note
#' See also "segment.gapsize"
#'@source
segment.length <- function(ring.segment){
  
  require(dplyr)
  
  if(c("dist.ref") %in% names(ring.segment)){
    
    ring.length.sum <-
      ring.segment %>% 
      dplyr::group_by(clst) %>% 
      dplyr::arrange(theta) %>% 
      dplyr::mutate(d.theta = theta - dplyr::lag(theta))%>% 
      dplyr::reframe(npts = n(),
                     #det  = all(range(theta) == c(0.0001, 6.2832)),
                     .p   = which(d.theta > 0.001),
                     Edge.01 = theta[.p],
                     Edge.02 = theta[(.p-1)])
    
    ring.padID <- ring.length.sum$clst
    ring.subID <- dplyr::setdiff(unique(ring.segment$clst), ring.padID)
    
    #'[No Need to Pad Ring Segments]
    temp.keep <- ring.segment %>% dplyr::filter(., clst %in% ring.subID)
    keep.sum <-
      temp.keep %>% 
      dplyr::group_by(clst) %>% 
      dplyr::arrange(theta) %>% 
      dplyr::mutate(x    = dist * cos(theta),
                    y    = dist * sin(theta),
                    dx   = x - dplyr::lag(x),
                    dy   = y - dplyr::lag(y),
                    x.sq = dx^2,
                    y.sq = dy^2,
                    segment.l = sqrt(x.sq + y.sq)) %>% 
      dplyr::summarise(delta     = (dist - dist.ref) %>% mean(),
                       segment.l = sum(segment.l, na.rm = TRUE))
    
    #'[Padded Ring Segments]
    .pi = round(pi, digits = 4)
    temp.pad <- ring.segment %>% dplyr::filter(., clst %in% ring.padID)
    temp.pad <- 
      temp.pad %>%
      dplyr::left_join(., ring.length.sum) %>% 
      dplyr::mutate(det    = theta < (Edge.01 - 0.0009),
                    .theta = theta + 2*.pi*det) %>% 
      # dplyr::filter(., det == TRUE)
      dplyr::mutate(.keep    = c("none"),
                    clst     = clst,
                    theta    = .theta,
                    dist     = dist,
                    dist.ref = dist.ref)
    pad.sum <-
      temp.pad %>% 
      dplyr::group_by(clst) %>% 
      dplyr::arrange(theta) %>% 
      dplyr::mutate(x    = dist * cos(theta),
                    y    = dist * sin(theta),
                    dx   = x - dplyr::lag(x),
                    dy   = y - dplyr::lag(y),
                    x.sq = dx^2,
                    y.sq = dy^2,
                    segment.l = sqrt(x.sq + y.sq)) %>% 
      dplyr::summarise(delta     = (dist - dist.ref) %>% mean(),
                       segment.l = sum(segment.l, na.rm = TRUE))
    
  }else{
    
    ring.length.sum <-
      ring.segment %>% 
      dplyr::group_by(clst) %>% 
      dplyr::arrange(theta) %>% 
      dplyr::mutate(d.theta = theta - dplyr::lag(theta))%>% 
      dplyr::reframe(npts = n(),
                     #det  = all(range(theta) == c(0.0001, 6.2832)),
                     .p   = which(d.theta > 0.001),
                     Edge.01 = theta[.p],
                     Edge.02 = theta[(.p-1)])
    
    ring.padID <- ring.length.sum$clst
    ring.subID <- dplyr::setdiff(unique(ring.segment$clst), ring.padID)
    
    #'[No Need to Pad Ring Segments]
    temp.keep <- ring.segment %>% dplyr::filter(., clst %in% ring.subID)
    keep.sum <-
      temp.keep %>% 
      dplyr::group_by(clst) %>% 
      dplyr::arrange(theta) %>% 
      dplyr::mutate(x    = dist * cos(theta),
                    y    = dist * sin(theta),
                    dx   = x - dplyr::lag(x),
                    dy   = y - dplyr::lag(y),
                    x.sq = dx^2,
                    y.sq = dy^2,
                    segment.l = sqrt(x.sq + y.sq)) %>% 
      dplyr::summarise(segment.l = sum(segment.l, na.rm = TRUE))
    
    #'[Padded Ring Segments]
    .pi = round(pi, digits = 4)
    temp.pad <- ring.segment %>% dplyr::filter(., clst %in% ring.padID)
    temp.pad <- 
      temp.pad %>%
      dplyr::left_join(., ring.length.sum) %>% 
      dplyr::mutate(det    = theta < (Edge.01 - 0.0009),
                    .theta = theta + 2*.pi*det) %>% 
      # dplyr::filter(., det == TRUE)
      dplyr::mutate(.keep = c("none"),
                    clst  = clst,
                    theta = .theta,
                    dist  = dist)
    pad.sum <-
      temp.pad %>% 
      dplyr::group_by(clst) %>% 
      dplyr::arrange(theta) %>% 
      dplyr::mutate(x    = dist * cos(theta),
                    y    = dist * sin(theta),
                    dx   = x - dplyr::lag(x),
                    dy   = y - dplyr::lag(y),
                    x.sq = dx^2,
                    y.sq = dy^2,
                    segment.l = sqrt(x.sq + y.sq)) %>% 
      dplyr::summarise(segment.l = sum(segment.l, na.rm = TRUE))
    
  }
  
  #'[Results]
  length.sum <- rbind(keep.sum, pad.sum)
  return(length.sum)
  
}



## Compute angular positions for tree ring segments
#'@description 
#' Consider the angular grouping range (column = Ang.pos),
#' classify angular positions within each tree ring segment.
#'@source
segment.AngularCompute <- function(df.mt, theta.range){
  
  require(dplyr)
  
  # Data
  df.clst <- df.mt
  
  theta.range
  n.theta = 360/theta.range
  Ang.pos <- replicate(length(dist), n.theta)
  i = 1
  while (TRUE) {
    
    .theta <- 2*pi - i*theta.range/180*pi
    temp.p <- which(df.clst$theta < .theta)
    Ang.pos[temp.p] <- (n.theta-i)
    
    # Loop
    i = i+1
    if(i == n.theta){break}
    
  }
  
  # Check
  # unique(Ang.pos) %>% table()
  
  df.clst$Ang.pos <- Ang.pos
  
  return(df.clst)
  
}



## Get trust clusters from tree ring segments
#'@description 
#' Take input from segment.AngularCompute,
#' checking which segment clusters have complete tree ring structure.
#'@source
segment.ClustTrust <- function(df.clst, theta.range){
  
  #'[Preparation:]
  n.theta = 360/theta.range
  
  clst.summary <- 
    df.clst %>% 
    group_by(clst, Ang.pos) %>%
    summarise(n.point = n(),
              mn.dist = mean(dist),
              .groups = "drop")
  
  clst.visibility <- 
    clst.summary %>%
    group_by(clst) %>%
    summarise(visibility = n(),
              dist = mean(mn.dist)) %>% 
    arrange(desc(dist))
  
  # Check
  # clst.visibility$visibility %>% range()
  
  # Visibility check to determine the trusted clusters
  trust.freq   <- n.theta
  if(max(clst.visibility$visibility) < trust.freq){
    trust.freq <- max(clst.visibility$visibility)
  }
  
  clst.p     <- which(clst.visibility$visibility >= trust.freq)
  clst.trust <- clst.visibility$clst[clst.p]
  
  return(clst.trust)
  
}



## Middle Search (Morphological Free Search) based on given segment clusters
#'@description 
#' Used *ONLY* within the *Middle Check* of the detection code section.
#' This offers the ability to morph and group neighboring segments
#' by assigning a targeted tree ring segment cluster.
#'@note
#' Free search without limited search range & too large grouped result
#' is REALLY *ERROR-Prone*. 
#' In case this function is used for other scenarios,
#' DO REMEMBER to CHECK the RESULTS GRAPHICALLY.
#'@source
segment.MiddleCheck <- function(target.segment.df,
                                ref.ring,
                                considered.segments,
                                npts.range,
                                segment.size,
                                segment.SUM,
                                Force.small = FALSE,
                                ReSample.Size = 0.0001){
  
  #'[Preparation:]
  #'[Environment]
  require(dplyr)
  
  #'[Standardized Parameters]
  ReSample.Size
  .digits <- log10(ReSample.Size) %>% abs()
  .M.Step <- ReSample.Size * 0.5
  .M.npts <- 10^.digits
  .pi     <- round(pi, digits = .digits)
  
  #'[target.segment.df]
  .mid.big <- 
    target.segment.df %>% 
    dplyr::mutate(theta = round(theta, digits = .digits))
  
  #'[target.segment.ID]
  morph.id <- target.segment.df$clst %>% unique()
  
  #'[Considered.segments within certain theta range ..ring.ck]
  ..ring.ck <- 
    considered.segments %>% 
    dplyr::mutate(theta = round(theta, digits = .digits))
  
  #'[Reference Tree Ring Structure]
  ring.2 <- 
    ref.ring %>% 
    dplyr::mutate(theta = round(theta, digits = .digits))
  
  #'[Reference Tree Ring Segment Summary]
  lookup.SUM <- segment.SUM
  
  #'[target.segment.sum]
  .mid.big.sum <-
    .mid.big  %>% 
    ring.join(., ring.segment = ref.ring) %>%   
    dplyr::mutate(., delta = abs(dist - dist.ref)) %>% 
    dplyr::group_by(clst) %>% 
    dplyr::arrange(delta) %>% 
    dplyr::summarise(npts      = n(),
                     m.delta   = mean(delta[1:round(npts.range*npts)]),
                     sd.delta  = sd  (delta[1:round(npts.range*npts)]),
                     # m.delta   = mean(delta[1:(npts*0.05)]),
                     # sd.delta  = sd  (delta[1:(npts*0.05)]),
                     # m.Ang.pos = npts / nrow(ring.2),
                     max.theta = max(theta),
                     min.theta = min(theta),
                     dist.l    = dist[which.min(theta)],
                     dist.r    = dist[which.max(theta)]) %>%
    as.data.frame()
  
  #'[Middle_Search_03:Consecutive While-loop Search]
  #'@description 
  #' Start the free search from suspicious Big-Segment.
  #' (1) morph.id: Targeted Tree ring segment ID*s*;
  #' (2) .mid.big: Targeted Tree ring segment structure data;
  #'@note
  #' *clst* Name in .mid.big & .mid.big.sum is *Unique*.
  #' It will only be the one that was targeted.
  #'@signal
  #' (1) mid.final
  mid.final <- FALSE
  while (TRUE) {
    
    #'[03_01:Update Ref Ring Structure (Ring.2) for Middle Search]
    #'@note
    #' (1) *mid.ring.2*: is the ref ring structure (**see padding).
    #'     This will be the free search backbone from targeted 
    #'     big-segment (1st morph.id / .mid.big.sum)
    #' (2) *.mid.big.sum*: indicates the targeted tree ring segment,
    #'     and its current neighbor grouping extension. Hence, 
    #'     - 1. Edge boundaries (min / max) & 
    #'     - 2. Edge dist range (dist.l, and dist.r)
    #'     are updated through the while-loop; & 
    #'     - 3. Targeted cluster ID remains.
    
    # Update Edge and corresponding dist
    mid.edge.01 <- .mid.big.sum$min.theta
    mid.edge.02 <- .mid.big.sum$max.theta
    mid.dist.01 <- .mid.big.sum$dist.l
    mid.dist.02 <- .mid.big.sum$dist.r
    
    # mid.ref as reference ring for Middle Search
    #'@note
    #' ring.2 / mid.ref: reference tree ring structure
    mid.ref <- ring.2
    .pos    <- which(mid.ref$theta < (mid.edge.02))
    mid.ref$theta[.pos] <- mid.ref$theta[.pos] + 2*.pi
    
    mid.ring.2  <- Section.Morph(Edge.01  = mid.edge.02,
                                 Edge.02  = mid.edge.01 + 2*.pi,
                                 dist.01  = mid.dist.02,
                                 dist.02  = mid.dist.01,
                                 ref.ring = mid.ref)
    
    .pos <- which(mid.ring.2$theta > max(ring.2$theta))
    mid.ring.2$theta[.pos] <- mid.ring.2$theta[.pos] - 2*.pi
    mid.ring.2 <- rbind(mid.ring.2, .mid.big[,names(mid.ring.2)])
    mid.ring.2$clst <- "ref"
    
    # Ensure Data Structure
    mid.ring.2$theta <- mid.ring.2$theta %>% round(digits = .digits)
    
    # Check
    # plot(mid.ring.2$theta,
    #      mid.ring.2$dist,
    #      cex = 0.3, col = "lightblue")
    # points(.mid.big$theta, .mid.big$dist, 
    #        cex = 0.5, col = "darkgreen")
    
    # points(mid.ring.2$theta,
    #        mid.ring.2$dist,
    #        cex = 0.3, col = "lightblue")
    
    #'[03_07:"mid.final" Adjustments]
    #'@special
    #' When *mid.final == TRUE*, "Further search" was triggered.
    #' As it changed ".mid.big.sum" table for further-search building
    #' new morphed-reference-ring structure "mid.ring.2" for checking
    #' possible neighboring segments, *after build of mid.ring.2*
    #' corrections on .mid.big.sum table should be made to
    #' justified the searched neighboring segments.
    #' For example, 
    #' the considered segments subset from ring.ck
    #' is based on mid.edge.01 / mid.edge.02. 
    #' And these Edge positions depends on ".mid.big.sum".
    if(mid.final == TRUE){
      
      # .mid.id  <- .mid.big.sum$clst
      # .mid.pos <- which(temp.section.sum$clst == .mid.id)
      
      # "..mid.big.sum" from previous check
      ..mid.big.sum
      
      # Resume ".mid.big.sum" from Further Check
      if(mid.edge.01 == ..mid.big.sum$max.theta){
        .mid.big.sum$min.theta <- ..mid.big.sum$min.theta
      }
      
      if(mid.edge.02 == ..mid.big.sum$min.theta){
        .mid.big.sum$max.theta <- ..mid.big.sum$max.theta
      }
      
      # Update mid.edge & mid.dist
      mid.edge.01 <- .mid.big.sum$min.theta
      mid.edge.02 <- .mid.big.sum$max.theta
      
      mid.dist.01 <- 
        mid.ring.2$dist[which(mid.ring.2$theta == mid.edge.01)]
      mid.dist.02 <- 
        mid.ring.2$dist[which(mid.ring.2$theta == mid.edge.02)]
      
    }
    
    #'[03_02:Middle Search Summary Table]
    #'@note
    #' Free search of Big-Segment within the constrained 
    #' segments (..ring.ck). 
    #' Align with its reference tree ring structure,
    #' the fitted neighboring segments will be picked up and 
    #' updates the while-loop search.
    
    # Alignment of tree ring segments and morphed structure
    temp.section.mid <-
      ..ring.ck %>% 
      segment.filter(Edge.01 = (mid.edge.01-2),
                     Edge.02 = (mid.edge.02+2),
                     ring.segments = .) %>% 
      ring.join(., ring.segment = mid.ring.2)
    
    head(temp.section.mid)
    
    if(purrr::is_empty(temp.section.mid$clst)){
      
      # Nothing in between
      temp.section.sum <- temp.ring.sc.sum %>% as.data.frame()
      temp.section.sum <- temp.section.sum[vector(),]
      
    }else if(Force.small == TRUE){
      
      # Something in between
      temp.section.sum <-
        temp.section.mid %>% 
        dplyr::mutate(., delta = abs(dist - dist.ref)) %>% 
        dplyr::group_by(clst) %>% 
        dplyr::arrange(delta) %>% 
        dplyr::summarise(npts      = n(),
                         m.delta   = mean(delta[1:round(npts.range*npts)]),
                         sd.delta  = sd  (delta[1:round(npts.range*npts)]),
                         max.theta = max(theta),
                         min.theta = min(theta),
                         dist.l    = dist[which.min(theta)],
                         dist.r    = dist[which.max(theta)]) %>% 
        dplyr::filter(., m.delta < 5) %>% 
        dplyr::left_join(., lookup.SUM, by = "clst") %>% 
        segment.clean(dist.thres = 3) %>% 
        segment.clean(dist.thres = 2, 
                      segment.size = segment.size) %>% 
        segment.overlap(segment.size  = segment.size) %>% 
        as.data.frame()
      
    }else{
      
      # Something in between
      temp.section.sum <-
        temp.section.mid %>% 
        dplyr::mutate(., delta = abs(dist - dist.ref)) %>% 
        dplyr::group_by(clst) %>% 
        dplyr::arrange(delta) %>% 
        dplyr::summarise(npts      = n(),
                         m.delta   = mean(delta[1:round(npts.range*npts)]),
                         sd.delta  = sd  (delta[1:round(npts.range*npts)]),
                         # m.delta   = mean(delta[1:(npts*0.05)]),
                         # sd.delta  = sd  (delta[1:(npts*0.05)]),
                         # m.Ang.pos = npts / nrow(ring.2),
                         max.theta = max(theta),
                         min.theta = min(theta),
                         dist.l    = dist[which.min(theta)],
                         dist.r    = dist[which.max(theta)]) %>% 
        dplyr::filter((npts<(0.2*.M.npts) & abs(dist.l-dist.r)>10) != TRUE) %>% 
        dplyr::filter(., m.delta < 5) %>% 
        dplyr::left_join(., lookup.SUM, by = "clst") %>% 
        # segment.clean(., dist.thres = 2) %>% 
        segment.clean(dist.thres = 3) %>% 
        segment.clean(dist.thres = 2, 
                      segment.size = segment.size) %>% 
        segment.overlap(segment.size  = segment.size,
                        overlap.check = "All") %>% 
        as.data.frame()
      
    }
    
    # Check
    temp.section.sum
    
    #'[03_03:Close-by search within Middle Search Section]
    #'@note
    #' use old code section,
    #' note that:
    #' (1) e1 is for search towards Edge.02
    #' (2) e2 is for search towards Edge.01
    
    # like section "Segment close-by"
    .ringi_closeby.e1 <- 
      temp.section.sum %>% 
      dplyr::filter(min.theta > .mid.big.sum$max.theta) %>% 
      segment.select(neighbor.thres = 0.1,
                     dist.check = TRUE,
                     dist.thres = 10)
    .ringi_closeby.e2 <- 
      temp.section.sum %>% 
      dplyr::filter(max.theta < .mid.big.sum$min.theta) %>% 
      segment.select(neighbor.thres = 0.1,
                     rev = TRUE,
                     dist.check = TRUE,
                     dist.thres = 10)
    
    #'Fix
    if(any(.ringi_closeby.e1$segment.l > segment.size)){
      
      .pos <- 
        which(.ringi_closeby.e1$segment.l > segment.size) %>% 
        tail(., n = 1)
      .ringi_closeby.e1 <- .ringi_closeby.e1[c(1:.pos),]
      
    }
    if(any(.ringi_closeby.e2$segment.l > segment.size)){
      
      .pos <- 
        which(.ringi_closeby.e2$segment.l > segment.size) %>% 
        tail(., n = 1)
      .ringi_closeby.e2 <- .ringi_closeby.e2[c(1:.pos),]
      
    }
    
    #'[03_04:Search Direction & Ending Control]
    #'@note 
    #' Argument "# Final search" ONLY for helping
    #' search again the other side of .mid.ID
    #' when morphing from different position sometimes 
    #' yields different results.
    #' However, when NOTHING is searched before, 
    #' Makes NO SENCE to go further.
    if(purrr::is_empty(.ringi_closeby.e1$clst) &
       purrr::is_empty(.ringi_closeby.e2$clst)){
      
      #'[03_04_01:No close-by Segments]
      #'@note
      #' This is a loop-control statement,
      #' when NO "close-by" segments, we need to assure that:
      #' (1) Arguments of The Further Search;
      #' (2) Direction of The Further Search;
      #' (3) Whether it is TRUELY Nothing Found.
      
      # Preparation
      .mid.id  <- .mid.big.sum$clst
      #'Fix: temp.section.sum might not always have .mid.id
      # .mid.pos <- which(temp.section.sum$clst == .mid.id)
      ..mid.big.sum <-
        temp.section.mid %>%  
        dplyr::filter(clst == .mid.id) %>% 
        dplyr::mutate(., delta = abs(dist - dist.ref)) %>% 
        dplyr::group_by(clst) %>% 
        dplyr::arrange(delta) %>% 
        dplyr::summarise(npts      = n(),
                         m.delta   = mean(delta[1:round(npts.range*npts)]),
                         sd.delta  = sd  (delta[1:round(npts.range*npts)]),
                         # m.delta   = mean(delta[1:(npts*0.05)]),
                         # sd.delta  = sd  (delta[1:(npts*0.05)]),
                         # m.Ang.pos = npts / nrow(ring.2),
                         max.theta = max(theta),
                         min.theta = min(theta),
                         dist.l    = dist[which.min(theta)],
                         dist.r    = dist[which.max(theta)]) %>%
        as.data.frame()
      
      #'[01:Decision...]
      #'@note
      #' when further search triggered, the ".mid.big" were modified.
      #' Once it was done before, no need to do another search.
      if(..mid.big.sum$min.theta == min(.mid.big$theta) &
         ..mid.big.sum$max.theta == max(.mid.big$theta)){
        
        message("... Nothing Found")
        break
        
      }else{
        
        message("# Final search ...")
        
      }
      
      #'[02:Further Search...or ENDs]
      if(..mid.big.sum$min.theta == min(.mid.big$theta)){
        
        # Signal
        mid.final <- TRUE
        
        # Report
        print("... Further Search (L)")
        
        # Safe ".mid.big"
        ..mid.big <- .mid.big
        
        # Modify ".mid.big"
        .mid.big <- 
          .mid.big %>% 
          dplyr::filter(theta > (..mid.big.sum$max.theta - .M.Step))
        
        # stop("Check From Here...")
        
      }else if(..mid.big.sum$max.theta == max(.mid.big$theta)){
        
        # Signal
        mid.final <- TRUE
        
        # Report
        print("... Further Search (R)")
        
        # Safe ".mid.big"
        ..mid.big <- .mid.big
        
        # Modify ".mid.big"
        .mid.big <- 
          .mid.big %>% 
          dplyr::filter(theta < (..mid.big.sum$min.theta + .M.Step))
        
      }else{
        
        # Resume ".mid.big"
        if(mid.final == TRUE){
          .mid.big  <- ..mid.big
          mid.final <- FALSE
        }
        
        # Loop-control
        break
        
      }
      
    }else{
      
      #'[03_04_02:Process close-by Segments]
      #'@note
      #' This section aims at post-processing close-by segment groups
      #' (both .ringi_closeby.e2 or .ringi_closeby.e1).
      #' As in middle search (while-loop), 
      #' segments were continuously grouped to initiate new search,
      #' the list of "morph.id" is elongated accordingly.
      #' Hence, there will be in either scenarios in this section:
      #' (1) the further search is not activated (*mid.final*),
      #' (2) close-by segments were further searched product.
      #' In either case above, 
      #' as new close-by segments were found (or one of them), 
      #' the "mid.final" should be turned off to allow another 
      #' "further search" to search neighboring segments.
      
      # Safe temp.section.sum
      .temp.section.sum <- temp.section.sum
      
      #'[Situation_01:.ringi_closeby.e2]
      #' Only when there are something close-by, then make sense
      while (purrr::is_empty(.ringi_closeby.e2$clst)!=TRUE) {
        
        .range <-
          max(.ringi_closeby.e2$max.theta) -
          min(.ringi_closeby.e2$min.theta)
        #' Condition
        #' (0) Distance Limitation
        cond.0 <- 
          (mid.edge.01 - .ringi_closeby.e2$max.theta[1]) > 0.3 |
          (mid.dist.01 - .ringi_closeby.e2$dist.r[1]) %>% abs() > 10 
        #' (1) If Segments far (>0.1) and small
        cond.1 <- 
          (mid.edge.01 - .ringi_closeby.e2$max.theta[1]) > 0.1 &
          sum(.ringi_closeby.e2$segment.l) < segment.size
        #' (2) If Segments are small, and qualified with ...
        cond.2 <-
          .range < 0.1 &
          nrow(.ringi_closeby.e2) < 3 &
          sum(.ringi_closeby.e2$segment.l) < 0.8*segment.size
        #' (3) Gap Jump
        cond.3 <-
          (mid.edge.01 - .ringi_closeby.e2$max.theta[1]) < 0.1 &
          (mid.dist.01 - .ringi_closeby.e2$dist.r[1]) %>% 
          abs() %>% round() > 5
        if(cond.0 | cond.1 | cond.2 | cond.3){
          
          .temp.clst <- dplyr::setdiff( temp.section.sum$clst,
                                        .ringi_closeby.e2$clst)
          temp.section.sum <-
            temp.section.sum %>% dplyr::filter(clst %in% .temp.clst)
          
          # .ringi_closeby.e2 <- .ringi_closeby.e2[vector(),]
          
          .ringi_closeby.e2 <-
            temp.section.sum %>%
            dplyr::filter(max.theta < .mid.big.sum$min.theta) %>%
            segment.select(neighbor.thres = 0.05,
                           rev = TRUE,
                           dist.check = TRUE,
                           dist.thres = 5)
          
        }else{
          
          #'[Method Addition]
          #'@note 
          #' Get Big-segments close-by together as well.
          
          # Resume temp.section.sum
          temp.section.sum <- .temp.section.sum
          
          # Check Any Big-Segments close-by (***)
          .e2.lb <- tail(.ringi_closeby.e2$min.theta, n = 1)
          .e2.ld <- tail(.ringi_closeby.e2$dist.l, n = 1)
          ..temp.big.sum <-
            temp.section.sum %>% 
            dplyr::filter(max.theta < .e2.lb) %>% 
            dplyr::filter(segment.l > segment.size) %>% 
            segment.select(neighbor.thres = 0.2,
                           rev = TRUE,
                           dist.check = TRUE,
                           dist.thres = 5)
          if(purrr::is_empty(..temp.big.sum$clst)!=TRUE){
            
            .dt = .e2.lb - ..temp.big.sum$max.theta[1]
            .dd = .e2.ld - ..temp.big.sum$dist.r[1]
            if(.dt > 0.3 | abs(.dd) > 10){
              ..temp.big.sum <- ..temp.big.sum[vector(),]
            }
            
          }
          
          # Fix Get the segments in-between
          if(purrr::is_empty(..temp.big.sum$clst)!=TRUE){
            
            ...temp.big.sum <-
              temp.section.sum %>% 
              dplyr::filter(max.theta < .e2.lb) %>%
              dplyr::filter(min.theta > ..temp.big.sum$max.theta[1])
            ..temp.big.sum <- rbind(..temp.big.sum, ...temp.big.sum)
            
          }
          
          # Combine the Results:
          .ringi_closeby.e2 <- 
            rbind(.ringi_closeby.e2, ..temp.big.sum)
          
          break
          
        }
        
      }
      
      #'[Situation_02:.ringi_closeby.e1]
      #' Only when there are something close-by, then make sense
      while (purrr::is_empty(.ringi_closeby.e1$clst)!=TRUE) {
        
        .range <-
          max(.ringi_closeby.e1$max.theta) -
          min(.ringi_closeby.e1$min.theta)
        #' Condition
        #' (0) Distance Limitation
        cond.0 <- 
          (.ringi_closeby.e1$min.theta[1] - mid.edge.02) > 0.3 |
          (.ringi_closeby.e1$dist.l[1] - mid.dist.02) %>% abs() > 10
        #' (1) If Segments far (>0.1) and small
        cond.1 <- 
          (.ringi_closeby.e1$min.theta[1] - mid.edge.02) > 0.1 &
          sum(.ringi_closeby.e1$segment.l) < segment.size
        #' (2) If Segments are small, and qualified with ...
        cond.2 <-
          .range < 0.1 &
          nrow(.ringi_closeby.e1) < 3 &
          sum(.ringi_closeby.e1$segment.l) < 0.8*segment.size
        #' (3) Gap Jump
        cond.3 <-
          (.ringi_closeby.e1$min.theta[1] - mid.edge.02) < 0.1 &
          (.ringi_closeby.e1$dist.l[1] - mid.dist.02) %>% 
          abs() %>% round() > 5
        if(cond.0 | cond.1 | cond.2 | cond.3){
          
          .temp.clst <- dplyr::setdiff(temp.section.sum$clst,
                                       .ringi_closeby.e1$clst)
          temp.section.sum <-
            temp.section.sum %>% dplyr::filter(clst %in% .temp.clst)
          
          # .ringi_closeby.e1 <- .ringi_closeby.e1[vector(),]
          
          .ringi_closeby.e1 <-
            temp.section.sum %>%
            dplyr::filter(min.theta > .mid.big.sum$max.theta) %>%
            segment.select(neighbor.thres = 0.05,
                           dist.check = TRUE,
                           dist.thres = 5)
          
        }else{
          
          #'[Method Addition]
          #'@description Get Big-segments close-by together as well.
          
          # Resume temp.section.sum
          temp.section.sum <- .temp.section.sum
          
          # Check Any Big-Segments close-by (***)
          .e1.rb <- tail(.ringi_closeby.e1$max.theta, n = 1)
          .e1.rd <- tail(.ringi_closeby.e1$dist.r, n = 1)
          ..temp.big.sum <-
            temp.section.sum %>% 
            dplyr::filter(min.theta > .e1.rb) %>% 
            dplyr::filter(segment.l > segment.size) %>% 
            segment.select(neighbor.thres = 0.2,
                           dist.check = TRUE,
                           dist.thres = 5)
          if(purrr::is_empty(..temp.big.sum$clst)!=TRUE){
            
            .dt = ..temp.big.sum$min.theta[1] - .e1.rb
            .dd = ..temp.big.sum$dist.r[1] - .e1.rd
            if(.dt > 0.3 | abs(.dd) > 10){
              ..temp.big.sum <- ..temp.big.sum[vector(),]
            }
            
          }
          
          # Fix Get the segments in-between
          if(purrr::is_empty(..temp.big.sum$clst)!=TRUE){
            
            ...temp.big.sum <-
              temp.section.sum %>% 
              dplyr::filter(min.theta > .e1.rb) %>%
              dplyr::filter(max.theta < ..temp.big.sum$min.theta[1])
            ..temp.big.sum <- rbind(..temp.big.sum, ...temp.big.sum)
            
          }
          # Combine the Results:
          .ringi_closeby.e1 <- 
            rbind(.ringi_closeby.e1, ..temp.big.sum)
          
          break
          
        }
        
      }
      
      # Resume ".mid.big"
      if(mid.final == TRUE){
        .mid.big  <- ..mid.big
        mid.final <- FALSE
      }
      
      #'[Situation_03_Control:]
      if(purrr::is_empty(.ringi_closeby.e1$clst) &
         purrr::is_empty(.ringi_closeby.e2$clst)){break}
      
    }
    
    #'[03_05:Morph Section Updated & Smoothed]
    if(mid.final == FALSE){
      
      #'@special
      #' It seems that it is better 
      #' to compute a search from one side then move to the other.
      if(sum(.ringi_closeby.e1$npts) > sum(.ringi_closeby.e2$npts)){
        .ringi_closeby <- .ringi_closeby.e1
      }else{
        .ringi_closeby <- .ringi_closeby.e2
      }
      .morph.id  <- append( morph.id, .ringi_closeby$clst)
      .mid.sp.df <- ..ring.ck %>% dplyr::filter(clst %in% .morph.id)
      
      # spline degree of freedom decision
      .temp.mid.range <- .mid.sp.df$theta %>% range() %>% dist()
      .temp.mid.dist  <- .mid.sp.df$dist  %>% range() %>% dist()
      .temp.mid.npts  <- .mid.sp.df %>% nrow()
      if(.temp.mid.range > 3 | .temp.mid.dist > 100){
        .temp.mid.df <- 40
      }else if(.temp.mid.range > 2.0 | .temp.mid.dist > 50){
        .temp.mid.df <- 20
      }else if(.temp.mid.range > 0.5 | .temp.mid.npts > 1000){
        .temp.mid.df <- 10
      }else{
        .temp.mid.df <- 5
      }
      
      .mid.sp <- smooth.spline(x  = .mid.sp.df$theta,
                               y  = .mid.sp.df$dist,
                               df = .temp.mid.df)
      temp.mid.sp <-
        predict(.mid.sp, 
                x = seq(.mid.sp.df$theta %>% min(),
                        .mid.sp.df$theta %>% max(),
                        by = ReSample.Size)) %>%
        as.data.frame()
      names(temp.mid.sp) <- c("theta", "dist")
      temp.mid.sp$theta <- temp.mid.sp$theta %>% round(digits = .digits)
      
      # Check
      # plot(.mid.sp.df$theta, .mid.sp.df$dist, cex = 0.3)
      # lines(.mid.sp, col = "darkgreen")
      
      # Morph & update the loop-control
      # Update .mid.big
      .mid.big      <- temp.mid.sp
      .mid.big$clst <- .mid.big.sum$clst
      morph.id      <- append(morph.id, .morph.id) %>% unique()
      
    }
    
    #'[03_06:Update .mid.big.sum by updated .mid.big]
    #'@note
    #' The updated .mid.big comes from:
    #' (1) Further Search Section (mid.final == TRUE);
    #' (2) Smoothed Spline Section (mid.final == FALSE)
    #' Update .mid.big.sum
    .mid.big.sum$max.theta <- max(.mid.big$theta)
    .mid.big.sum$min.theta <- min(.mid.big$theta)
    .mid.big.sum$dist.r    <- .mid.big$dist[which.max(.mid.big$theta)]
    .mid.big.sum$dist.l    <- .mid.big$dist[which.min(.mid.big$theta)]
    
    # if(mid.final == TRUE){stop("Check")}
    
  }
  
  # Output Results
  MidSearch.list <- list()
  MidSearch.list$temp.section.sum <- 
    temp.section.mid %>% 
    dplyr::filter(clst %in% morph.id) %>% 
    dplyr::mutate(., delta = abs(dist - dist.ref)) %>% 
    dplyr::group_by(clst) %>% 
    dplyr::summarise(npts      = n(),
                     m.delta   = mean(delta[1:round(npts.range*npts)]),
                     sd.delta  = sd  (delta[1:round(npts.range*npts)]),
                     max.theta = max(theta),
                     min.theta = min(theta),
                     dist.l    = dist[which.min(theta)],
                     dist.r    = dist[which.max(theta)]) %>% 
    dplyr::left_join(., lookup.SUM, by = "clst")
  MidSearch.list$mid.big.sum      <- .mid.big.sum
  MidSearch.list$morph.id         <- morph.id
  MidSearch.list$mid.big          <- .mid.big
  MidSearch.list$mid.ring.2       <- mid.ring.2
  
  return(MidSearch.list)
  
}



segment.DetHold <- function(target.ring,
                            target.morph,
                            hold.id,
                            hold.morph,
                            check.segments,
                            ReSample.Size,
                            npts.range = 0.5){
  
  #'[Preparation] 
  require(dplyr)
  
  # Digits for round-up based on ReSample.Size .digits
  .digits = log10(ReSample.Size) %>% abs()
  
  # Target Tree Ring Structure
  ring.keep   <- target.ring
  .hold.morph <- target.morph
  ring.ck     <- check.segments
  
  # Ensure Data Structure
  ring.ck    $theta <-  ring.ck$theta    %>% round(digits = .digits)
  .hold.morph$theta <- .hold.morph$theta %>% round(digits = .digits)
  
  #'[Det_Check: How many unique cluster in ring.keep]
  #' How to define "Number of ring.keep" (see: [Special](b))
  .keep.rowlimt <- 30
  
  #' Preparation
  #'(1) Size (npts) def
  #'(2) Segment cluster N.
  ..keep.sum <- 
    ring.keep %>% 
    dplyr::group_by(clst) %>% 
    dplyr::summarise(npts = n())
  .keep.size <- ..keep.sum$npts %>% quantile(0.8) %>% as.numeric()
  .keep.nrow <- ..keep.sum %>% nrow()
  
  k = 1
  hold.det      <- vector()
  hold.sum.det  <- vector()
  hold.sum.list <- list()
  while (k <= length(hold.id)) {
    
    # Report
    print(k)
    
    # Found ID
    .keep.id  <- ring.keep$clst %>% unique()
    
    # Difference in Clusters
    diff.clst <- dplyr::setdiff(.keep.id, hold.id[[k]])
    
    #'(1) [hold.det]
    .hold.det <- (length(diff.clst)) / (length(.keep.id))
    
    #'[Adjust .hold.det if necessary]
    if(.hold.det > 0.1){
      
      #'(1.1) Compute .ck.sum
      .ck.sum <- 
        ring.ck %>%
        dplyr::filter(clst %in% hold.id[[k]]) %>% 
        dplyr::mutate(theta = round(theta, digits = .digits)) %>% 
        ring.join(., .hold.morph) %>% 
        dplyr::mutate(delta = dist - dist.ref) %>% 
        dplyr::group_by(clst) %>% 
        dplyr::summarise(npts = n(),
                         det  = mean(delta[1:round(npts.range*npts)]),
                         sd   = sd  (delta[1:round(npts.range*npts)]),
                         theta.l = min(theta),
                         theta.r = max(theta)) 
      ..ck.sum <-
        .ck.sum %>% 
        dplyr::filter((abs(det) + sd) > 3)
      
      #'[Scenarios:]
      #'@note
      #' ..ck.sum indicates segments far from ring.keep 
      #'(1) ..ck.sum EMPTY: no conflicts occur, det ~ 0.0001;
      #'(2) ..ck.sum Exist: how many & how big are clusters;
      #'(3)               : Special| Most (>0.8) segments are far from Morph
      
      #'(1.2) Update hold.det & diff.clst
      if(purrr::is_empty(..ck.sum$clst)){
        
        message("... Better Searched Adjustment")
        print(.ck.sum)
        
        # Difference in Clusters
        diff.clst <- ..ck.sum$clst
        .hold.det <- 0.0001
        
      }else if(nrow(..ck.sum) > 0.8 * nrow(.ck.sum)){
        
        message("... Different Structure Confirmed")
        
        # Difference in Clusters (Remaine the same)
        diff.clst
        .hold.det
      
      }else{
        
        #' ..ck.sum Exist: how many & how big are clusters
        #'@note
        #' *Def. .det = Unique.keep / keep*
        #' Two parts determine .det:
        #'[keep]
        #' (1) Segment Size;
        #' (2) Number of segments in ring.keep.
        #' Hence, the situations can be denoted as:
        #' (a) Any segments in "Unique.keep" are *Large in Size*; +
        #'     Segments of ring.keep is *Small in Number* of segments
        #' (b) ALL segments in "Unique.keep" are *Small in Size*; +
        #'     Segments of ring.keep is *Small in Number* of segments
        #' (c) Any segments in "Unique.keep" are *Large in Size*; +
        #'     Segments of ring.keep is *Large in Number* of segments
        #' (d) Any segments in "Unique.keep" are *Small in Size*; +
        #'     Segments of ring.keep is *Large in Number* of segments
        #'[Def. Unique.keep]
        #' .ck.sum is considering how similar it might be between considered
        #' structures (ex. hold.ID[[k]]) and morphed structure by ring.keep.
        #' However, considering ..ck.sum by sd > 3 to capture unique segments
        #' on ring keep will miss some unique segment structures within 
        #' gaps from considered structures.
        #' Hence, instead of ..ck.sum sd > 3, *..ck.sum by sd <= 3* considered.
        #' The corresponding theta ranges will be used to capture ring.keep
        #' segments as *Unique.keep*.
        #'[Special]
        #' (a) Computation below *WILL Ignore* segment cross 0 & 2pi region.
        #' (b) define "Number of ring.keep": nrow(ring.keep) > 20;
        #'     Assume 50% positions have been recognized in ring.keep,
        #'     under Resample.Size = 0.0001, total npts = 31416;
        #'     if nrow = 10, average segment Size ~ 3141 in ring.keep, 
        #'     the dist.clst can ONLY be <= 1 (Too little);
        #'     Hence,
        #'     if nrow = 30, average segment Size ~ 1000 in ring.keep, 
        #'     the dist.clst can ONLY be <= 3 (little, but better)
        
        #' Convert ..ck.sum by sd <= 3
        ..ck.sum <- .ck.sum %>% dplyr::filter((abs(det) + sd) <= 3)
        
        # Remove small segment that across 0 & 2pi
        ck.pos <- which(..ck.sum$theta.l == ReSample.Size &
                        ..ck.sum$theta.r == round(2*pi, digits = .digits))
        if(purrr::is_empty(ck.pos) != TRUE){..ck.sum <- ..ck.sum[-ck.pos,]}
        
        #' Update diff.clst by ..ck.sum
        #' Capture similar segments between ring.keep & hold.id[[k]]
        .k = 1
        ..kept.ck <- list()
        while (.k <= nrow(..ck.sum)) {
          
          ..temp.kept.ck <-
            ring.keep %>% 
            dplyr::filter(.,theta > ..ck.sum[.k,]$theta.l &
                            theta < ..ck.sum[.k,]$theta.r)
          ..kept.ck[[.k]] <- ..temp.kept.ck
          .k = .k + 1
          
        }
        ..kept.ck
        
        #' Segment Structure Captured by ..ck.sum
        ..kept.ck <- ..kept.ck %>% dplyr::bind_rows()
        kept.clst <- ..kept.ck$clst %>% unique()
        
        #' Unique Ring.keep Segments (..keep.ck)
        ..keep.ck <-
          ring.keep %>%
          dplyr::filter(!clst %in% kept.clst)
        
        #' Unique Ring.keep Segments Summary
        ..keep.ck.sum <-
          ..keep.ck %>%
          dplyr::group_by(clst) %>%
          dplyr::summarise(npts = n())
        
        # Update diff.clst
        diff.clst <- ..keep.ck.sum$clst
        
        # How to define "Number of ring.keep" (see: [Special](b))
        .keep.rowlimt
        
        if(any(..keep.ck.sum$npts > .keep.size) &
           .keep.nrow < .keep.rowlimt){
          
          #'@note
          #' Adjustments no need; Naturally det > 0.1
          .hold.det <- (length(diff.clst)) / (length(.keep.id))
          
        }else if(all(..keep.ck.sum$npts < .keep.size) &
                 .keep.nrow < .keep.rowlimt){
          
          #'@note
          #' Adjustments needed; Unique.keep / ring.keep
          .hold.det <- nrow(..keep.ck) / nrow(ring.keep)
          
        }else if(any(..keep.ck.sum$npts > .keep.size) &
                 .keep.nrow > .keep.rowlimt){
          
          #'@note
          #' Adjustments needed; Unique.keep / ring.keep
          .hold.det <- nrow(..keep.ck) / nrow(ring.keep)
          
        }else if(all(..keep.ck.sum$npts < .keep.size) &
                 .keep.nrow > .keep.rowlimt){
          
          #'@note
          #' if the condition still comes here, 
          #' .det value should simply be updated by ..keep.ck.sum
          .hold.det <- (length(diff.clst)) / (length(.keep.id))
          
        }
        
      }
      
    }
    
    #'(1.3) hold.det
    hold.det  <- append(hold.det, .hold.det)
    
    #'(2) [hold.sum.list]
    #'(2.1) Ensure Data Structure
    hold.diff       <- ring.ck %>% filter(clst %in% diff.clst)
    hold.diff$theta <- hold.diff$theta %>% round(digits = .digits)
    .hold.ref       <- hold.morph[[k]]
    .hold.ref$theta <- .hold.ref$theta %>% round(digits = .digits)
    
    #'(2.2) Compute & Summarize
    hold.sum <- 
      hold.diff %>% 
      ring.join(., .hold.ref) %>% 
      dplyr::mutate(delta = dist - dist.ref) %>% 
      dplyr::group_by(clst) %>% 
      dplyr::summarise(npts = n(),
                       det  = mean(delta[1:round(npts.range*npts)]),
                       sd   =   sd(delta[1:round(npts.range*npts)])) 
    hold.sum.list[[k]] <- hold.sum
    
    #'(3) [k-loop control]
    k = k + 1
    
  }
  
  # Result Output
  result.list <- list()
  result.list$hold.det      <- hold.det
  result.list$hold.sum.list <- hold.sum.list
  
  return(result.list)
  
}



segment.GroupRefine <- function(segment.df){
  
  .check <- segment.df$theta %>% table()
  .pos   <- which(.check > 1)
  if(purrr::is_empty(.pos)!=TRUE){
    
    # Preparation
    diff_values <- diff(sort(unique(segment.df$theta)))
    .digits <- max(nchar(sub(".*\\.", "", diff_values)))
    
    # replicates in segment.df
    ck.num <- .check[.pos] %>% names()
    .ckpos <- which(segment.df$theta %in% ck.num)
    
    # Separate the overlapped segment sections
    .ring.check  <- segment.df[ .ckpos,]
    .ring.remain <- segment.df[-.ckpos,]
    
    ## Save the clusters with Number IDs
    ..keep.id <- 
      .ring.check$clst %>% 
      unique() %>% 
      grep("^[0-9]+$", ., value = TRUE)
    .ring.mkeep <- 
      .ring.check %>% 
      dplyr::filter(clst %in% ..keep.id) %>% 
      dplyr::group_by(theta) %>% 
      dplyr::summarise(dist = mean(dist)) %>% 
      dplyr::mutate(clst  = ..keep.id[1],
                    theta = round(theta, digits = .digits))
    
    ## Merge any other morph segments
    .ring.morph <-
      .ring.check %>% 
      dplyr::filter(!clst  %in% ..keep.id) %>% 
      dplyr::filter(!theta %in% as.character(.ring.mkeep$theta)) %>% 
      dplyr::group_by(theta) %>% 
      dplyr::summarise(dist = mean(dist)) %>% 
      dplyr::mutate(clst  = "clst.M",
                    theta = round(theta, digits = .digits))
    
    ## Merge .ring.mkeep & .segment.df
    segment.df <- 
      rbind(.ring.mkeep, .ring.morph) %>% 
      dplyr::select(names(.ring.remain)) %>% 
      rbind(., .ring.remain)
    
  }
  
  return(segment.df)
  
}



# Not Possible...
morph.check <- function(ring.morph){
  
  # Sort the ring.morph
  .hold.morph <- temp.morph %>% dplyr::arrange(theta)
  
  de1.morph <- .hold.morph$dist - dplyr::lag(.hold.morph$dist)
  
  j = 1
  j.npts = 1000 # 100 npts = 0.01 theta
  morph.step = 0.0001
  de.morph    <- vector()
  cos.de.morph <- vector()
  while (TRUE) {
    
    l.de <- mean(de1.morph[c((1 + (j-1)):(j.npts + (j-1)))], na.rm = T)
    r.de <- mean(de1.morph[c((j.npts + j)):(2*j.npts + j)],  na.rm = T)
    
    de.morph <- append(de.morph, (r.de - l.de))
    
    m.lde <- matrix(c(morph.step, l.de), nrow = 1, ncol = 2)
    m.rde <- matrix(c(morph.step, r.de), nrow = 1, ncol = 2)
    cos.de <- 
      rowSums(m.lde * m.rde)/
      (sqrt(rowSums(m.rde^2))*sqrt(rowSums(m.lde^2)))
    cos.de.morph <- append(cos.de.morph, cos.de)
    
    # j-loop control
    j = j + 1
    if((2*j.npts + j) > length(de1.morph)){break}
  }
  
  plot(de1.morph)
  plot(de.morph)
  plot(cos.de.morph)
  
  
}

