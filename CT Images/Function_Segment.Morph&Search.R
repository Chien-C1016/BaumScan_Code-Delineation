

segment.Morph.Search <- function(ring.segments,
                                 rev.control,
                                 theta.range,
                                 segment.size,
                                 ReSample.Size,
                                 search.range = 25,
                                 sapw.control = NULL,
                                 .plotcheck = FALSE){
  
  #'[Segment_Grouping:]
  # 1. Preparation ----
  # Targeted Stem Tree Ring Segments
  ring.ReSample.rev <- ring.segments
  
  # 1.1_Parameters: ----
  #'[Preparation_01]
  ReSample.Size <- ReSample.Size
  .digits = log10(ReSample.Size) %>% abs()
  .M.Step = ReSample.Size * 0.5
  .pi = round(pi, digits = .digits) %>% as.character() %>% as.numeric()
  
  #'[Preparation_02]
  #'@description 
  #' (1) segment.size: How large should a big segment be.
  #' (2) search.range: How much pixel range should the reference search begin.
  #'                   If there is not enough segments to be considered, 
  #'                   the threshold will be automatically enlarged.
  #' (3) ang.npts: How many points were expected for each degree range.
  segment.size <- segment.size
  search.range 
  ang.npts     <- ((1/180*.pi) %>% round(., digits = .digits))/ReSample.Size
  
  #'[Preparation_03]
  # Trusted Tree Ring Segments with Complete Structure
  clst.trust <- 
    ring.ReSample.rev %>% 
    segment.ClustTrust(theta.range = theta.range)
  if(rev.control == TRUE){clst.trust <- clst.trust %>% rev()}
  length.clst.trust <- length(clst.trust)
  
  #'[Preparation_04]
  # sapwood condition control
  #'@note
  #' This defines how many "true" segments should be included to
  #' construct the reference tree ring.
  if(purrr::is_empty(sapw.control)){
    sapw.control <- append(0.3, rep(0.6, (length(clst.trust)-1)))
  }
  
  #'[Segment Grouping within whole stem structure]
  #'@description 
  #' ii-loop:
  #' ii controls the input of the Trust Tree ring Cluster.
  #'@Target
  #' ring.ReSample, ring.ReSample.rev
  #' Work on ring.ReSample.rev, ring.ReSample remain as a safe.
  #' (1) ring.1: only for initializing the ii-loop (other usages were removed)
  #' (2) ring.2: the reference ring, extracted from clst.trust.
  
  # -------------------------------------------------------------------------- #
  # ___II_START -------------------------------------------------------------
  # -------------------------------------------------------------------------- #
  
  #'[ii-loop Controls loop segment grouping within clst.trust]
  ii <- 1
  ring.1 <- ring.ReSample.rev %>% dplyr::filter(., clst %in% clst.trust[ii])
  ring.k.list <- list(ring.1[,c("clst", "theta", "dist")]) # Segments Keep
  ring.m.list <- list(ring.1[,c("clst", "theta", "dist")]) # Morphed structure
  ring.b.list <- list(ring.1[,c("clst", "theta", "dist")]) # References Base
  while (TRUE){
    
    # 1.2_Slice Segments ----
    
    #'control: Length of the trust clst
    length.clst.trust
    
    #'[Slice Ring Segments & Select Reference Tree Ring]
    #'@note
    #' This while-loop slice the segment data frame (ring.ReSample.rev)
    #' by the trust tree ring objects (clst.trust).
    #' It has the ability to move ii forward based on the general segment 
    #' theta coverage within the searching range (within clst.trust ii - ii+1)
    while (ii <= length.clst.trust) {
      
      message(ii)
      
      # ENSURE THETA STRUCTURE
      ring.ReSample.rev$theta <- 
        ring.ReSample.rev$theta %>% round(., digits = .digits)
      
      # Reset
      ring.ReSampled <- ring.ReSample.rev
      
      # Select Target Region
      if(length.clst.trust != 1 & length.clst.trust != ii){
        
        ring.ReSampled <- segment.slice(ring.ReSampled,
                                        cluster   = clst.trust[ii],
                                        slice.dir = "outer")
        ring.ReSampled <- segment.slice(ring.ReSampled,
                                        cluster  = clst.trust[ii+1],
                                        slice.dir = "inner")
        
        # Ring.2 as Morphing Base
        ring.2 <- 
          ring.ReSample.rev %>% 
          dplyr::filter(., clst == clst.trust[ii]) %>% 
          dplyr::select(clst, theta, dist) %>% 
          dplyr::mutate(clst = "ref")
        
      }else{
        
        ring.ReSampled <- segment.slice(ring.ReSampled,
                                        cluster  = clst.trust[ii],
                                        slice.dir = "outer")
        
        # Ring.2 as Morphing Base
        ring.2 <- 
          ring.ReSample.rev %>% 
          dplyr::filter(., clst == clst.trust[ii]) %>% 
          dplyr::select(clst, theta, dist) %>% 
          dplyr::mutate(clst = "ref")
        
      }
      
      # The control to next "Trust Cluster Segments"
      #'@note 
      #' if there are not enough segments, move the search to 
      #' the next "Trust Tree Ring Segment Cluster (ii)"
      #'@Fix
      #' Ensure at least a pair of numbers for range(), dist()
      dist.resample <-
        if(nrow(ring.ReSampled) > 1){
          range(ring.ReSampled$theta) %>% diff()}else{0}
      
      if(dist.resample < 2*.pi*0.7){
        
        # Output Complete-Ring Result
        # ring.k.list <- append(ring.k.list, list(ring.2))
        # ring.m.list <- append(ring.m.list, list(ring.2))
        ring.b.list <- append(ring.b.list, list(ring.2))
        
        ii = ii + 1
        if(ii < length.clst.trust){
          
          message("No Ring Segments in between Ring.ID: ", 
                clst.trust[ii],
                " and ",
                clst.trust[ii+1],
                ",\n", "Move to Next Searching Range")
          
        }else if(ii > length.clst.trust){
          
          message("No Ring Segments Left for Searching...")
          
        }else{
          
          message("Move to Final Trust Tree Ring for Searching...")
          
        }
        
      }else{break}
      
    }
    
    #'ii-loop control
    if(ii > length.clst.trust){break}
    
    # PlotCheck
    # For someone wants to check the sliced tree ring segments:
    # plot(ring.ReSampled$theta, ring.ReSampled$dist, cex = 0.3)
    
    # 1.3_Loop-Parameters (I|J|L) ----
    #'[Lookup Table of Sliced Segments for grouping consideration]
    #'@note
    #' dist.group parameters can be updated based on different scenarios.
    Check.ID <- ring.ReSampled$clst %>% unique()
    ring.ck <- 
      ring.ReSampled %>% 
      dplyr::filter(., clst %in% Check.ID) %>% 
      dplyr::select(., c("clst", "theta", "dist")) %>% 
      dist.group(., 
                 tol.degree   = 1.0,
                 dist.thres   = 1.0,
                 shrink_by    = 0.01,
                 Morph.step   = ReSample.Size,
                 Morph.method = "spline") %>%
      ring.join(., ring.segment = ring.2)
    
    lookup.SUM <-
      ring.ck %>%
      segment.length() %>% 
      dplyr::arrange(delta)
    
    ring.SUM <-
      lookup.SUM %>% 
      dplyr::filter(., segment.l > segment.size) %>%
      dplyr::arrange(delta) %>% 
      dplyr::filter(., delta < search.range) 
    
    # Fix
    if(purrr::is_empty(ring.SUM$clst)){
      ring.SUM <- lookup.SUM %>% slice(1)
    }
    
    lookup.SUM$delta <- NULL
    # big.segment.ID   <- ring.SUM$clst
    
    #'[Parameters01]
    #'@description
    #' Loop Controls:
    #' (1) i: which segment within ring.SUM for grouping its structure
    #' (2) l: Holding; for any possible tree ring structure from ring.SUM
    #' (3) j: output of each final grouped tree ring structure within ii-loop
    #' Holding Control:
    #' (1) hold.id
    #' (2) hold.keep
    #' (3) hold.morph
    #' (4) hold.value
    #' (5) hold.dist
    #' (6) hold.number
    #' (7) hold.clst.ck
    #' (8) Full.segment: Signal; indicate 1st tree ring candidate; 
    #'     Holding starts.
    #' Segment Grouping Result bins within each ii-loop:
    #' (1) group.ID
    #' (2) ring.keep.list
    #' (3) ring.morph.list
    #' (4) ring.base.list
    #' Determination Signal (End of Holding)
    #' (1) .update
    
    # i-loop
    i  = 1
    ni = nrow(ring.SUM)
    
    # j-loop
    j = 1
    group.ID        <- vector()
    ring.keep.list  <- list()
    ring.morph.list <- list()
    ring.base.list  <- list()
    .update         <- FALSE
    .update.2       <- FALSE
    
    # l-loop
    l = 1
    hold.id         <- list()
    hold.keep       <- list()
    hold.morph      <- list()
    hold.value      <- vector()
    hold.dist       <- vector()
    hold.number     <- 3
    hold.clst.ck    <- vector() #ring.SUM$clst[1:hold.number]
    Full.segment    <- FALSE
    
    #'[Parameters02]
    #'@description 
    #' 
    #' (1) Big.Group
    #' (2) S.Flip
    #' (3) Tail.check
    #' (4) Small.Back
    #' (5) Ring.i.tempclst
    #' (6) .plotcheck
    Big.Group  <- FALSE
    S.Flip     <- FALSE
    Tail.check <- FALSE
    Small.Back <- FALSE
    .plotcheck <- .plotcheck
    Ring.i.tempclst <- as.integer(99999)
    

    # ------------------------------------------------------------------------ #
    # ___I_START ------------------------------------------------------------- 
    # ------------------------------------------------------------------------ #
    
    while (TRUE) {
      
      # 2. Find and Morph Segments ----
      
      # For i-Check
      # if(i == 9){stop("I Check")}
      #'@note
      #' Use saved & un-padded ring.2 and ring.ck to re-run 
      # ring.2  <- .ring.2
      # ring.ck <- .ring.ck
      
      #' [Preparation]
      # Updated ring.2 (reference ring)
      ring.2
      .ring.2 <- ring.2
      head(ring.2)
      
      # Restore ring.ck
      ring.ck  <- ring.ck[,c("clst", "theta", "dist")]
      .ring.ck <- ring.ck
      head(ring.ck)
      
      # Restore lookup.sum
      # cluster = 99999 for representing another side of Ring.ID.i
      lookup.SUM <- lookup.SUM %>% dplyr::filter(clst != Ring.i.tempclst)
      
      # Select the Ring Segment I
      Ring.ID.i = ring.SUM$clst[i] %>% as.integer()
      ring.i <- ring.ck %>% dplyr::filter(., clst %in% Ring.ID.i)
      
      #' [Check Padding]
      #'@description 
      #' Whether the tree ring segment 
      #' passing through 0 and 2pi in polar coordinates.
      #'@note
      #' (1) ReSample.Size = 0.0001
      #' (2) .pi2 = 6.2832
      .pi2 = round(2*.pi, digits = .digits)
      if(all(range(ring.i$theta) == c(ReSample.Size, .pi2)) == TRUE){
        
        temp.k     <- ring.i[order(ring.i$theta),]
        temp.theta <- temp.k$theta
        temp.diff  <- temp.theta %>% diff() %>% abs
        .p <- which(temp.diff > 10*ReSample.Size)
        Edge.01 <- temp.theta[.p]
        Edge.02 <- temp.theta[(.p+1)]
        Padding <- FALSE
        dist.01 <- temp.k$dist[.p]
        dist.02 <- temp.k$dist[(.p+1)]
        
      }else{
        
        Edge.01 <- max(ring.i$theta)
        Edge.02 <- min(ring.i$theta)
        Padding <- TRUE 
        dist.01 <- ring.i$dist[which.max(ring.i$theta)]
        dist.02 <- ring.i$dist[which.min(ring.i$theta)]
        
      }
      
      Padding
      
      if(Padding == TRUE){
        
        # Padding ring.ck
        ring.ck <- segment.padding(ring.segments = ring.ck, 
                                   Padding.from = Edge.02)
        
        # Padding ring.2
        ring.2  <- segment.padding(ring.segments = ring.2, 
                                   Padding.from = Edge.02)
        
        # Show
        ring.i  <- ring.ck %>% dplyr::filter(., clst %in% Ring.ID.i)
        
        # Padding Edge.02, dist.02 remains
        Edge.02 <- Edge.02 + 2 * .pi
        
      }
      
      # Fix: Ensure Data Structure
      Edge.01 <- Edge.01 %>% round(digits = .digits)
      Edge.02 <- Edge.02 %>% round(digits = .digits)
      
      #' [Sectional Morph Preparation]
      # Morphing in between the section
      ring.m <- Section.Morph(Edge.01  = Edge.01,
                              Edge.02  = Edge.02,
                              dist.01  = dist.01,
                              dist.02  = dist.02,
                              ref.ring = ring.2)
      
      # Filter the ring segments in between 2 edges
      temp.1st <-
        ring.ck %>% 
        segment.filter(Edge.01 = Edge.01,
                       Edge.02 = Edge.02,
                       ring.segments = .) %>% 
        ring.join(., ring.segment = ring.m) %>% 
        dplyr::mutate(., delta = abs(dist - dist.Morph))
      
      head(temp.1st)
      
      # Summaries the filtered ring segments
      # 873 = ((5/360*2*.pi) %>% round(., digits = 4))/0.0001 # Step
      temp.1st.sum <-
        temp.1st %>% 
        dplyr::group_by(clst) %>% 
        dplyr::arrange(delta) %>% 
        dplyr::summarise(npts      = n(),
                         m.delta   = mean(delta[1:round(npts*0.1)]),
                         sd.delta  = sd  (delta[1:round(npts*0.1)]),
                         max.theta = max(theta),
                         min.theta = min(theta),
                         dist.l    = dist[which.min(theta)],
                         dist.r    = dist[which.max(theta)]) %>% 
        dplyr::filter(., m.delta < 5) %>% 
        dplyr::left_join(., lookup.SUM, by = "clst") %>% 
        # segment.clean(., dist.thres = 2) %>% 
        segment.clean(., dist.thres = 3) %>% 
        segment.clean(., dist.thres = 2, segment.size = segment.size) %>% 
        segment.overlap(., segment.size = segment.size) %>% 
        as.data.frame()
      temp.1st.sum
      
      #'[Plot Check]
      if(.plotcheck == TRUE){
        
        Temp.checkID <-
          temp.1st %>% 
          dplyr::group_by(clst) %>% 
          dplyr::arrange(delta) %>% 
          dplyr::summarise(npts      = n(),
                           m.delta   = mean(delta[1:round(npts*0.1)]),
                           sd.delta  = sd  (delta[1:round(npts*0.1)]),
                           max.theta = max(theta),
                           min.theta = min(theta),
                           dist.l    = dist[which.min(theta)],
                           dist.r    = dist[which.max(theta)]) %>% 
          dplyr::filter(., m.delta < 20) %>% 
          dplyr::left_join(., lookup.SUM, by = "clst") %>% 
          dplyr::select(clst) %>% 
          unlist() %>% 
          unname()
        Temp.ck.df <-
          temp.1st %>%
          dplyr::filter(., clst %in% Temp.checkID)
        
        test2 <- rbind(ring.i    [,c("clst", "theta", "dist")],
                       Temp.ck.df[,c("clst", "theta", "dist")])
        plot(test2$theta,
             test2$dist,
             col = as.factor(test2$clst),
             cex = 0.3)
        
      }
      
      #'[Find_n_Morph:]
      #'@description 
      #' k only controls the sequence of saving k-loop results within Ring.i,
      #' Loop is controlled by how many "big.segment" are in such search region
      #' Only consider the 1st position of the summary
      temp.border  <- c(Edge.01, Edge.02)
      .border.save <- temp.border
      
      # Morph & Search (i-loop) Targets:
      temp.ring.sc     <- temp.1st
      temp.ring.sc.sum <- temp.1st.sum
      big.segment.pos  <- which(temp.ring.sc.sum$segment.l > segment.size)
      
      
      if(purrr::is_empty(temp.ring.sc.sum$clst) != TRUE){
        
        # -------------------------------------------------------------------- #
        # ___K_START --------------------------------------------------------- 
        # -------------------------------------------------------------------- #
        
        k = 1
        small.segment.search = 2 * pi / 10
        neighbor.thres = small.segment.search / 10
        ring.keep  <- list()
        morph.list <- list()
        while (TRUE) {
          
          # k-Check
          # if(k == 2){stop("K Check")}
          # if(is_empty(which(is.na(ring.2$dist)))!=TRUE){stop("ring.2 Check")}
          
          # Update Check
          ring.m
          temp.ring.sc
          temp.ring.sc.sum
          
          # 2.1_Find Segments ----
          #'[00_Group_Close_By_Segments:]
          #'[Preparation]
          #'@note 
          #' to make sure small segments are trust worthy or not
          #' A certain theta range should not be decided by 
          #' extremely small size segments which might be just noises.
          #' The "temp.small" only make sense 
          #' when general search of big.clst is not empty./ 
          #' or it is really close to the starting/ ending edges.
          
          # Big-Segments
          big.segment.pos <- which(temp.ring.sc.sum$segment.l > segment.size)
          
          ## (01)_Tail Check ----
          #'[00_01_Tail Check]
          #'@description 
          #' This section aims for big segments that close to 
          #' either edges (Edge01 & Edge02), especially for those with
          #' (1) relative small theta distance,
          #' (2) mid-range to huge radius distance difference
          #' This situation looks like a huge gap-jump, 
          #' causing non-continuous Tree ring structure. 
          #' Hence,
          #' we use "Middle Check" to search independently 
          #' from such big segment and the corresponding Edge
          #' to make sure the result's continuity.
          #'@seealso 
          #' "_03_Middle Check"
          #'@condition
          #' Only when there are Big Segments & Close to either Edges
          if(purrr::is_empty(big.segment.pos) != TRUE){
            
            if(max(temp.ring.sc.sum$max.theta[big.segment.pos]) > 
               (Edge.02 - 1.5)){
              
              # Coding Check 
              # stop("back ring.2 adjustment...")
              
              big.segment.pos
              .temp.big.sum <- temp.ring.sc.sum[big.segment.pos,]
              .last.sum  <- .temp.big.sum[which.max(.temp.big.sum$max.theta),]
              
              temp.last <- 
                temp.ring.sc %>% dplyr::filter(clst == .last.sum$clst)
              # last.dist  <- temp.last$dist[which.max(temp.last$theta)]
              # last.delta <- abs(dist.02 - last.dist)
              
              # -------------------------------------------------------------- #
              # Huge-Slope Adjustment
              .temp.last <- 
                temp.last %>% 
                dplyr::arrange(desc(theta)) %>%
                dplyr::select(c("clst", "theta", "dist")) %>% 
                ring.join(., ring.m) %>% 
                dplyr::mutate(delta = abs(dist.Morph - dist))
              
              .pos <- c(1:round(.last.sum$npts*0.1))
              .morph.delta <- .temp.last$delta[.pos] %>% mean()
              .morph.sd    <- .temp.last$delta[.pos] %>%   sd()
              
              # Check
              # plot(.temp.last$theta[.pos], .temp.last$dist[.pos])
              # points(ring.m$theta, ring.m$dist, cex = 0.3, col = "lightblue")
              # -------------------------------------------------------------- #
              
              if(round((Edge.02 - .last.sum$max.theta), digits = 1) < 0.1 &
                 ((.morph.delta + .morph.sd) > 5)){
                
                # stop("Tail Check...")
                
                # Report
                message("# Tail Check ...")
                
                # Signal
                Tail.check <- TRUE
                
                print(.last.sum)
                
                #'[Add-in S-Flip]
                .gap.l <- .last.sum$min.theta - Edge.01
                .gap.r <- Edge.02 - .last.sum$max.theta
                S.Flip <- FALSE
                if(.gap.r < .gap.l){
                  
                  message("# Search-Flip from Tail #")
                  S.Flip <- TRUE
                  
                  # Check
                  print("# Section TAIL Special")
                  
                  big.segment.pos <- tail(big.segment.pos, n = 1)
                  
                }
                
              }
              
            }
          }
          
          
          
          ## (02)_Segment close-by ----
          #'[00_02 Ring.i "eats" the segments "close-by"]
          #'@description 
          #' Check from which side of Edges to make "close-by" search.
          #' Get targeted segment grouping if ANY.
          #'@export
          #' ringi_closeby
          
          #'[Preparation]
          #'@note 
          #' (1) close-by within the search range.(= 0.1)
          #' (2) The Segments should also "close-by" to each other.(= 0.05)
          #' (3) Total eaten-segment length should > 'segment.size' or 'npts'
          #'@import 
          #' (1) .ringi_closeby.e1: Segments Summary from Edge01
          #' (2) .ringi_closeby.e2: Segments Summary from Edge02 (Close to 2pi)
          .ringi_closeby.e1 <- segment.select(ring.summary = temp.ring.sc.sum,
                                              neighbor.thres = 0.1,
                                              dist.check = TRUE,
                                              dist.thres = 10)
          .ringi_closeby.e2 <- segment.select(ring.summary = temp.ring.sc.sum,
                                              neighbor.thres = 0.1,
                                              rev = TRUE,
                                              dist.check = TRUE,
                                              dist.thres = 10)
          
          #' While-loop close-by segments search:
          #'@condition
          #' (1) Only when there are something close-by, then make sense
          #' (2) If Tail.check Signals, 
          #'     it should be handled in section 03_Middle Check
          while (Tail.check != TRUE){
            
            # Control
            if(purrr::is_empty(.ringi_closeby.e1$clst) &
               purrr::is_empty(.ringi_closeby.e2$clst)){
              
              # Report
              message("# No close-by segments found")
              break
              
            }
            
            # Preparation
            big.segment.pos <- which(temp.ring.sc.sum$segment.l > segment.size)
            if(purrr::is_empty(big.segment.pos)!=TRUE){
              .gap.l <- 
                min(temp.ring.sc.sum$min.theta[big.segment.pos]) - Edge.01
            }else{
              .gap.l <- Inf
            }
            
            #'[Close-by Tree ring segments:]
            #'[01:Close-by Segments Search at Edge01] 
            #'@description 
            #' Check from which side to make "close-by" search.
            #' Close-by morph grouping 1st from Edge01, then Edge02.
            #'@signal
            #' (1) s.control
            #' (2) rev.control == FALSE
            rev.control    <- FALSE
            #'@condition
            #' # No need to continue when:
            #' (1) segments not close enough (by "segment.select"), or
            #' (2) close-search from Edge.01 is too far
            if(purrr::is_empty(.ringi_closeby.e1$clst) != TRUE){
              
              if((.ringi_closeby.e1$min.theta[1] - Edge.01) < 0.1){
                
                .ringi_closeby <- .ringi_closeby.e1
                s.control      <- TRUE
                
                
              }else{
                
                s.control      <- FALSE
                .ringi_closeby.e1 <- .ringi_closeby.e1[vector(),]
                
              }
              
            }else{
              
              s.control      <- FALSE
              
            }
            
            
            
            #'[02:Any-Big from Edge.01 & Edge.02]
            #'@special
            #' Additional Check for segment groups from Edge.01 & Edge.02
            #' If the any-Big segments exist,
            #' The search should follow such "Big-Segment".
            #'@signal
            #' (1) any.big.e1
            #' (2) any.big.e2
            #' (3) s.control == FALSE, when any.big.e* == TRUE
            
            # Any-Big from Edge.01
            if(any(.ringi_closeby.e1$segment.l>segment.size) |
               .gap.l < 0.1){
              
              # No need to continue when:
              # (2) The big-segment exist and close to Edge.01
              any.big.e1 <- TRUE
              s.control  <- FALSE
              
            }else{
              
              # Signal
              # (close.segments found from e1 are small segments)
              # (if s.control  == TRUE)
              any.big.e1 <- FALSE
              
            }
            
            # Any-Big from Edge.02
            if(any(.ringi_closeby.e2$segment.l>segment.size)){
              
              # Signal
              # (Ask for handling from edge.02)
              any.big.e2 <- TRUE
              s.control  <- FALSE
              
            }else{
              
              any.big.e2 <- FALSE
              
            }
            
            #'[03:small.edge.01 | small.edge.02]
            #'@description 
            #' If there is no small.edge.01 found, (s.control == FALSE)
            #' The search from tail should be activated.
            #'@signal
            #' (1) s.control;
            #' (2) rev.control: ONLY controls statements of small segments;
            #' (3) S.Flip: General Function across whole Tree Ring Grouping.
            #'@condition
            #' # No need to continue when:
            #' (1) segments not close enough (by "segment.select"), or
            #' (2) close-search from Edge.02 is too far; or
            #' (3) ** There are big segments from either Edge.
            if(purrr::is_empty(.ringi_closeby.e2$clst) != TRUE){
              
              if((Edge.02 - .ringi_closeby.e2$max.theta[1]) < 0.1 &
                 any.big.e1 == FALSE &    # No Big-Segments close to Edge.01
                 any.big.e2 == FALSE){    # No Big-Segments close to Edge.02
                
                # Compare the current "close-by" from e1 & e2
                # choose the larger one
                if(sum(.ringi_closeby.e2$npts) > sum(.ringi_closeby.e1$npts)){
                  
                  .ringi_closeby <- .ringi_closeby.e2
                  s.control      <- TRUE
                  rev.control    <- TRUE
                  S.Flip         <- TRUE
                  message("# Search-Flip from Tail #")
                  
                  # Check
                  print("# Section I")
                  
                }
                
              }else{
                
                # segments not close enough
                .ringi_closeby.e2 <- .ringi_closeby.e2[vector(),]
                
              }
              
            }
            
            #'[04:Fix Existence of ".ring_closeby"]
            if(exists(".ringi_closeby") != TRUE){
              .ringi_closeby <- .ringi_closeby.e1[vector(),]
            }
            
            #'[Neighboring Morphing Search Based on ringi_closeby]
            #'@description 
            #' 
            #'@export
            #' temp.small.list >>> ringi_closeby (list of grouped segments)
            #'@condition
            #' Only when segments are close enough, then make sense
            #' If s.control == FALSE, ringi_closeby will be EMPTY.
            if(s.control == TRUE){
              
              # ringi_closeby
              ringi_closeby  <- .ringi_closeby
              
              #'[01:Preparation of Neighboring Morphing Search]
              # While search bins
              temp.small.list <- list()
              
              # Save the found Big.Segment Summary
              temp.small.list[[1]] <- ringi_closeby
              
              # while-loop control
              # Compute and find small segments in between
              temp.small.i <- 2
              
              #'[02:Simplified tree ring structure from considered Edge] 
              # Parameters
              #' (1) .Degree_Control: How large is the angular range used;
              #' (2) .size.ctrl: How many points based on .Degree_Control.
              .Degree_Control = 10  # 10 means 10 degrees
              .size.ctrl <- .Degree_Control * ang.npts
              
              # Considered tree ring structure from ring.i
              #'[Add-in]: S.Flip
              .temp.total <- morph.list %>% dplyr::bind_rows()
              .temp.total <- rbind(.temp.total, ring.i)
              if(S.Flip == TRUE){
                
                temp.ring.i <- 
                  .temp.total %>% 
                  dplyr::filter(theta >= Edge.02 & theta < (Edge.02 + 0.5))
                
                # Size control
                if(nrow(temp.ring.i) > .size.ctrl){
                  temp.ring.i <- temp.ring.i[order(-temp.ring.i$theta),]
                  temp.ring.i <- temp.ring.i[c((nrow(temp.ring.i)-.size.ctrl):
                                                 nrow(temp.ring.i)),]
                }
                
              }else{
                
                temp.ring.i <- 
                  .temp.total %>% 
                  dplyr::filter(theta >= (Edge.01 - 0.5) & theta < Edge.02)
                
                # Size control
                if(nrow(temp.ring.i) > .size.ctrl){
                  temp.ring.i <- temp.ring.i[order(temp.ring.i$theta),]
                  temp.ring.i <- temp.ring.i[c((nrow(temp.ring.i)-.size.ctrl):
                                                 nrow(temp.ring.i)),]
                }
                
              }
              
              #'[03:While-loop for further close-by segments search]
              while (TRUE) {
                
                temp.section.sum <- temp.small.list %>% dplyr::bind_rows()
                temp.small.ID  <- temp.section.sum$clst
                temp.sec.small <- 
                  ring.ck %>% 
                  dplyr::filter(clst %in% temp.small.ID) %>% 
                  rbind(temp.ring.i)
                
                # spline degree of freedom decision
                .thres.ctrl = 0.5
                .npts.ctrl  = .thres.ctrl / ReSample.Size
                temp.small.range <- 
                  max(temp.section.sum$max.theta) - 
                  min(temp.section.sum$min.theta)
                temp.small.npts  <- temp.section.sum$npts %>% sum()
                if(temp.small.range > .thres.ctrl | 
                   temp.small.npts  > .npts.ctrl ){
                  temp.small.df <- 10}else{temp.small.df <- 5}
                
                # smooth.spline
                temp.small.smsp  <- smooth.spline(x  = temp.sec.small$theta,
                                                  y  = temp.sec.small$dist,
                                                  df = temp.small.df)
                temp.sp <-
                  predict(temp.small.smsp,
                          x = temp.sec.small$theta) %>%
                  as.data.frame()
                names(temp.sp) <- c("theta", "dist")
                temp.sp$theta  <- temp.sp$theta %>% round(digits = .digits)
                
                # Create Morphed tree ring structure 
                # between Edge and close-by segments accordingly (S.Flip)
                #'@note
                #' Search towards the other Edges far from close-by segments
                #' (1) S.Flip != TRUE: closed-by segments - Edge.02
                #' (2) S.Flip == TRUE: closed-by segments - Edge.01
                #' [Add-in]: S.Flip
                if(S.Flip != TRUE){
                  
                  # Morph in between Edge.02 and close-by
                  temp.sec.edge01 <- temp.sp$theta %>% max()
                  temp.sec.dist01 <- temp.sp$dist[which.max(temp.sp$theta)]
                  temp.rm <- Section.Morph(Edge.01  = temp.sec.edge01,
                                           Edge.02  = Edge.02,
                                           dist.01  = temp.sec.dist01,
                                           dist.02  = dist.02,
                                           ref.ring = ring.2)
                  
                  # Alignment of tree ring segments and morphed structure
                  temp.section.r <-
                    ring.ck %>% 
                    segment.filter(Edge.01 = temp.sec.edge01,
                                   Edge.02 = Edge.02,
                                   ring.segments = .) %>% 
                    ring.join(., ring.segment = temp.rm)
                  
                }else{
                  
                  # Morph in between Edge.01 and close-by
                  temp.sec.edge02 <- temp.sp$theta %>% min()
                  temp.sec.dist02 <- temp.sp$dist[which.min(temp.sp$theta)]
                  temp.rm <- Section.Morph(Edge.01  = Edge.01,
                                           Edge.02  = temp.sec.edge02,
                                           dist.01  = dist.01,
                                           dist.02  = temp.sec.dist02,
                                           ref.ring = ring.2)
                  
                  # Alignment of tree ring segments and morphed structure
                  temp.section.r <-
                    ring.ck %>% 
                    segment.filter(Edge.01 = Edge.01,
                                   Edge.02 = temp.sec.edge02,
                                   ring.segments = .) %>% 
                    ring.join(., ring.segment = temp.rm)
                  
                }
                
                head(temp.section.r)
                
                # Section.summary of in-between tree ring segments
                if(purrr::is_empty(temp.section.r$clst)){
                  
                  # Nothing in between
                  temp.section.sum <- temp.ring.sc.sum %>% as.data.frame()
                  temp.section.sum <- temp.section.sum[vector(),]
                  
                }else{
                  
                  # Something in between
                  temp.section.sum <-
                    temp.section.r %>% 
                    dplyr::mutate(., delta = abs(dist - dist.Morph)) %>% 
                    dplyr::group_by(clst) %>% 
                    dplyr::arrange(delta) %>% 
                    dplyr::summarise(npts     = n(),
                                     m.delta  = mean(delta[1:round(npts*0.1)]),
                                     sd.delta = sd  (delta[1:round(npts*0.1)]),
                                     # m.Ang.pos = npts / nrow(ring.2),
                                     max.theta = max(theta),
                                     min.theta = min(theta),
                                     dist.l    = dist[which.min(theta)],
                                     dist.r    = dist[which.max(theta)]) %>% 
                    dplyr::filter(., m.delta < 5) %>% 
                    dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                    # segment.clean(., dist.thres = 2) %>% 
                    segment.clean(dist.thres = 3) %>% 
                    segment.clean(dist.thres = 2, 
                                  segment.size = segment.size) %>% 
                    segment.overlap(., segment.size = segment.size) %>% 
                    as.data.frame()
                  
                }
                
                # Select segments that are closed to close-by segments
                temp.section.sum
                temp.section.sum <- 
                  temp.section.sum %>% 
                  segment.select(neighbor.thres = 0.1,
                                 rev = rev.control,
                                 dist.check = TRUE,
                                 dist.thres = 10)
                
                # Only when there are something close-by, then make sense
                if(purrr::is_empty(temp.section.sum$clst) != TRUE){
                  
                  # .range
                  .range <-
                    max(temp.section.sum$max.theta)-
                    min(temp.section.sum$min.theta)
                  
                  # Only when segments are close enough, then make sense
                  if(S.Flip != TRUE){
                    
                    #'[Conditions:]
                    #' Close-by (Theta)
                    cond.1 <- 
                      (temp.section.sum$min.theta[1] - temp.sec.edge01) > 0.1
                    #' Close-by (dist)
                    cond.2 <- 
                      abs(temp.section.sum$dist.l[1] - temp.sec.dist01) > 10
                    #' Noise
                    cond.3 <-
                      .range < 0.1 &
                      nrow(temp.section.sum) < 3 &
                      sum(temp.section.sum$segment.l) < 0.8*segment.size
                    if(cond.1 | cond.2 | cond.3){
                      temp.section.sum <- temp.section.sum[vector(),]
                    }
                    
                  }
                  if(S.Flip == TRUE){
                    
                    #'[Conditions:]
                    #' Close-by (Theta)
                    cond.1 <- 
                      (temp.sec.edge02 - temp.section.sum$max.theta[1]) > 0.1
                    #' Close-by (dist)
                    cond.2 <- 
                      abs(temp.section.sum$dist.r[1] - temp.sec.dist02) > 10
                    #' Noise
                    cond.3 <-
                      .range < 0.1 &
                      nrow(temp.section.sum) < 3 &
                      sum(temp.section.sum$segment.l) < 0.8*segment.size 
                    if( cond.1 | cond.2 | cond.3){
                      temp.section.sum <- temp.section.sum[vector(),]
                    }
                    
                  }
                }
                
                # while-loop control
                #'@note 
                #' The search will stop if the "temp.section.sum" is empty.
                #' Otherwise, the selected segments will be collected in bin,
                #' and another above morph-search will be processed again,
                #' considering all selected segments in bin.
                #' (1) temp.small.list: bin;
                #' (2) temp.small.i: while-loop step.
                if(purrr::is_empty(temp.section.sum$clst)){break}
                
                temp.small.list[[temp.small.i]] <- temp.section.sum
                temp.small.i <- temp.small.i + 1
                
                # Plot Check
                # plot(temp.sec.small$theta, temp.sec.small$dist)
                # lines(temp.small.smsp, col = "darkgreen")
                # points(ring.m$theta, 
                #        ring.m$dist, 
                #        cex = 0.3, col = "lightblue")
                
              }
              
              ringi_closeby <- temp.small.list %>% dplyr::bind_rows()
              
            }else{
              
              ringi_closeby <- .ringi_closeby[vector(),]
              
            }
            
            #'[Post-processing after close-by search]
            #'@description
            # The simplest way is to pass results to following codes
            # by updating the "temp.ring.sc.sum"
            #'[In_ringi_closeby]
            #' (1) If There are any "big-Segment",
            #'     it will be handled from that position forward,
            #'     by "big-Segment" search.
            #' (2) If There are no "big-Segment",
            #'     The ringi_closeby will be handled as "Small.Back".
            
            #'[01:Big-Segment Check]
            #' Only when ringi_closeby exist
            close.big <- FALSE
            if(purrr::is_empty(ringi_closeby$clst)){
              
              .temp.range <- 0
              
            }else{
              
              #'Big-Segment Check
              #'@note:
              #' The Goal is to handle the continuous segments 
              #' However, the "segment.select" does not have 
              #' the possibility to find out a jump in a continuous series
              #' Especially when there are "big-segments" involved.
              #' Hence, we can only keep small segments in close-by search. 
              #' (to make sure such big segment will be checked)
              if(any(ringi_closeby$segment.l > segment.size)){
                
                .pos <- which(ringi_closeby$segment.l > segment.size)
                .ringi_closeby <- ringi_closeby[c(1:(.pos[1])),]
                ringi_closeby  <- ringi_closeby[c(1:(.pos[1]-1)),]
                
                # Signal
                close.big <- TRUE
                
              }
              
              # Only when segments are big enough, then make sense
              # 1000 = 0.1 theta range
              # old: ringi_closeby$npts %>% sum() >= 1000
              .temp.range <- 
                max(ringi_closeby$max.theta) - min(ringi_closeby$min.theta)
              
            }
            
            #'[02:close-by segment check after previous adjustment]
            #'@description 
            #' If there are big-segments and therefore 
            #' yield really small "ringi_closeby" (small ".temp.range")
            #' We forward close-by table including the 1st big-segments
            #' to the following search: "Big-segment Search".
            #' Otherwise, the small close-by segment search finishes.
            #'@signal
            #' Small.Back
            #' temp.ring.sc.sum
            
            print(.temp.range)
            
            if(.temp.range >= 0.1 |
               ringi_closeby$segment.l %>% sum() >= segment.size){
              
              # Signal for "Small-Segment Search"
              Small.Back <- TRUE
              temp.ring.sc.sum <- ringi_closeby
              
              # Loop-control 
              break
              
            }else{
              
              #'@special
              #' If code ends here, means the first search yield nothing,
              #' there will be 2 possibilities:
              #' (1) Big-Segment exist;
              #' (2) Not Enough structural length 
              #'     (< segment.size | small .temp.range)
              #' Two different solutions:
              #' (1) Big-Segment exist: End, pass to "Big-Segment Search"
              #' (2) Not Enough structural length: 
              #'     Search from another .ringi_closeby.e* (03:LOOP-BACK)
              if(close.big == TRUE){
                
                message("# Attempt Failed ... Big-Segment")
                
                # Signal for "Big-Segment Search"
                Small.Back <- FALSE
                temp.ring.sc.sum <- .ringi_closeby
                
                # Loop-control 
                break
                
              }
              
              
              #'[03:LOOP-BACK]
              #'@special
              # (1) Check which search has performed, clear it.
              # (2) Re-set all control Signals, and
              #     search from the other Edge.
              #'@note
              #' The search always first choose the closest segments to Edges.
              #' either .ringi_closeby.e1 or .ringi_closeby.e2.
              #' Therefore, the loop-back and outer while structure
              #' is for changing considered close-by group when 
              #' the one, originally targeted, 
              #' ended up with big segments and failed.
              
              # Report
              message("# Attempt Failed ... Loop-back")
              
              if(S.Flip == TRUE){
                
                # Clear edge.02 search
                .ringi_closeby.e2 <- .ringi_closeby.e2[vector(),]
                
              }else{
                
                # Clear edge.01 search
                .ringi_closeby.e1 <- .ringi_closeby.e1[vector(),]
                
              }
              
              # Re-set control signals
              s.control   <- FALSE
              S.Flip      <- FALSE
              rev.control <- FALSE
              
            }
            
          }
          
          
          
          ## (03)_Summary Table Adjustments ----
          #'[00_03 Acquire "Big-Segment" Positions]
          #'Fix
          if(Tail.check != TRUE){
            big.segment.pos <- which(temp.ring.sc.sum$segment.l > segment.size)}
          
          #' [Add-in: temp Search-flip from edge.01 to edge.02]
          #'@note
          #' When the "Big-segment" is too far from the edge,
          #' It is less meaningful...
          #' Switch the search might help...
          if(purrr::is_empty(big.segment.pos)!=TRUE & Tail.check != TRUE){
            
            #'@description 
            #' If code comes in this section,
            #' This means the "close-by search" did not work.
            #' (Small.Back == FALSE)
            #' Hence,
            #' This section determines segment grouping direction.
            #' If No big-segments left, 
            #' the code will handle it by "small.middle search".
            #' ** small.middle search will be rough without constraints.
            #' (1) .gap.l: gap between Edge01 and temp.ring.sc.sum
            #' (2) .gap.r: gap between temp.ring.sc.sum and Edge02
            #'@signal
            #' (1) S.Flip: indicates the code to perform segment search
            #'             from Edge02 towards Edge01.
            .gap.l <- min(temp.ring.sc.sum$min.theta[big.segment.pos]) - Edge.01
            .gap.r <- Edge.02 - max(temp.ring.sc.sum$max.theta[big.segment.pos])
            S.Flip <- FALSE
            if(.gap.r < .gap.l){
              
              message("# Search-Flip from Tail #")
              S.Flip <- TRUE
              
              # Check
              print("# Section II")
              
              #'@note
              #' The Search for small segments from edges.01.02 are done...
              #' The S.Flip will ONLY be activated if "Small.Back" = FALSE
              #' Sequentially, S.Flip only for:
              #' (1) Big-segments close to "Edge.02"
              #' (2) Small-segments close to "Edge.02" and in-between.(later)
              #' Hence, 
              #' same logic should also apply constrain for (1).
              #' ".gap.r < .gap.l": big-segment is closer to Edge02.
              big.segment.pos <- tail(big.segment.pos, n = 1)
              
            }
            
          }else{
            
            # Signals
            .gap.l <- 0
            .gap.r <- 0
            
          }
          
          
          
          ## (04)_Middle Check ----
          #'[00_04_Middle Morph Check based on Big-segments]
          #'@description 
          #' If the processing step is taken into the "Middle Check Loop",
          #' this means that the "Big-Segment" is considered suspicious.
          #' The solution offer here is:
          #' (1) conducting a "Free" segment grouping
          #'     from such suspicious "Big-Segment" ,and then 
          #'     check whether it still finds the Edge.01 and Edge.02.
          #' (2) Normally, the "Reverse-Middle Search" is triggered to
          #'     enhance the accuracy of "Middle Search".
          #'     The "Reverse" means: when search starts from A segment,
          #'     yielded as vector(E, F, A, B, C, D);
          #'     depends on the current k-loop search direction, the 
          #'     reverse search starts from the other side of the 1st results:
          #'     [S.Flip == F]: D;
          #'     [S.Flip == T]: E; 
          #' (3) If the Middle Check is carried out, the corresponding ring.2
          #'     (reference ring base) will also be updated based on the middle
          #'     search results.
          #'@conditions:
          #' (1) The selected big.segment is "far" from both Edges 
          #'     (.gap.l, .gap.r)  > 0.1 ;
          #' (2) Based on the search direction of k-loop, the corresponding 
          #'     "distance to pith" has huge gap to the Edge;
          #'     (.gap.dist > 5)
          #' (3) Tail Check == TRUE.
          
          #'[Preparation]
          #'[01:Check Gap distance between Big-segment & corresponding Edge]
          #'@note ".gap.dist"
          #' This aims at checking whether the considered big-segment 
          #' has a huge distance (dist to pith) difference compared to 
          #' corresponding Edge (Check: S.Flip).
          if(purrr::is_empty(big.segment.pos)){
            
            .gap.dist <- 0
            
          }else{
            
            # Target Segment:
            .mid.big.sum <- temp.ring.sc.sum[big.segment.pos[1],]
            if(S.Flip == TRUE){
              
              .gap.dist <- abs(.mid.big.sum$dist.r - dist.02)
              
            }else{
              
              .gap.dist <- abs(dist.01 - .mid.big.sum$dist.l)
              
            }
            
          }
          
          #'[02:Ensure Data Structure]
          ring.i$theta <- ring.i$theta %>% round(digits = .digits)
          
          
          # Middle.Check ------------------------------------------------------#
          #'[Middle Check]
          #'@note
          #' Massive section: line 1497 - 4185
          #'@description 
          #' This process activated when:
          #' (1) theta gap (.gap.l, .gap.r) considered is too big;
          #' (2) dist gap (.gap.dist) considered is too big;
          #' (3) Tail.check == TRUE.
          if(min(.gap.l, .gap.r) > 0.1 |
             .gap.dist > 5 |
             Tail.check == TRUE){
            
            # Check
            # stop("Middle Check")
            
            #'[01:Tail Check Adjustments]
            if(Tail.check == TRUE){
              
              # Signal
              Tail.check <- FALSE
              
              # Report
              message("# Middle Check ... TAIL")
              
              .mid.big.sum <- 
                temp.ring.sc.sum[tail(big.segment.pos, n = 1),]
              
            }else{
              
              # Report
              message("# Middle Check ...")
              
              .mid.big.sum <- temp.ring.sc.sum[big.segment.pos[1],]
              
            }
            
            #'[02:Considered Segments Control by Reverse-search (S.Flip)]
            #'@description 
            #' This is to limit the free search range for "Middle Search".
            #' (1) ..ring.ck: ring.ck controls segments to be considered 
            #'     in between 2 Edges. "..ring.ck" is similar, but only 
            #'     consider in from Ring.ID.i towards corresponding Edge.
            #' (2) S.Flip: indicates whether it is Edge.01 or Edge.02 to be 
            #'     the relevant Edge.
            if(S.Flip == TRUE){
              ..ring.ck <-
                ring.ck %>% 
                dplyr::filter(clst != Ring.ID.i |
                                (clst == Ring.ID.i & theta > Edge.01))
            }else{
              ..ring.ck <-
                ring.ck %>% 
                dplyr::filter(clst != Ring.ID.i |
                                (clst == Ring.ID.i & theta < Edge.02))
            }
            
            #'[03:Consecutive While-loop Search]
            #'@description 
            #' Start the free search from suspicious Big-Segment.
            #' (1) morph.id: Targeted Tree ring segment ID*s*;
            #' (2) .mid.big: Targeted Tree ring segment structure data;
            #'@note
            #' *clst* Name in .mid.big & .mid.big.sum is *Unique*.
            #' It will only be the one that was targeted.
            mid.big.id   <- .mid.big.sum$clst
            mid.big.df   <- ring.ck %>% filter(clst %in% mid.big.id)
            npts.range   <- 0.2
            middle.check <- 
              segment.MiddleCheck(target.segment.df   = mid.big.df,
                                  considered.segments = ..ring.ck,
                                  ref.ring            = ring.2,
                                  npts.range          = npts.range,
                                  segment.SUM         = lookup.SUM,
                                  segment.size        = segment.size,
                                  ReSample.Size       = ReSample.Size)
            
            # Use Results directly from "middle.check"
            .temp.section.sum <- middle.check$temp.section.sum
            .mid.big.sum      <- middle.check$mid.big.sum
            morph.id          <- middle.check$morph.id
            .mid.big          <- middle.check$mid.big
            mid.ring.2        <- middle.check$mid.ring.2
            
            
            
            # Reverse.Middle.Check --------------------------------------------#
            #'[04:Reverse_Middle.check]
            #'@description
            #' To Ensure the accuracy of the Middle Search results. 
            
            #'[04_01:Decision of Reverse-Search]
            #'@signal *reverse.check*
            # preparation
            .range <- 
              middle.check$mid.big$theta %>% 
              range() %>% dist() %>% as.numeric()
            .dist  <- 
              middle.check$mid.big$dist %>% 
              range() %>% dist() %>% as.numeric()
            
            # Decision
            reverse.check <- FALSE
            if(.range <= 2 & .dist < 100){
              
              # Preparation of reverse search
              reverse.sum <- 
                middle.check$temp.section.sum %>% 
                dplyr::filter(clst %in% middle.check$morph.id) %>% 
                dplyr::filter(segment.l > segment.size)
              
              # S.Flip Adjustments
              if(S.Flip == FALSE){
                
                ..mid.big.sum <-
                  reverse.sum %>% 
                  dplyr::arrange(desc(min.theta)) %>% 
                  dplyr::slice(1)
                
              }else{
                
                ..mid.big.sum <-
                  reverse.sum %>% 
                  dplyr::arrange(min.theta) %>% 
                  dplyr::slice(1)
                
              }
              
              if(..mid.big.sum$clst != middle.check$mid.big.sum$clst){
                
                #'@note
                # No need to do reverse-check,
                # IF the targeted segment cluster is the same.
                
                # Signal
                reverse.check <- TRUE
                
              }
              
            }
            
            #'[04_02:Determination of Middle-Search Results....]
            #'@note
            #' (1) .temp.section.sum: grouped segments summary;
            #' (2) .mid.big.sum: grouped Middle-search segment summary;
            #' (3) morph.id: grouped segment IDs;
            #' (4) .mid.big: grouped Middle-search segment structure;
            #' (5) mid.ring.2: final Middle-Search reference base.
            if(reverse.check == TRUE){
              
              #'[_01:Conduct Reverse Middle-Search]
              #'@note
              #' Perform the Middle-search AGAIN from the searched neighbor in
              #' REVERSED-Direction based on "S.Flip" Signal.
              #' The aim is to ensure the searched segments did not perform a
              #' "Jump" to the neighboring tree ring structure.
              
              # In-put for fun "segment.MiddleCheck"
              mid.big.id <- ..mid.big.sum$clst
              mid.big.df <- ring.ck %>% filter(clst %in% mid.big.id)
              npts.range <- 0.2
              
              # Run Middle Check Again...
              middle.ckrev <- 
                segment.MiddleCheck(target.segment.df   = mid.big.df,
                                    considered.segments = ..ring.ck,
                                    ref.ring      = ring.2,
                                    npts.range    = npts.range,
                                    segment.SUM   = lookup.SUM,
                                    segment.size  = segment.size,
                                    ReSample.Size = ReSample.Size)
              
              #'[_02:Decision of Compiling Reverse Middle-Search Results]
              #'@import 
              #' (1) middle.check
              #' (2) middle.ckrev
              #'@eval
              #' (1) Both (middle.check, middle.ckrev) results are the same;
              #' (2) reversed-search found the middle.check targeted ID;
              #' (3) all reversed-search were in the result of middle.check;
              #' (4) A jump has happened in middle.check.
              #'@export
              #' (1) result.update == TRUE: *morph.id*;
              #' (2) result.update == FALSE: Results from middle.check
              
              # Signal-off
              reverse.check <- FALSE
              result.update <- FALSE
              if(identical(sort(middle.ckrev$morph.id), 
                           sort(middle.check$morph.id))){
                
                # Remain Middle Check Results from "middle.check"
                # Signal-off
                result.update <- FALSE
                
              }else if(middle.check$mid.big.sum$clst %in% 
                       middle.ckrev$morph.id){
                
                # Combine the search product "morph.id"
                morph.id <- append(middle.ckrev$morph.id, 
                                   middle.check$morph.id) %>% unique()
                
                # Signal
                result.update <- TRUE
                
              }else if(all(middle.ckrev$morph.id %in% 
                           middle.check$morph.id)){
                
                # Remain Middle Check Results from "middle.check"
                # Signal-off
                result.update <- FALSE
                
              }else{
                
                #'@special
                #' This means that the search from "middle.check"
                #' performs a *JUMP* to reverse-searched results "middle.ckrev"
                #' Hence, we keep the reverse searched segment.(*middle.ckrev*)
                
                # Preparation of reverse search Decisions
                #'[Middle-Search]
                save.sum <- 
                  middle.check$temp.section.sum %>% 
                  dplyr::filter(clst %in% middle.check$morph.id)
                save.diff.id <- 
                  dplyr::setdiff(middle.check$morph.id, middle.ckrev$morph.id)
                save.sum.diff <- 
                  middle.check$temp.section.sum %>% 
                  dplyr::filter(clst %in% save.diff.id)
                
                #'[Targeted Segment for Middle Check]
                save.id.sum <-
                  save.sum %>% 
                  dplyr::filter(clst %in% middle.check$mid.big.sum$clst)
                
                #'[Reversed_Middle-Search]
                rev.sum <-
                  middle.ckrev$temp.section.sum %>% 
                  dplyr::filter(clst %in% middle.ckrev$morph.id)
                rev.diff.id <- 
                  dplyr::setdiff(middle.ckrev$morph.id, middle.check$morph.id)
                rev.sum.diff <-
                  middle.ckrev$temp.section.sum %>% 
                  dplyr::filter(clst %in% rev.diff.id)
                
                #'[Intersected]
                # overlap.id <- dplyr::intersect(morph.id, save.morph.id)
                # overlap.sum <-
                #   temp.section.sum %>% 
                #   dplyr::filter(clst %in% overlap.id)
                
                #'[Preparation]
                rev.range <- 
                  middle.ckrev$mid.big$theta %>% 
                  range() %>% dist() %>% as.numeric()
                rev.dist  <- 
                  middle.ckrev$mid.big$dist %>% 
                  range() %>% dist() %>% as.numeric()
                
                #'[Condition Check]
                #'@note
                #' (0) cond.0: reverse results are large (1st reverse if def.);
                #' (1) cond.1: both small, but reverse one is larger;
                #' (2) cond.2: Non-flip cond., rev. results close to Edge.01;
                #' (3) cond.3: Flip cond., reversed results close to Edge.02
                cond.0 <- rev.range > 2 | rev.dist > 100
                cond.1 <- 
                  sum(save.sum.diff$npts) < 1000 & 
                  sum(save.sum.diff$npts) < sum(rev.sum.diff$npts)
                # Fix
                cond.2 <- FALSE
                cond.3 <- FALSE
                if(purrr::is_empty(rev.sum.diff$clst)!=TRUE){
                  cond.2 <-
                    S.Flip == FALSE &
                    min(rev.sum.diff$min.theta) < save.id.sum$max.theta
                  cond.3 <-
                    S.Flip == TRUE &
                    max(rev.sum.diff$min.theta) > save.id.sum$min.theta
                }
                
                if(cond.0 | cond.1 | cond.2 | cond.3){
                  
                  # stop("Condition Check ... Reverse-Search")
                  
                  # Signal Report
                  message("# Middle-Check... Reverse-Results Acquired.")
                  result.update <- FALSE
                  
                  # Use Results directly from "middle.check"
                  .temp.section.sum <- middle.ckrev$temp.section.sum
                  .mid.big.sum      <- middle.ckrev$mid.big.sum
                  morph.id          <- middle.ckrev$morph.id
                  .mid.big          <- middle.ckrev$mid.big
                  mid.ring.2        <- middle.ckrev$mid.ring.2
                  
                }else{
                  
                  #' [Update morph.id]
                  #' Searched Results without segments from reversed search.
                  morph.id <- save.diff.id
                  
                  # Signal Report
                  message("# Middle-Check... Reverse-Results Dedined.")
                  result.update <- TRUE
                  
                }
                
              }
              
              #'[_03:Update Reverse Middle-Search Results]
              #'@note
              #' Decision & modification of *morph.id* controls how to "update"
              #' relevant outputs from "Middle Check".
              if(result.update == TRUE){
                
                #' [Update .mid.big]
                .mid.sp.df <- ..ring.ck %>% dplyr::filter(clst %in% morph.id)
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
                temp.mid.sp$theta  <- 
                  temp.mid.sp$theta %>% round(digits = .digits)
                
                # Check
                # plot(.mid.sp.df$theta, .mid.sp.df$dist, cex = 0.3)
                # lines(.mid.sp, col = "darkgreen")
                
                # Morph & update the loop-control
                # Update .mid.big
                .mid.big      <- temp.mid.sp
                .mid.big$clst <- .mid.big.sum$clst
                
                #' [Update .mid.big.sum]
                .mid.big.sum$max.theta <- max(.mid.big$theta)
                .mid.big.sum$min.theta <- min(.mid.big$theta)
                .mid.big.sum$dist.r <- .mid.big$dist[which.max(.mid.big$theta)]
                .mid.big.sum$dist.l <- .mid.big$dist[which.min(.mid.big$theta)]
                
                #' [Update mid.ring.2]
                mid.edge.01 <- .mid.big.sum$min.theta
                mid.edge.02 <- .mid.big.sum$max.theta
                mid.dist.01 <- .mid.big.sum$dist.l
                mid.dist.02 <- .mid.big.sum$dist.r
                
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
                
              }
              
            }
            
            #'[Plot Check]
            if(.plotcheck == TRUE){
              
              points(.mid.big$theta, .mid.big$dist, cex = 0.3, col = "darkred")
              
            }
            
            # Re-write the ring.2 -------------------------------------------- #
            #'[05:Update the reference tree ring structure "ring.2"]
            #'@note
            #' This ONLY Section Happens when Middle Check if activated.
            #' As the Middle check provide the morph-search from segment,
            #' its *path of search (mid.ring.2)* indicated more Robust idea of 
            #' how the tree ring structure at current location should be.
            
            # Dev.Check
            # stop("check ring.2 morph")
            
            #'[05_01:Ring.i Search-Check]
            #'@export
            #' (1) ...ring.i
            #' (2) *..ring.keep*: segments kept by k-loop within specific range*
            #'@note
            #' Fix: ring.i based on "S.Flip"
            #' As while-loop has search function from Edge.02 (S.Flip = T),
            #' the middle-check also adapts such idea.
            #' Noted that in "#' Fix:...", 
            #' as sometimes Middle-search Result (.mid.big) is large,
            #' eventually meeting the other Edges accidentally, 
            #' the influence on ...ring.keep.sum by padded ring.i is handled.
            
            ...ring.i <-
              ring.i %>% 
              dplyr::filter(clst != "S.Flip") %>% 
              {if(S.Flip != TRUE){dplyr::filter(., theta < (Edge.02 - .M.Step))
              }else{dplyr::filter(., theta > (Edge.01 + .M.Step))}}
            
            ..ring.keep <-
              ring.keep %>% 
              dplyr::bind_rows() %>% 
              rbind(...ring.i) %>% 
              dplyr::filter(.$theta < (max(.mid.big$theta) + 1.5)) %>%
              dplyr::filter(.$theta > (min(.mid.big$theta) - 1.5))
            
            #' Fix: Morph-Base has lacking one side but not TRUE
            if(purrr::is_empty(..ring.keep$clst)!=TRUE){
              
              if(Edge.01 > (min(.mid.big$theta) - 1.5) &
                 min(..ring.keep$theta) > max(.mid.big$theta)){
                
                # Fix
                if(identical(unique(..ring.keep$clst), Ring.ID.i)){
                  # Adjust cluster name
                  ..ring.keep$clst <- Ring.i.tempclst #"Ring.i"
                  .pos = which(lookup.SUM$clst == Ring.ID.i)
                  .ring.i.lookup <- 
                    data.frame(clst = Ring.i.tempclst, 
                               segment.l = lookup.SUM$segment.l[.pos])
                  lookup.SUM <- rbind(lookup.SUM, .ring.i.lookup)
                }
                
                ...ring.i <- 
                  ring.i %>% 
                  dplyr::filter(clst != "S.Flip") %>% 
                  dplyr::filter(theta < (Edge.02 - .M.Step))
                
                ..ring.keep <- 
                  ..ring.keep %>% 
                  dplyr::filter(clst != Ring.ID.i) %>% 
                  rbind(...ring.i)
                
              }
              if(Edge.02 < (max(.mid.big$theta) + 1.5) &
                 max(..ring.keep$theta) < min(.mid.big$theta)){
                
                # Fix
                if(identical(unique(..ring.keep$clst), Ring.ID.i)){
                  # Adjust cluster name
                  ..ring.keep$clst <- Ring.i.tempclst  #"Ring.i"
                  .pos = which(lookup.SUM$clst == Ring.ID.i)
                  .ring.i.lookup <- 
                    data.frame(clst = Ring.i.tempclst, 
                               segment.l = lookup.SUM$segment.l[.pos])
                  lookup.SUM <- rbind(lookup.SUM, .ring.i.lookup)
                }
                
                ...ring.i <- 
                  ring.i %>% 
                  dplyr::filter(clst != "S.Flip") %>% 
                  dplyr::filter(theta > (Edge.01 + .M.Step))
                
                ..ring.keep <- 
                  ..ring.keep %>% 
                  dplyr::filter(clst != Ring.ID.i) %>% 
                  rbind(...ring.i)
                
              }
              
            }
            
            # Ensure Data Structure
            ..ring.keep$clst <- ..ring.keep$clst %>% as.integer()
            head(..ring.keep)
            
            # Check
            # points(..ring.keep$theta,
            #        ..ring.keep$dist,
            #        cex = 0.5, col = "red3")
            # points(mid.ring.2$theta,
            #        mid.ring.2$dist,
            #        cex = 0.3, col = "lightblue")
            
            #'[05_02:Segments from ring.keep to structure new reference]
            #'@export
            #' (1) ..keep.sum.list -> *...ring.keep.sum*
            #' (2) ..keep.sum.value
            #'@note
            #' Align "ring.keep" by reference structure (mid.ring.2),
            #' update the metrics of the summary based on mid.ring.2.
            #'@details 
            #' In this section, segments in "ring.keep" were used to 
            #' spline the updated structure to 
            #' merge the original ring.2 & mid.ring.2. 
            
            ...ring.keep <- 
              ..ring.keep %>% 
              ring.join(., ring.segment = mid.ring.2) %>% 
              dplyr::mutate(., delta = abs(dist - dist.ref))
            head(...ring.keep)
            
            if(purrr::is_empty(...ring.keep$clst)){
              
              # Nothing in between
              ..ring.keep.sum <- temp.ring.sc.sum %>% as.data.frame()
              ..ring.keep.sum <- ..ring.keep.sum[vector(),]
              
              # Fix
              ..keep.sum.ref <- data.frame(clst = "FAKE",
                                           npts = 0,
                                           m.delta   = Inf,
                                           sd.delta  = 0,
                                           max.theta = 0,
                                           min.theta = 0,
                                           dist.l    = 0,
                                           dist.r    = 0)
              
              ..keep.sum.list  <- list()
              ..keep.sum.list[[1]] <- ..ring.keep.sum[vector(),]
              ..keep.sum.list[[2]] <- ..ring.keep.sum[vector(),]
              ..keep.sum.value <- c(Inf, Inf)
              
            }else{
              
              # Update ..ring.keep.sum
              ..ring.keep.sum <-
                ...ring.keep %>%
                dplyr::group_by(clst) %>% 
                dplyr::arrange(delta) %>% 
                dplyr::summarise(npts   = n(),
                                 m.delta=mean(delta[1:round(npts.range*npts)]),
                                 sd.delta= sd(delta[1:round(npts.range*npts)]),
                                 max.theta = max(theta),
                                 min.theta = min(theta),
                                 dist.l    = dist[which.min(theta)],
                                 dist.r    = dist[which.max(theta)]) %>% 
                dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                dplyr::arrange(min.theta) %>% 
                as.data.frame()
              
              # Fix npts might == 1
              ..fix.pos <- which(..ring.keep.sum$npts == 1)
              ..ring.keep.sum$sd.delta[..fix.pos] <- 0
              
              # Scenarios of ring.keep for adjustments:
              #'@special
              #' The ring.keep have 3 different scenarios with mid.big:
              #' (1) None over-lap;
              #' (2) Over-lap 1 side, we need such position;
              #' (3) Over-lap 2 sides.
              
              # Save original data (Backup)
              ...ring.keep.sum <- ..ring.keep.sum
              
              # While-loop elements:
              ..keep.sum.list  <- list()
              ..keep.sum.value <- numeric()
              .temp.keep.sum   <- ..ring.keep.sum[vector(),]
              while (TRUE) {
                
                # Acquire close-by segments from both sides
                ..keep.sum.ref.e1 <-
                  ...ring.keep.sum %>% 
                  dplyr::filter(min.theta > min(max(.mid.big$theta), 
                                                (Edge.02 - .M.Step))) %>% 
                  segment.select(neighbor.thres = 0.1,
                                 rev = FALSE,
                                 dist.check = TRUE,
                                 dist.thres = 5)
                
                ..keep.sum.ref.e2 <-
                  ...ring.keep.sum %>% 
                  dplyr::filter(max.theta < max(min(.mid.big$theta), 
                                                (Edge.01 + .M.Step))) %>% 
                  segment.select(neighbor.thres = 0.1,
                                 rev = TRUE,
                                 dist.check = TRUE,
                                 dist.thres = 5)
                
                # Direction Check: e1 & e2
                if(purrr::is_empty(..keep.sum.ref.e1$clst) != TRUE){
                  ..keep.sum.ref <- ..keep.sum.ref.e1
                }else{..keep.sum.ref <- ..keep.sum.ref.e2}
                
                # loop-control
                if(purrr::is_empty(..keep.sum.ref$clst)){
                  
                  print("# Kept Segments (L) EMPTY...")
                  
                  # (1) Output Summary
                  if(purrr::is_empty(.temp.keep.sum$clst) != TRUE){
                    
                    ..keep.sum.list[[1]] <- .temp.keep.sum
                    
                    ..keep.sum.list[[2]] <-
                      ..ring.keep.sum %>% 
                      dplyr::filter(max.theta < max(min(.mid.big$theta), 
                                                    (Edge.01 + .M.Step))) %>% 
                      segment.select(neighbor.thres = 0.1,
                                     rev = TRUE,
                                     dist.check = TRUE,
                                     dist.thres = 5)
                    
                  }else{
                    
                    ..keep.sum.list[[1]] <-
                      ..ring.keep.sum %>% 
                      dplyr::filter(min.theta > min(max(.mid.big$theta), 
                                                    (Edge.02 - .M.Step))) %>% 
                      segment.select(neighbor.thres = 0.1,
                                     rev = FALSE,
                                     dist.check = TRUE,
                                     dist.thres = 5)
                    
                    ..keep.sum.list[[2]] <-
                      ..ring.keep.sum %>% 
                      dplyr::filter(max.theta < max(min(.mid.big$theta), 
                                                    (Edge.01 + .M.Step))) %>% 
                      segment.select(neighbor.thres = 0.1,
                                     rev = TRUE,
                                     dist.check = TRUE,
                                     dist.thres = 5)
                    
                  }
                  
                  # loop-control
                  break
                  
                }
                
                # Check ..keep.sum.ref
                .range <- 
                  max(..keep.sum.ref$max.theta) - 
                  min(..keep.sum.ref$min.theta)
                if(sum(..keep.sum.ref$segment.l) < 0.8*segment.size |
                   .range < 0.2){
                  
                  # Adjust ...ring.keep.sum (While-loop)
                  ..clst.keep <- dplyr::setdiff(...ring.keep.sum$clst,
                                                ..keep.sum.ref  $clst)
                  ...ring.keep.sum <- 
                    ...ring.keep.sum %>% 
                    dplyr::filter(clst %in% ..clst.keep)
                  
                }else{
                  
                  if(purrr::is_empty(.temp.keep.sum$clst)){
                    
                    if(identical(..keep.sum.ref.e1, ..keep.sum.ref)){
                      
                      print("# Kept Segments (R) Found...")
                      ..search.rb <- 
                        max(min(.mid.big$theta), (Edge.01 + .M.Step))
                      
                      if(purrr::is_empty(..keep.sum.ref.e2$clst) != TRUE){
                        
                        # While-loop update
                        .temp.keep.sum <- ..keep.sum.ref
                        ...ring.keep.sum <- 
                          ...ring.keep.sum %>% 
                          dplyr::filter(max.theta < ..search.rb)
                        
                      }else{
                        
                        print("# Kept Segments (L) EMPTY...")
                        print("# Search ENDs...")
                        
                        ..keep.sum.list[[1]] <- ..keep.sum.ref
                        
                        ..keep.sum.list[[2]] <-
                          ..ring.keep.sum %>% 
                          dplyr::filter(max.theta < ..search.rb) %>% 
                          segment.select(neighbor.thres = 0.1,
                                         rev = TRUE,
                                         dist.check = TRUE,
                                         dist.thres = 5)
                        # loop-control
                        break
                        
                      }
                      
                    }else{
                      
                      # Report
                      print("# Kept Segments (R) EMPTY...")
                      
                    }
                    
                    if(identical(..keep.sum.ref.e2, ..keep.sum.ref)){
                      
                      # Report
                      print("# Kept Segments (L) Found...")
                      ..search.lb <- 
                        min(max(.mid.big$theta), (Edge.02 - .M.Step))
                      
                      if(purrr::is_empty(..keep.sum.ref.e1$clst) != TRUE){
                        
                        # .temp.keep.sum <- ..keep.sum.ref
                        # ...ring.keep.sum <- 
                        #   ...ring.keep.sum %>% 
                        #   dplyr::filter(min.theta > max(.mid.big$theta))
                        stop("not possible ... check while-loop")
                        
                      }else{
                        
                        # Report
                        print("# Search ENDs...")
                        
                        ..keep.sum.list[[1]] <-
                          ..ring.keep.sum %>% 
                          dplyr::filter(min.theta > ..search.lb) %>% 
                          segment.select(neighbor.thres = 0.1,
                                         rev = FALSE,
                                         dist.check = TRUE,
                                         dist.thres = 5)
                        
                        ..keep.sum.list[[2]] <- ..keep.sum.ref
                        
                        # loop-control
                        break
                        
                      }
                      
                    }
                    
                  }else{
                    
                    # Report
                    print("# Kept Segments (L) Found...")
                    print("# Search ENDs...")
                    
                    ..keep.sum.list[[1]] <- .temp.keep.sum
                    
                    ..keep.sum.list[[2]] <- ..keep.sum.ref
                    
                    
                    # loop-control
                    break
                    
                  }
                  
                }
                
              }
              
              # (1) Output mean.delta indicators
              # (difference between segments and reference tree ring structure)
              if(is_empty(..keep.sum.list[[1]]$clst)){
                ..keep.sum.value[1] <- Inf
              }else{
                
                #'@special
                .pos <- which.min(..keep.sum.list[[1]]$min.theta)
                if(..keep.sum.list[[1]]$m.delta[.pos] < 2 &
                   ..keep.sum.list[[1]]$segment.l[.pos] > 0.8*segment.size){
                  
                  ..keep.sum.value[1] <- ..keep.sum.list[[1]]$m.delta[.pos]
                  
                }else{
                  
                  ..rb    <- min(..keep.sum.list[[1]]$min.theta) + 0.2
                  ..clst  <- ..keep.sum.list[[1]]$clst
                  ..theta <- ..keep.sum.list[[1]]$min.theta 
                  ..temp.keep.1 <- 
                    ...ring.keep %>% 
                    dplyr::filter(clst %in% ..clst) %>% 
                    dplyr::filter(theta < ..rb) %>% 
                    dplyr::filter(theta > min(..theta)) %>%
                    dplyr::mutate(delta = abs(dist - dist.ref)) %>% 
                    dplyr::arrange(delta) %>% 
                    dplyr::slice(c(1:round(npts.range*n())))
                  ..keep.sum.value[1] <- ..temp.keep.1$delta %>% mean()
                  
                }
                
              }
              if(is_empty(..keep.sum.list[[2]]$clst)){
                ..keep.sum.value[2] <- Inf
              }else{
                
                #'@special
                .pos <- which.max(..keep.sum.list[[2]]$min.theta)
                if(..keep.sum.list[[2]]$m.delta[.pos] < 2 &
                   ..keep.sum.list[[2]]$segment.l[.pos] > 0.8*segment.size){
                  
                  ..keep.sum.value[2] <- ..keep.sum.list[[2]]$m.delta[.pos]
                  
                }else{
                  
                  ..lb    <- max(..keep.sum.list[[2]]$max.theta) - 0.2
                  ..clst  <- ..keep.sum.list[[2]]$clst
                  ..theta <- ..keep.sum.list[[2]]$max.theta
                  ..temp.keep.2 <- 
                    ...ring.keep %>% 
                    dplyr::filter(clst %in% ..clst) %>% 
                    dplyr::filter(theta > ..lb) %>% 
                    dplyr::filter(theta < max(..theta)) %>%
                    dplyr::mutate(delta = abs(dist - dist.ref)) %>% 
                    dplyr::arrange(delta) %>% 
                    dplyr::slice(c(1:round(npts.range*n())))
                  ..keep.sum.value[2] <- ..temp.keep.2$delta %>% mean()
                  
                }
                
              }
              
              # Check
              ..keep.sum.list
              ..keep.sum.value
              
            }
            
            ...ring.keep.sum <- ..keep.sum.list %>% dplyr::bind_rows()
            
            #' Revision of the "Ring.ID.i"
            #'@special
            #' In case there is *ONLY 1 row* "...ring.keep.sum" available &
            #' row ID == *Ring.ID.i*, increase the acceptance criteria by 5.
            if(purrr::is_empty(...ring.keep.sum$clst) != TRUE){
              
              if(Ring.ID.i %in% ...ring.keep.sum$clst){
                
                .ring.i.pos <- which(...ring.keep.sum$clst == Ring.ID.i)
                
                if(nrow(...ring.keep.sum) == 1){
                  
                  .pos <- which(..keep.sum.value < Inf)
                  if(...ring.keep.sum$m.delta < 5){
                    
                    # Report
                    message("# Ring.ID.i Adjustments...")
                    print("Larger allowance on at the beginning of ring.i ")
                    ...ring.keep.sum$m.delta <- 0
                    ..keep.sum.value[.pos]   <- 0
                    
                  }
                  
                }else if(...ring.keep.sum$npts[.ring.i.pos] < 3){
                  
                  ...ring.keep.sum$sd.delta[.ring.i.pos] <- 0
                  ..keep.sum.value[is.na(..keep.sum.value)] <- 0
                  
                }
                
                
              }
              
            }else if(purrr::is_empty(...ring.keep.sum$clst) &
                     nrow(..ring.keep.sum) == 1){
              
              #'*Normally, this should not happen*
              #' The section means the given summary is the one without clst ID.
              #' The origin of such summary should be traced before continuing.
              #' Additionally, the possible solutions is prepared below...
              
              # Report
              message("# Ring.ID.i Adjustments...")
              print(" Large Middle.Search-Base Covers both Edges")
              
              stop("Check the Situation ...")
              
              if(S.Flip == TRUE){
                
                ..keep.sum.list[[1]] <-
                  ..ring.keep %>% 
                  dplyr::mutate(., delta = abs(dist - dist.ref)) %>% 
                  dplyr::filter(.$theta < (max(.mid.big$theta) + 1.5)) %>%
                  dplyr::filter(.$theta > (Edge.01 + .M.Step)) %>% 
                  dplyr::arrange(delta) %>% 
                  dplyr::group_by(clst) %>% 
                  dplyr::summarise(npts    = n(),
                                   ..npts  = round(npts.range*npts),
                                   m.delta = mean(delta[1:..npts]),
                                   sd.delta  = sd(delta[1:..npts]),
                                   max.theta = max(theta),
                                   min.theta = min(theta),
                                   dist.l    = dist[which.min(theta)],
                                   dist.r    = dist[which.max(theta)]) %>%
                  dplyr::select(-..npts) %>%  
                  dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                  dplyr::arrange(min.theta) %>% 
                  as.data.frame()
                
                ..keep.sum.value[1] <- ..keep.sum.list[[1]]$m.delta
                
              }else{
                
                ..keep.sum.list[[2]] <-
                  ..ring.keep %>% 
                  dplyr::mutate(., delta = abs(dist - dist.ref)) %>% 
                  dplyr::filter(.$theta < (Edge.02 - .M.Step)) %>%
                  dplyr::filter(.$theta > (min(.mid.big$theta) - 1.5)) %>% 
                  dplyr::arrange(delta) %>% 
                  dplyr::group_by(clst) %>% 
                  dplyr::summarise(npts    = n(),
                                   ..npts  = round(npts.range*npts),
                                   m.delta = mean(delta[1:..npts]),
                                   sd.delta  = sd(delta[1:..npts]),
                                   max.theta = max(theta),
                                   min.theta = min(theta),
                                   dist.l    = dist[which.min(theta)],
                                   dist.r    = dist[which.max(theta)]) %>%
                  dplyr::select(-..npts) %>%  
                  dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                  dplyr::arrange(min.theta) %>% 
                  as.data.frame()
                
                ..keep.sum.value[2] <- ..keep.sum.list[[2]]$m.delta
                
              }
              
              # Fix
              ...ring.keep.sum <- ..keep.sum.list %>% dplyr::bind_rows()
              
            }
            
            # Check
            ..ring.keep.sum
            ...ring.keep.sum
            
            #'[05_03:Segment structure data by ...ring.keep.sum]
            #'@note
            #' (1) For some reason, to extract structure data 
            #'     from ..ring.keep by ...ring.keep.sum, as ..ring.keep 
            #'     is occupied & idea is to provide a backbone structure to 
            #'     merge with middle.check's reference structure (mid.ring.2),
            #'     *..ring.i* is used instead of "...ring.keep".
            #' (2) ..ring.i.sum: segment components within .mid.big that 
            #'     overlapped with considered ring.keep base (..ring.keep). 
            
            ..ring.i <- 
              ..ring.keep %>% 
              dplyr::filter(clst %in% ...ring.keep.sum$clst)
            
            if(is_empty(..ring.i$clst)){
              
              ..ring.i.sum <- .mid.big[vector(),]
              
            }else{
              
              ..ring.i$clst <- "ref"
              
              #'@note
              #' Make Sure "Delta" = ring.keep/ring.i - .mid.big
              ..ring.i.sum <-
                .mid.big %>% 
                ring.join(., ..ring.i) %>% 
                dplyr::filter(is.na(dist.ref) != TRUE) %>% 
                dplyr::mutate(delta = dist.ref - dist)
              
            }
            
            # Check
            # points(..ring.i.sum$theta, 
            #        ..ring.i.sum$dist, cex = 0.5, col = "red")
            
            #'[05_04:Theta-Gap between .mid.big & ring.keep]
            #'@note
            #' (1) Negative gap values mean overlap with ..ring.i (ring.keep)
            #' (2) ..keep.sum.list: indicates theta gap's directional info.
            if(purrr::is_empty(..keep.sum.list[[1]]$clst)){
              ..keep.e1.lb <- Edge.02
            }else{
              ..keep.e1.lb <- min(..keep.sum.list[[1]]$min.theta)
            }
            if(purrr::is_empty(..keep.sum.list[[2]]$clst)){
              ..keep.e2.rb <- Edge.01
            }else{
              ..keep.e2.rb <- max(..keep.sum.list[[2]]$max.theta)
            }
            .gap.l <- min(.mid.big$theta) - ..keep.e2.rb
            .gap.r <- ..keep.e1.lb - max(.mid.big$theta)
            
            #'[05_05:Decision01]
            #'@description 
            #' (1) Differences to Edges < 3 pixels;
            #' (2) Theta to either edges are "close (<0.5 radius)".
            #' (3) If ..ring.i.sum is *EMPTY*, 
            #'     this means .mid.big locates within Edge.01 & Edge.02.
            #'     (Nothing Overlapped with-in and with Ring.i)
            #'@note
            #' Decision-aim here is to target structure differences < 3 &
            #' close to either edges from ring.keep.
            #' Noted that, theta < 0.5 radius includes "Negative" values,
            #' indicating co-agreements with ring.keep "might" be the case.
            #'@conditions
            #' (1) middle.I  :Small radius and theta differences;
            #' (2) middle.II :
            #' (3) middle.III:Same Segments were found by updated ring.2.
            
            middle.I   <- FALSE 
            middle.II  <- FALSE 
            middle.III <- FALSE 
            small.mid  <- FALSE
            if(min(..keep.sum.value) < 3){ # min(.gap.l, .gap.r) < 0.5
              
              message("# Middle-Search Aligned with kept segments")
              print(paste(c("... Segments :",
                            ...ring.keep.sum$clst), collapse = " "))
              print(..keep.sum.list)
              print(..keep.sum.value)
              middle.I  <- TRUE
              
            }else if(purrr::is_empty(..ring.i.sum$clst)){
              
              #'@note
              #' ..ring.i.sum *EMPTY*, means the middle search does not
              #' overlap with the ring.keep segments, yielded by
              #' (1) Strong segmentation, they *STILL* belong each other;
              #' (2) *Disagreements*, we use *mid.ring2* 
              #'     to update reference ring structure (ring.2) 
              #'     to avoid a jump to neighboring tree rings.
              
              message("# Reshape reference ring at ring.ID ",
                      .mid.big.sum$clst, " position")
              
              # Ensure data structure
              ring.2 <- 
                ring.2 %>% 
                dplyr::arrange(theta) %>% 
                dplyr::mutate(theta = theta %>% round(digits = .digits))
              
              .pos <- which(ring.2$theta > (min(.mid.big$theta) - .M.Step) &
                              ring.2$theta < (max(.mid.big$theta) + .M.Step))
              .Edge.01 = 
                (min(.mid.big$theta) - ReSample.Size) %>% 
                round(., digits = .digits)
              .Edge.02 = 
                (max(.mid.big$theta) + ReSample.Size) %>% 
                round(., digits = .digits)
              .dist.01 = ring.2$dist[(range(.pos)[1] - 1)]
              .dist.02 = ring.2$dist[(range(.pos)[2] + 1)]
              .morph.ref <- Section.Morph(Edge.01  = .Edge.01,
                                          Edge.02  = .Edge.02,
                                          dist.01  = .dist.01,
                                          dist.02  = .dist.02,
                                          ref.ring = .mid.big)
              
              ring.2$dist[.pos] <- .morph.ref$dist
              ring.2$clst <- "ref"
              
              # Ensure ring.2 Structure
              ring.2$theta <- ring.2$theta %>% round(digits = .digits)
              
              temp.rm <- Section.Morph(Edge.01  = Edge.01,
                                       Edge.02  = Edge.02,
                                       dist.01  = dist.01,
                                       dist.02  = dist.02,
                                       ref.ring = ring.2)
              
              # Update searched segments by updated ring.2 & morphed structure
              print("Updating temp.ring.sum & big.segment.pos")
              
              temp.section.r <-
                ring.ck %>% 
                segment.filter(Edge.01 = Edge.01,
                               Edge.02 = Edge.02,
                               ring.segments = .) %>% 
                ring.join(., ring.segment = temp.rm)
              
              head(temp.section.r)
              
              if(purrr::is_empty(temp.section.r$clst)){
                
                # Nothing in between
                temp.ring.sc.sum <- temp.ring.sc.sum %>% as.data.frame()
                temp.ring.sc.sum <- temp.section.sum[vector(),]
                
              }else{
                
                # Something in between
                temp.ring.sc.sum <-
                  temp.section.r %>% 
                  dplyr::mutate(., delta = abs(dist - dist.Morph)) %>% 
                  dplyr::group_by(clst) %>% 
                  dplyr::arrange(delta) %>% 
                  dplyr::summarise(npts    = n(),
                                   ..npts  = round(npts.range*npts),
                                   m.delta = mean(delta[1:..npts]),
                                   sd.delta  = sd(delta[1:..npts]),
                                   max.theta = max(theta),
                                   min.theta = min(theta),
                                   dist.l    = dist[which.min(theta)],
                                   dist.r    = dist[which.max(theta)]) %>% 
                  dplyr::select(-..npts) %>% 
                  dplyr::filter(., m.delta < 5) %>% 
                  dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                  # segment.clean(., dist.thres = 2) %>% 
                  segment.clean(dist.thres = 3) %>% 
                  segment.clean(dist.thres = 2, 
                                segment.size = segment.size) %>% 
                  segment.overlap(., segment.size = segment.size) %>% 
                  as.data.frame()
                
              }
              
              # Check
              temp.ring.sc.sum
              
              # Update big.segment.pos
              big.segment.pos <-
                which(temp.ring.sc.sum$segment.l > segment.size)
              
              #'@special
              #' In this section after,
              #' [01_From "Big-Segment"] is going to search after the
              #' assigned "Big-segment" positions (big.segment.pos).
              #' However, 
              #' It might miss segments that found "closely-grouped" together.
              #'@solution
              #' If updated Possible Big-segments, which will be handled in
              #' next section 01_From "Big-Segment" for Big-Segment search,
              #' the 1st (depends on S.Flip control) Big-segment is within
              #' *morph.id*, this MEANS that the found ".mid.big" is the 
              #' "big-segments" the 01_From "Big-Segment" is looking for.
              #' Hence (*If TRUE*), 
              #' (1) jump Big-Segment-search; and 
              #' (2) handle joint Big-segments;
              #' (3) S.Flip MATTERS.
              #' (4) *middle.III == TRUE*
              #' Otherwise (*If FALSE*),
              #' updated "big.segment.pos" & "temp.ring.sc.sum" will be used 
              #' in 01_From "Big-Segment" section with updated "ring.2".
              
              if(purrr::is_empty(big.segment.pos)!=TRUE){
                
                .ck.edge01 <- temp.ring.sc.sum$clst[big.segment.pos[1]]
                .ck.edge02 <- temp.ring.sc.sum$clst[tail(big.segment.pos, n=1)]
                if((S.Flip != TRUE & .ck.edge01 %in% morph.id) |
                   (S.Flip == TRUE & .ck.edge02 %in% morph.id) ){
                  
                  middle.III <- TRUE
                  message("... Same joint Segments Found")
                  
                }
                
              }else{
                
                # Development Description:
                #'@note
                #' This else ONLY serve for the possible future...
                #' However, normally, if there are no "Big-segments" available,
                #' the small segment search will be initialized by 
                #' [02_Find Small-Segments] + [03_Check "Gap" in-between]
                
                # stop("Check Middle Update with no Big-Segments")
                
              }
              
            }
            
            #'[05_05:Decision02]
            #'@description 
            #' (1) ...ring.keep.sum:
            #'     is the kept segments in ring.keep for re-building reference
            #'     structure (ring.2 -> ring.2 + .mid.ring.2)
            #' (2) ..ring.i.sum:
            #'     is the .mid.big that overlapped with considered segments
            #'     (...ring.keep.sum). If EMPTY, means that .mid.big exist only
            #'     within the gap between Edge.01 & Edge.02.
            #'@conditions
            #' (1) ..ring.i.sum
            #' ..ring.i.sum EMPTY: middle.III
            #'                     middle.I (*Non-Overlap*)
            #'                     - min(..keep.sum.value) < 3 | 
            #'                     - min(.gap.l, .gap.r) < 0.5
            #' ..ring.i.sum Exist: middle.I (*Overlap*)
            #'                     - min(..keep.sum.value) < 3 | 
            #'                     - min(.gap.l, .gap.r) < 0.5
            #'                     
            #' (2) Condition middle.I
            #' min(..keep.sum.value): diffs between ring.keep & .mid.ring2;
            #' min(.gap.l, .gap.r): How close are the .mid.big to both
            #'                      edges of ring.keep.
            #'
            #' (3) Condition middle.III
            #' confirmed above, simply adapted into the code below.
            
            if(purrr::is_empty(..ring.i.sum$clst) != TRUE |
               middle.I == TRUE | middle.III == TRUE){
              
              # Check
              # stop("EMPTY ..ring.i.sum | middle.I | middle.III")
              # ..ring.i.sum$delta %>% hist()
              
              if(purrr::is_empty(..ring.i.sum$clst)){
                
                #'[Scenario01 ..ring.i.sum EMPTY]
                #' .mid.big exist only within gap between Edge.01 & Edge.02.
                #'@note
                #' (1) min(..keep.sum.value): differences 
                #'     between ring.keep & .mid.ring2, 
                #'     *how likely the same are ring.keep & .mid.ring2* 
                #'     tree ring structures.
                #' (2) min(.gap.l, .gap.r): 
                #'     ..ring.i.sum EMPTY -> *Non-Overlap*.
                #'     Hence, if code is not accepted by 1st diff. < 3 cond.,
                #'     The reference ring should be adjusted as middle search 
                #'     results may not belong to current tree ring structure.
                
                if(min(..keep.sum.value) < 3){
                  
                  # Report
                  message("... Kept segments agree with middle-check:")
                  print(...ring.keep.sum)
                  print(..keep.sum.value)
                  
                  # FAKE DATA
                  ..ring.i.sum <- data.frame(clst     = "FAKE",
                                             theta    = 0,
                                             dist     = 0,
                                             dist.ref = 0,
                                             delta    = 0)
                  
                }else if(middle.III == TRUE){
                  
                  # Report
                  # message("... Same joint Segments Found") # Line 2655
                  print(temp.ring.sc.sum[big.segment.pos,])
                  
                  # FAKE DATA
                  ..ring.i.sum <- data.frame(clst     = "FAKE",
                                             theta    = 0,
                                             dist     = 0,
                                             dist.ref = 0,
                                             delta    = 0)
                  
                  # Signal off
                  middle.III <- FALSE
                  
                }else{
                  
                  # Check
                  # stop("Check min(.gap.l, .gap.r) < 0.5 Agreement ...")
                  
                  #'@note
                  #' Build up summary table for adjusting reference tree ring:
                  #' (1) ..ring.i.sum; & 
                  #' (2) corresponding morph section (S.Flip)
                  
                  # Report
                  message("... Kept segments reject middle-check:")
                  print(...ring.keep.sum)
                  print(..keep.sum.value)
                  
                  # Check ..ring.i.sum conditions
                  stop("Check Empty ..ring.i.sum & ..keep.sum.value < 3")
                  plot(test2$theta, # von Line 648
                       test2$dist,
                       col = as.factor(test2$clst),
                       cex = 0.3)
                  points(..ring.i$theta, ..ring.i$dist, cex = 0.3, col = "red")
                  
                  # Fix
                  if(S.Flip == TRUE){
                    
                    ..rb <- min(..keep.sum.list[[1]]$min.theta)
                    ..ring.i.sum <- 
                      ...ring.keep %>% 
                      dplyr::filter(clst %in% ..keep.sum.list[[1]]$clst) %>% 
                      dplyr::filter(theta < (..rb + 0.2)) %>% 
                      dplyr::filter(theta > ..rb) %>%
                      dplyr::mutate(delta = dist - dist.ref)
                    
                  }else{
                    
                    ..lb <- max(..keep.sum.list[[2]]$max.theta)
                    ..ring.i.sum <- 
                      ...ring.keep %>% 
                      dplyr::filter(clst %in% ..keep.sum.list[[2]]$clst) %>% 
                      dplyr::filter(theta > (..lb - 0.2)) %>% 
                      dplyr::filter(theta < ..lb) %>%
                      dplyr::mutate(delta = dist - dist.ref)
                    
                  }
                  
                }
                
              }
              
              if(mean(abs(..ring.i.sum$delta)) > 3){
                
                message("# Reshape ring.Morph at ring.ID ",
                        .mid.big.sum$clst, " position")
                
                # stop("Please develop this section ...")
                
                ..ring.keep.overlap <-
                  ..ring.keep %>% 
                  dplyr::filter(theta < (max(..ring.i.sum$theta)+.M.Step)) %>% 
                  dplyr::filter(theta > (min(..ring.i.sum$theta)-.M.Step))
                
                #'@special
                #' Here, "S.Flip" does not matter,
                #' but the position that overlaps does.
                #' Hence, Use the signal to determine the ring.Morph
                middle.L <- FALSE
                middle.R <- FALSE
                if(any(..ring.i.sum$theta <= Edge.01)){middle.L <- TRUE}
                if(any(..ring.i.sum$theta >= Edge.02)){middle.R <- TRUE}
                
                if(middle.L == TRUE & middle.R != TRUE){
                  
                  ...ring.i.sum <-
                    ..ring.i.sum %>% 
                    dplyr::arrange(desc(theta)) %>% 
                    dplyr::filter(theta < Edge.02) %>% 
                    dplyr::filter(theta > (max(..ring.i.sum$theta)-0.2))
                  .dist.shift <- ...ring.i.sum$delta %>% mean()
                  
                  # stop("check middle.L = TRUE & middle.R != TRUE")
                  
                  # Check
                  # plot(..ring.keep.overlap$theta,
                  # ..ring.keep.overlap$dist)
                  # points(...ring.i.sum$theta,
                  #        ...ring.i.sum$dist - .dist.shift,
                  #        col = "red3")
                  # points(mid.ring.2$theta,
                  #        mid.ring.2$dist + .dist.shift,
                  #        cex = 0.3, col = "lightblue")
                  
                  # Update ring.rm
                  # (1) The Section from .mid.big
                  .ring.rm <- 
                    .mid.big %>% 
                    dplyr::filter(theta > (Edge.01 + .M.Step)) %>% 
                    dplyr::mutate(dist = dist + .dist.shift)
                  
                  # points(.ring.rm$theta,
                  #        .ring.rm$dist, cex = 0.3, col = "red3")
                  # points(.mid.big$theta, 
                  #        .mid.big$dist, cex = 0.3, col = "darkred")
                  # points(mid.ring.2$theta,
                  #        mid.ring.2$dist,
                  #        cex = 0.3, col = "lightblue")
                  
                  
                  # (2) Morph between .mid.big and Edge.02
                  .Edge.01  = max(.ring.rm$theta) %>% round(digits = .digits)
                  .Edge.02  = Edge.02
                  .dist.01  = .ring.rm$dist[which.max(.ring.rm$theta)]
                  .dist.02  = dist.02
                  
                  # (3) Fix
                  #'@special
                  #' For middle.check not overlap but
                  #' disagreed with kept segments.
                  if(Edge.01 < min(.ring.rm$theta)){
                    
                    ..Edge.01  = Edge.01
                    ..Edge.02  = min(.ring.rm$theta) %>% round(digits = .digits)
                    ..dist.01  = dist.01
                    ..dist.02  = .ring.rm$dist[which.min(.ring.rm$theta)]
                    
                    ..ring.rm <- 
                      Section.Morph(Edge.01  = ..Edge.01,
                                    Edge.02  = ..Edge.02,
                                    dist.01  = ..dist.01,
                                    dist.02  = ..dist.02,
                                    ref.ring = ring.2)
                    
                    .ring.rm <- rbind(.ring.rm, ..ring.rm)
                    
                  }
                  
                }
                
                if(middle.L == FALSE & middle.R == TRUE){
                  
                  ...ring.i.sum <-
                    ..ring.i.sum %>% 
                    dplyr::arrange(theta) %>% 
                    dplyr::filter(theta > Edge.01) %>% 
                    dplyr::filter(theta < (min(..ring.i.sum$theta)+0.2))
                  .dist.shift <- ...ring.i.sum$delta %>% mean()
                  
                  # Check
                  # plot(..ring.keep.overlap$theta,
                  #      ..ring.keep.overlap$dist)
                  # points(...ring.i.sum$theta,
                  #        ...ring.i.sum$dist + .dist.shift,
                  #        col = "red3")
                  # points(mid.ring.2$theta,
                  #        mid.ring.2$dist  + .dist.shift,
                  #        cex = 0.3, col = "lightblue")
                  
                  # Update ring.rm
                  # (1) The Section from .mid.big
                  .ring.rm <- 
                    .mid.big %>% 
                    dplyr::filter(theta < (Edge.02 - .M.Step)) %>% 
                    dplyr::mutate(dist = dist + .dist.shift)
                  
                  # points(.ring.rm$theta,
                  #        .ring.rm$dist, cex = 0.3, col = "red3")
                  
                  # (2) Morph between .mid.big and Edge.02
                  .Edge.01  = Edge.01
                  .Edge.02  = min(.ring.rm$theta) %>% round(digits = .digits)
                  .dist.01  = dist.01
                  .dist.02  = .ring.rm$dist[which.min(.ring.rm$theta)]
                  
                  # (3) Fix
                  #'@special
                  #' For middle.check not overlap but
                  #' disagreed with kept segments.
                  if(Edge.02 > max(.ring.rm$theta)){
                    
                    ..Edge.01  = max(.ring.rm$theta) %>% round(digits = .digits)
                    ..Edge.02  = Edge.02 
                    ..dist.01  = .ring.rm$dist[which.max(.ring.rm$theta)]
                    ..dist.02  = dist.02
                    
                    ..ring.rm <- 
                      Section.Morph(Edge.01  = ..Edge.01,
                                    Edge.02  = ..Edge.02,
                                    dist.01  = ..dist.01,
                                    dist.02  = ..dist.02,
                                    ref.ring = ring.2)
                    
                    .ring.rm <- rbind(.ring.rm, ..ring.rm)
                    
                  }
                  
                }
                
                if(middle.L == TRUE & middle.R == TRUE){
                  
                  # Signal PASS Through the NEXT Section
                  # for "ring.rm" Build
                  
                  #'[Method.02]
                  #'@note
                  #' Better than considering from the "large-enough"
                  #' close-by segments from both Edges.
                  #'@usage 
                  #' To use the same function with above,
                  #' Generate "Fake" .Edge.01 & .Edge.02
                  #' In order to pass through the below "#Fix" section
                  
                  # (1) "Fake" .Edge.01 & .Edge.02
                  .Edge.01 <-  Inf
                  .Edge.02 <- -Inf
                  
                  # .ring.rm
                  .ring.rm <- .mid.big
                  .ring.rm$theta <- .ring.rm$theta %>% round(digits = .digits)
                  
                }
                
                # Check
                if(middle.L != TRUE & middle.R != TRUE){
                  stop("Check middle.L & middle.R ... ")
                }
                
                # Fix
                if(.Edge.01 < .Edge.02){
                  
                  .morph.ref <- Section.Morph(Edge.01  = .Edge.01,
                                              Edge.02  = .Edge.02,
                                              dist.01  = .dist.01,
                                              dist.02  = .dist.02,
                                              ref.ring = ring.2)
                  ring.rm <- rbind(.ring.rm, .morph.ref)
                  ring.rm$clst <- "Morph"
                  
                }else{
                  
                  ring.rm <- Section.Morph(Edge.01  = Edge.01,
                                           Edge.02  = Edge.02,
                                           dist.01  = dist.01,
                                           dist.02  = dist.02,
                                           ref.ring = .ring.rm)
                  
                }
                ring.rm$theta <- ring.rm$theta %>% round(digits = .digits)
                
                # Check
                # points(ring.rm$theta,
                #        ring.rm$dist, cex = 0.3, col = "lightblue")
                
                print("Updating temp.ring.sum & big.segment.pos")
                
                temp.section.r <-
                  ring.ck %>% 
                  segment.filter(Edge.01 = Edge.01,
                                 Edge.02 = Edge.02,
                                 ring.segments = .) %>% 
                  ring.join(., ring.segment = ring.rm)
                
                head(temp.section.r)
                
                if(purrr::is_empty(temp.section.r$clst)){
                  
                  # Nothing in between
                  temp.ring.sc.sum <- temp.ring.sc.sum %>% as.data.frame()
                  temp.ring.sc.sum <- temp.section.sum[vector(),]
                  
                }else{
                  
                  # Something in between
                  temp.ring.sc.sum <-
                    temp.section.r %>% 
                    dplyr::mutate(., delta = abs(dist - dist.Morph)) %>% 
                    dplyr::group_by(clst) %>% 
                    dplyr::arrange(delta) %>% 
                    dplyr::summarise(npts      = n(),
                                     ..npts    = round(0.2*npts),
                                     m.delta   = mean(delta[1:..npts]),
                                     sd.delta  = sd  (delta[1:..npts]),
                                     max.theta = max(theta),
                                     min.theta = min(theta),
                                     dist.l    = dist[which.min(theta)],
                                     dist.r    = dist[which.max(theta)]) %>% 
                    dplyr::select(-..npts) %>% 
                    dplyr::filter(., m.delta < 5) %>% 
                    dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                    segment.clean(dist.thres = 2) %>% 
                    segment.overlap(., segment.size = segment.size) %>% 
                    dplyr::arrange(min.theta) %>% 
                    as.data.frame()
                  
                }
                
                # Check
                temp.ring.sc.sum
                
                # Update
                temp.section.sum <- temp.ring.sc.sum
                
                #'@note
                #' No Need to arrange temp.section.sum & temp.ring.sc.sum
                #' The Big-Segment Search is based on ordered summary, and
                #' The Small-Segment Search relies on the fun() itself.
                # temp.ring.sc.sum [Not RUN]
                # if(S.Flip != TRUE){
                #   temp.ring.sc.sum <- 
                #     temp.ring.sc.sum %>% dplyr::arrange(min.theta)
                # }else{
                #   temp.ring.sc.sum <- 
                #     temp.ring.sc.sum %>% dplyr::arrange(desc(min.theta))
                # }
                
                # Update big.segment.pos
                big.segment.pos <-
                  which(temp.ring.sc.sum$segment.l > segment.size)
                
                # Fix
                if(purrr::is_empty(big.segment.pos)!=TRUE){
                  
                  gap.l <- 
                    temp.ring.sc.sum$min.theta[big.segment.pos[1]] - Edge.01
                  gap.r <- Edge.02 - 
                    temp.ring.sc.sum$max.theta[tail(big.segment.pos, n = 1)]
                  if(gap.r < gap.l){S.Flip <- TRUE}else{S.Flip <- FALSE}
                  
                }
                
                # Signal for Final Morph Control
                # middle.II <- TRUE
                message("# [Special] Re-Structure ring.2 ... Middle")
                
                # middle.II <- TRUE
                ..ring.2 <- ring.2
                ring.2 <- 
                  ring.i %>% 
                  dplyr::filter(clst %in% Ring.ID.i) %>% 
                  rbind(., dplyr::bind_rows(morph.list)) %>% 
                  rbind(ring.rm) %>% 
                  segment.GroupRefine()
                ring.2$theta <- ring.2$theta %>% round(digits = .digits)
                ring.2$clst  <- "ref"
                
                if(nrow(ring.2) != nrow(.ring.2)){stop("check ring.2 Update")}
                
                # Just in-case Nothing should leave this section
                # Signal
                Big.Group.Mid <- FALSE
                
              }else if(mean(abs(..ring.i.sum$delta)) <= 3){
                
                message("# Acquire Searched-segments at ring.ID ",
                        .mid.big.sum$clst, " position")
                
                temp.mid.sum <-
                  temp.ring.sc %>% 
                  dplyr::filter(clst %in% morph.id) %>% 
                  # ring.join(., ..ring.i) %>% 
                  # dplyr::filter(is.na(dist.ref)==TRUE) %>% 
                  dplyr::select(names(ring.ck)) %>% 
                  ring.join(., ring.segment = ring.m) %>% 
                  dplyr::mutate(., delta = abs(dist - dist.Morph)) %>% 
                  dplyr::group_by(clst) %>% 
                  dplyr::arrange(delta) %>% 
                  dplyr::summarise(npts      = n(),
                                   ..npts    = round(0.2*npts),
                                   m.delta   = mean(delta[1:..npts]),
                                   sd.delta  = sd  (delta[1:..npts]),
                                   max.theta = max(theta),
                                   min.theta = min(theta),
                                   dist.l    = dist[which.min(theta)],
                                   dist.r    = dist[which.max(theta)]) %>%
                  dplyr::select(-..npts) %>%  
                  dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                  dplyr::arrange(min.theta) %>% 
                  as.data.frame()
                
                # temp.ring.sc.sum
                if(S.Flip != TRUE){
                  temp.ring.sc.sum <- 
                    temp.mid.sum %>% dplyr::arrange(min.theta)
                }else{
                  temp.ring.sc.sum <- 
                    temp.mid.sum %>% dplyr::arrange(desc(min.theta))
                }
                
                #'[Method]
                #'@description 
                #' (1) Use "Big.Group" to wrap all searched segments into
                #'     A single "Big-Segment".
                #' (2) Create temp.ring.sc.sum for further small segment search.
                
                # Signal
                Big.Group.Mid <- TRUE
                
                # Output "Big-Segment"
                temp.sec.big.sum <- .mid.big.sum
                temp.sec.big     <- 
                  ring.ck %>%
                  dplyr::filter(clst %in% temp.ring.sc.sum$clst)
                
                ..mid.big <-
                  .mid.big %>% 
                  dplyr::filter(theta > (min(temp.sec.big$theta) - .M.Step)) %>% 
                  dplyr::filter(theta < (max(temp.sec.big$theta) + .M.Step))
                temp.mid.sum <- temp.ring.sc.sum # With Direction
                
                # Re-build "temp.ring.sc.sum" for small-search
                if(S.Flip != TRUE){
                  
                  .Edge.01  = Edge.01
                  .Edge.02  = min(temp.sec.big$theta)
                  .dist.01  = dist.01
                  .dist.02  = temp.sec.big$dist[which.min(temp.sec.big$theta)]
                  
                }else{
                  
                  .Edge.01  = max(temp.sec.big$theta)
                  .Edge.02  = Edge.02
                  .dist.01  = temp.sec.big$dist[which.max(temp.sec.big$theta)]
                  .dist.02  = dist.02
                  
                }
                
                temp.rm <- Section.Morph(Edge.01  = .Edge.01,
                                         Edge.02  = .Edge.02,
                                         dist.01  = .dist.01,
                                         dist.02  = .dist.02,
                                         ref.ring = ring.2)
                
                temp.section.r <-
                  ring.ck %>% 
                  segment.filter(Edge.01 = .Edge.01,
                                 Edge.02 = .Edge.02,
                                 ring.segments = .) %>% 
                  ring.join(., ring.segment = temp.rm)
                
                
                head(temp.section.r)
                
                if(purrr::is_empty(temp.section.r$clst)){
                  
                  # Nothing in between
                  temp.ring.sc.sum <- data.frame(npts         = integer(),
                                                 m.delta      = numeric(),
                                                 sd.delta     = numeric(),
                                                 max.theta    = numeric(),
                                                 min.theta    = numeric(),
                                                 dist.l       = numeric(),
                                                 dist.r       = numeric(),
                                                 segment.size = numeric())
                  
                }else{
                  
                  # Something in between
                  temp.ring.sc.sum <-
                    temp.section.r %>% 
                    dplyr::mutate(., delta = abs(dist - dist.Morph)) %>% 
                    dplyr::group_by(clst) %>% 
                    dplyr::arrange(delta) %>% 
                    dplyr::summarise(npts      = n(),
                                     ..npts    = round(0.2*npts),
                                     m.delta   = mean(delta[1:..npts]),
                                     sd.delta  = sd  (delta[1:..npts]),
                                     max.theta = max(theta),
                                     min.theta = min(theta),
                                     dist.l    = dist[which.min(theta)],
                                     dist.r    = dist[which.max(theta)]) %>%
                    dplyr::select(-..npts) %>%   
                    dplyr::filter(., m.delta < 5) %>% 
                    dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                    # segment.clean(., dist.thres = 2) %>% 
                    segment.clean(dist.thres = 3) %>% 
                    segment.clean(dist.thres = 2, 
                                  segment.size = segment.size) %>% 
                    segment.overlap(., segment.size = segment.size) %>% 
                    as.data.frame()
                  
                }
                # Check
                temp.ring.sc.sum
                
                # Update
                temp.section.sum <- temp.ring.sc.sum
                
                # Signal
                big.segment.pos <- vector()
                
              }else{
                
                stop("UNKNOWN Condition ... CHECK")
                
              }
              
            }
            
          }else{
            
            # Signal
            small.mid     <- FALSE # Method I
            Big.Group.Mid <- FALSE # Method II
            
          }
          
          
          
          ## (05)_From "Big-Segment" ----
          #' [01_From "Big-Segment" Positions to search "Forward"]
          #' @description
          #' if Big-Segment" does not exist,
          #' the rest will be handled by "Nothing" or "No-Big-Segments"
          
          if(purrr::is_empty(big.segment.pos) != TRUE){
            
            #' [01_02_Forward Search]
            # Save the found Big.Segment Summary
            temp.big.list <- list()
            
            # Fix: Ensure Data Structure
            temp.ring.sc.sum <- 
              temp.ring.sc.sum %>% dplyr::arrange(min.theta)
            
            # use mean() might cause problems 
            # when the segments end slightly before such position
            # In such case, search from "Edge.02" might be a better choice.
            if(Edge.01 < .border.save %>% quantile(0.45) %>% as.numeric() &
               S.Flip == FALSE){ 
              temp.p <- 1 
            }else{
              temp.p <- big.segment.pos %>% length()
            }
            
            temp.big.list[[1]] <- temp.ring.sc.sum[big.segment.pos[temp.p],]
            
            # 01_01_Compute and find small segments in between
            temp.big.i <- 2
            Small.Forw <- FALSE
            temp.section.sum <- temp.ring.sc.sum
            while (TRUE) {
              
              if(Small.Forw == FALSE){
                temp.sec.big.ID <- 
                  temp.section.sum$clst[big.segment.pos[temp.p]]
              }else{
                temp.sec.big.ID <- temp.big.list %>% dplyr::bind_rows()
                temp.sec.big.ID <- temp.sec.big.ID$clst
                temp.sec.big.ID <- 
                  append(temp.sec.big.ID, 
                         temp.section.sum$clst[nrow(temp.section.sum)])
              }
              
              temp.sec.big <- dplyr::filter(ring.ck, clst %in% temp.sec.big.ID)
              
              if(Small.Forw == TRUE){
                temp.big.smsp <- smooth.spline(x  = temp.sec.big$theta,
                                               y  = temp.sec.big$dist,
                                               df = 5)
                temp.sec.big  <-
                  predict(temp.big.smsp,
                          x = temp.sec.big$theta) %>%
                  as.data.frame()
                names(temp.sec.big) <- c("theta", "dist")
                temp.sec.big$theta  <- 
                  temp.sec.big$theta %>% round(digits = .digits)
              }
              
              temp.sec.edge02 <- temp.sec.big$theta %>% min()
              temp.sec.dist02 <- 
                temp.sec.big$dist[which.min(temp.sec.big$theta)]
              temp.rm <- Section.Morph(Edge.01  = Edge.01,
                                       Edge.02  = temp.sec.edge02,
                                       dist.01  = dist.01,
                                       dist.02  = temp.sec.dist02,
                                       ref.ring = ring.2)
              
              temp.section.r <-
                ring.ck %>% 
                segment.filter(Edge.01 = Edge.01,
                               Edge.02 = temp.sec.edge02,
                               ring.segments = .) %>% 
                ring.join(., ring.segment = temp.rm)
              
              head(temp.section.r)
              
              if(purrr::is_empty(temp.section.r$clst)){
                
                # Nothing in between
                temp.section.sum <- temp.ring.sc.sum %>% as.data.frame()
                temp.section.sum <- temp.section.sum[vector(),]
                
              }else{
                
                # Something in between
                temp.section.sum <-
                  temp.section.r %>% 
                  dplyr::mutate(., delta = abs(dist - dist.Morph)) %>% 
                  dplyr::group_by(clst) %>% 
                  dplyr::arrange(delta) %>% 
                  dplyr::summarise(npts      = n(),
                                   ..npts    = round(0.2*npts),
                                   m.delta   = mean(delta[1:..npts]),
                                   sd.delta  = sd  (delta[1:..npts]),
                                   max.theta = max(theta),
                                   min.theta = min(theta),
                                   dist.l    = dist[which.min(theta)],
                                   dist.r    = dist[which.max(theta)]) %>%
                  dplyr::select(-..npts) %>%   
                  dplyr::filter(., m.delta < 5) %>% 
                  dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                  # segment.clean(., dist.thres = 2) %>% 
                  segment.clean(dist.thres = 3) %>% 
                  segment.clean(dist.thres = 2, 
                                segment.size = segment.size) %>% 
                  segment.overlap(., segment.size = segment.size) %>% 
                  dplyr::arrange(min.theta) %>% 
                  as.data.frame()
                
              }
              
              # loop-control
              temp.section.sum
              big.segment.pos <- 
                which(temp.section.sum$segment.l > segment.size)
              
              if(purrr::is_empty(big.segment.pos)){break}
              
              if(Edge.01 < .border.save %>% quantile(0.45) %>% as.numeric() &
                 S.Flip == FALSE){
                
                temp.p <- 1 
                .temp.big.list <- temp.section.sum[big.segment.pos,]
                
              }else{
                
                temp.p <- big.segment.pos %>% length()
                .temp.big.list <- temp.section.sum[big.segment.pos[temp.p],]
                if((temp.sec.edge02 - .temp.big.list$max.theta) > 0.2 &
                   Small.Forw == FALSE){
                  
                  if((temp.sec.edge02 - max(temp.section.sum$max.theta)) < 0.2){
                    Small.Forw <- TRUE
                    .temp.big.list <- temp.section.sum[vector(),]
                  }
                  
                }else{
                  Small.Forw <- FALSE
                }
                
              }
              
              if(purrr::is_empty(.temp.big.list$clst) != TRUE){
                temp.big.list[[temp.big.i]] <- .temp.big.list
                temp.big.i <- temp.big.i + 1
              }
              
            }
            
            
            
            #' [01_03_Group Big-Segments Together - "close-by"]
            # temp.sec.big.sum <- 
            #   temp.big.list %>% 
            #   dplyr::bind_rows() %>% 
            #   dplyr::arrange(min.theta) %>% 
            #   segment.overlap(overlap.check = "All") %>% 
            #   segment.select(neighbor.thres = 0.2)
            
            #' [Add-in: temp Search-flip from edge.01 to edge.02]
            if(S.Flip == TRUE){
              
              temp.sec.big.sum <- 
                temp.big.list %>% 
                dplyr::bind_rows() %>% 
                dplyr::arrange(min.theta) %>% 
                segment.overlap(overlap.check = "All") %>% 
                segment.select(neighbor.thres = 0.2, rev = TRUE)
              
              # Resume the temp.sec.big
              temp.sec.big <- 
                ring.ck %>% dplyr::filter(clst %in% temp.sec.big.sum$clst)
              
              # Rebuild temp.section.sum
              temp.sec.edge01 <- temp.sec.big$theta %>% max()
              temp.sec.dist01 <- 
                temp.sec.big$dist[which.max(temp.sec.big$theta)]
              temp.rm <- Section.Morph(Edge.01  = temp.sec.edge01,
                                       Edge.02  = Edge.02,
                                       dist.01  = temp.sec.dist01,
                                       dist.02  = dist.02,
                                       ref.ring = ring.2)
              
              temp.section.r <-
                ring.ck %>% 
                segment.filter(Edge.01 = temp.sec.edge01,
                               Edge.02 = Edge.02,
                               ring.segments = .) %>% 
                ring.join(., ring.segment = temp.rm)
              
              head(temp.section.r)
              
              if(purrr::is_empty(temp.section.r$clst)){
                
                # Nothing in between
                temp.section.sum <- temp.ring.sc.sum %>% as.data.frame()
                temp.section.sum <- temp.section.sum[vector(),]
                
              }else{
                
                # Something in between
                temp.section.sum <-
                  temp.section.r %>% 
                  dplyr::mutate(., delta = abs(dist - dist.Morph)) %>% 
                  dplyr::group_by(clst) %>% 
                  dplyr::arrange(delta) %>% 
                  dplyr::summarise(npts      = n(),
                                   ..npts    = round(0.2*npts),
                                   m.delta   = mean(delta[1:..npts]),
                                   sd.delta  = sd  (delta[1:..npts]),
                                   max.theta = max(theta),
                                   min.theta = min(theta),
                                   dist.l    = dist[which.min(theta)],
                                   dist.r    = dist[which.max(theta)]) %>%
                  dplyr::select(-..npts) %>%   
                  dplyr::filter(., m.delta < 5) %>% 
                  dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                  # segment.clean(., dist.thres = 2) %>% 
                  segment.clean(dist.thres = 3) %>% 
                  segment.clean(dist.thres = 2, 
                                segment.size = segment.size) %>% 
                  segment.overlap(., segment.size = segment.size) %>% 
                  as.data.frame()
                
              }
              
            }else{
              
              temp.sec.big.sum <- 
                temp.big.list %>% 
                dplyr::bind_rows() %>% 
                dplyr::arrange(min.theta) %>% 
                segment.overlap(overlap.check = "All") %>% 
                segment.select(neighbor.thres = 0.2)
              
            }
            
            # PASS to "Big-segment Summary Check"
            
          }else{
            
            if(small.mid     != TRUE &
               Big.Group.Mid != TRUE){
              
              temp.sec.big     <- ring.ck[vector(),]
              temp.sec.big.sum <- temp.ring.sc.sum[vector(),]
              temp.section.sum <- temp.ring.sc.sum
              
            }
            
          }
          
          # Big-segment Summary Check
          if(purrr::is_empty(temp.sec.big.sum$clst) != TRUE){
            
            # Linkage Check_dist
            if(nrow(temp.sec.big.sum) != 1){
              
              .temp.big <-
                temp.ring.sc %>% 
                dplyr::filter(clst %in% temp.sec.big.sum$clst)
              
              .temp.big.sum <-
                .temp.big %>% 
                dplyr::group_by(clst) %>% 
                dplyr::summarise(t.min = min(theta),
                                 t.max = max(theta),
                                 d.lb  = dist[which.min(theta)],
                                 d.rb  = dist[which.max(theta)],
                                 # dx    = t.max - t.min,
                                 # dy    = d.rb - d.lb,
                                 slope = (d.rb - d.lb)/(t.max - t.min)) %>% 
                dplyr::arrange(t.min)
              
              #' [Add-in] : S.Flip
              #' Arrangement influence ALOT especially when "S.Flip" exist.
              if(S.Flip == TRUE){
                .temp.big.sum <- .temp.big.sum %>% dplyr::arrange(desc(t.min))
              }
              
              .pos <- which(c(0, !!diff(sign(.temp.big.sum$slope))) != 0)
              if(purrr::is_empty(.pos) != TRUE){
                
                # There's Signal Change in Big-Segment.Group
                #'@note 
                #' If the distance between this change position is too much,
                #' BETTER to morph-search separately.
                
                .pos <- .pos[1] # ONLY Take the 1st one
                
                #' [Add-in] : S.Flip
                #' Arrangement influence ALOT especially when "S.Flip" exist.
                if(S.Flip == TRUE){
                  .gap.delta <- 
                    .temp.big.sum$d.lb[(.pos-1)] - .temp.big.sum$d.rb[.pos]
                }else{
                  .gap.delta <- 
                    .temp.big.sum$d.lb[.pos] - .temp.big.sum$d.rb[(.pos-1)]
                }
                
                
                if(abs(.gap.delta) > 2){
                  
                  # The Change of distance between this 2 segments is too big
                  temp.sec.big.sum <- temp.sec.big.sum[c(1:(.pos-1)),]
                  
                  # Update "ALL" Big-Segment related Output
                  temp.sec.big <- 
                    ring.ck %>% 
                    dplyr::filter(clst %in% temp.sec.big.sum$clst)
                  temp.sec.big.sum
                  temp.section.sum
                  
                }
                
              }
              
            }
            
            if(nrow(temp.sec.big.sum) != 1){
              
              # Signal
              Big.Group <- TRUE
              
              .temp.big <-
                temp.ring.sc %>% 
                dplyr::filter(clst %in% temp.sec.big.sum$clst) 
              
              .temp.big.sum <-
                .temp.big %>% 
                dplyr::group_by(clst) %>% 
                dplyr::summarise(t.min = min(theta),
                                 t.max = max(theta),
                                 d.lb  = dist[which.min(theta)],
                                 d.rb  = dist[which.max(theta)]) %>% 
                dplyr::arrange(t.min)
              
              group.i     <- 1
              group.morph <- list()
              while (TRUE) {
                
                group.morph[[group.i]] <- 
                  Section.Morph(Edge.01 = .temp.big.sum$t.max[group.i],
                                Edge.02 = .temp.big.sum$t.min[group.i+1],
                                dist.01 = .temp.big.sum$d.rb [group.i],
                                dist.02 = .temp.big.sum$d.lb [group.i+1],
                                ref.ring = ring.2)
                
                # group.i-loop control
                group.i <- group.i + 1
                if((group.i+1) > nrow(temp.sec.big.sum)){break}
                
              }
              temp.group.morph      <- group.morph %>% dplyr::bind_rows()
              temp.group.morph$clst <- "Morph.B"
              
              # Plot Check
              if(.plotcheck == TRUE){
                
                points(temp.group.morph$theta,
                       temp.group.morph$dist,
                       cex = 0.3,
                       col = "pink")
                
              }
              
            }
            
          }
          
          
          
          ## (06)_Find Small-Segments ----
          #' [02_Find Small-Segments]
          # Fine tune "temp.section.sum" as beginning
          #'@seealso : 
          #'@description:
          #' Revise "temp.ring.sc.sum" & temp.small
          #' This is to prevent really tiny small segments 
          #' influence the searching.
          #' MOST IMPORTANTLY,
          #' this prevent the spline function is applied in a
          #' extremely wide theta range.
          #' Because now small segments are selected based on the 
          #' "1st" Big Segment
          #'@note:
          #' The small segments excluded in this phase should be picked again.
          
          #'@note:
          #' If no big-segment and close-by found, 
          #' the "temp.section.sum" will NOT be reset well for below.
          #' Hence, re-establish "temp.section.sum" is crucial
          if(Small.Back == FALSE & 
             purrr::is_empty(temp.sec.big.sum$clst)){
            temp.section.sum <- temp.ring.sc.sum
          }
          
          #' [02_00_Preparation]
          #' Forward Search Only make sense if 
          #' (1) there are segments close-by
          #' (2) Big-Segments Exists (search from edge.02 is done at beginning)
          #'@note:
          #  This Adjustment is needed only when small.back = F
          #  (No small segments close by Edge.01 & Edge.02)
          if(Small.Back == FALSE &
             purrr::is_empty(temp.sec.big.sum$clst) != TRUE){
            
            # Save
            .temp.section.sum <- temp.section.sum
            
            # Forward Small-Segments grouping 
            # (based on found 1st big segment position)
            if(S.Flip == TRUE){
              
              .ringi_closeby <- 
                temp.section.sum %>% 
                segment.select(neighbor.thres = 0.05, rev = FALSE)
              
            }else{
              
              .ringi_closeby <- 
                temp.section.sum %>% 
                segment.select(neighbor.thres = 0.05, rev = TRUE)
              
            }
            
            # Only when close-by exist, then make sense
            #'@note:
            #' .ringi_closeby == EMPTY acts as a signal 
            #' to trigger the forward-search of small segments
            if(purrr::is_empty(.ringi_closeby$clst) != TRUE){
              
              # In the last search, 
              # it will encounter the tail of the selected segment.
              #'@note:
              #' Hence, in case of "no big-segment",
              #  This is searched already at beginning.
              if(S.Flip == TRUE){
                
                # Reverse the searching from temp.big
                .temp.tail <- max(temp.sec.big$theta)
                
                #  If not close-by, no need to search...
                if((min(.ringi_closeby$min.theta) - .temp.tail) > 0.1){
                  .ringi_closeby <- .ringi_closeby[vector(),]
                }
                
              }else{
                
                .temp.tail <- min(temp.sec.big$theta)
                
                #  If not close-by, no need to search...
                if((.temp.tail - max(.ringi_closeby$max.theta)) > 0.1){
                  .ringi_closeby <- .ringi_closeby[vector(),]
                }
                
              }
            }
            
          }else{
            
            .ringi_closeby <- temp.section.sum[vector(),]
            
          }
          
          if(Small.Back == TRUE){ # Small.Back-control
            
            #' [02_01_Small-Segments Grouped close to Edge.01]
            # Already searched.
            
            # temp.small.sum  <- temp.section.sum
            # temp.small.clst <- temp.small.sum$clst
            # temp.small      <- temp.ring.sc %>% 
            #   dplyr::filter(., clst %in% temp.small.sum$clst)
            # Small.Back <- FALSE
            
          }else if(purrr::is_empty(.ringi_closeby$clst) != TRUE){
            
            #' [02_02_Small-Segments_Forward-Search]
            #  If no grouped segments, no need to search...
            
            # Get Big segments base for search
            # The search will NOT happen if there is NO big-segments.
            #'@note:
            #' TWO scenarios:
            #' (1) S.Flip == TRUE,  closing Edge.02;
            #' (1) S.Flip == FALSE, closing Edge.01;
            if(S.Flip == TRUE){
              
              temp.ring.big <- 
                temp.sec.big %>% 
                dplyr::filter(theta > (max(temp.sec.big$theta) - 0.5))
              
              if(nrow(temp.ring.big) > 10*ang.npts){ # 10 means 10 degrees
                
                temp.ring.big <- temp.ring.big[order(-temp.ring.big$theta),]
                temp.ring.big <- temp.ring.big[c(1:10*ang.npts),]
                
              }
              
              # Control
              rev.control <- FALSE
              
            }else{
              
              temp.ring.big <- 
                temp.sec.big %>% 
                dplyr::filter(theta < (min(temp.sec.big$theta) + 0.5))
              
              if(nrow(temp.ring.big) > 10*ang.npts){ # 10 means 10 degrees
                
                temp.ring.big <- temp.ring.big[order(temp.ring.big$theta),]
                temp.ring.big <- temp.ring.big[c(1:10*ang.npts),]
                
              }
              
              # Control
              rev.control <- TRUE
              
            }
            
            
            # Save the found Small.Segment Summary
            temp.small.list      <- list()
            temp.small.list[[1]] <- .ringi_closeby
            
            # Compute and find small segments in between
            temp.small.i <- 2
            
            while (TRUE) {
              
              temp.section.sum <- temp.small.list %>% dplyr::bind_rows()
              temp.small.ID    <- temp.section.sum$clst
              temp.sec.small   <- 
                ring.ck %>% 
                dplyr::filter(clst %in% temp.small.ID) %>% 
                rbind(temp.ring.big)
              
              # spline degree of freedom decision
              .temp.range <- 
                max(temp.section.sum$max.theta) - 
                min(temp.section.sum$min.theta)
              .temp.npts  <- temp.section.sum$npts %>% sum()
              if(.temp.range > 1.0 | .temp.npts > 1000){
                .temp.df <- 10}else{.temp.df <- 5}
              
              temp.small.smsp  <- smooth.spline(x  = temp.sec.small$theta,
                                                y  = temp.sec.small$dist,
                                                df = .temp.df)
              temp.sp <-
                predict(temp.small.smsp,
                        x = temp.sec.small$theta) %>%
                as.data.frame()
              names(temp.sp) <- c("theta", "dist")
              temp.sp$theta  <- temp.sp$theta %>% round(digits = 4)
              
              # Check
              # plot(temp.sec.small$theta, temp.sec.small$dist)
              # lines(temp.small.smsp, col = "darkgreen")
              
              #' [Add-in]:S.Flip
              if(S.Flip == TRUE){
                
                temp.sec.edge01 <- temp.sp$theta %>% max()
                temp.sec.dist01 <- temp.sp$dist[which.max(temp.sp$theta)]
                temp.sec.edge02 <- Edge.02
                temp.sec.dist02 <- dist.02
                temp.rm <- Section.Morph(Edge.01  = temp.sec.edge01,
                                         Edge.02  = temp.sec.edge02,
                                         dist.01  = temp.sec.dist01,
                                         dist.02  = temp.sec.dist02,
                                         ref.ring = ring.2)
                
              }else{
                
                temp.sec.edge01 <- Edge.01
                temp.sec.dist01 <- dist.01
                temp.sec.edge02 <- temp.sp$theta %>% min()
                temp.sec.dist02 <- temp.sp$dist[which.min(temp.sp$theta)]
                temp.rm <- Section.Morph(Edge.01  = Edge.01,
                                         Edge.02  = temp.sec.edge02,
                                         dist.01  = dist.01,
                                         dist.02  = temp.sec.dist02,
                                         ref.ring = ring.2)
                
              }
              
              
              temp.section.tail <-
                ring.ck %>% 
                segment.filter(Edge.01 = temp.sec.edge01,
                               Edge.02 = temp.sec.edge02,
                               ring.segments = .) %>% 
                ring.join(., ring.segment = temp.rm)
              
              
              head(temp.section.tail)
              
              if(purrr::is_empty(temp.section.tail$clst)){
                
                # Nothing in between
                temp.section.sum <- temp.ring.sc.sum %>% as.data.frame()
                temp.section.sum <- temp.section.sum[vector(),]
                
              }else{
                
                # Something in between
                temp.section.sum <-
                  temp.section.tail %>% 
                  dplyr::mutate(., delta = abs(dist - dist.Morph)) %>% 
                  dplyr::group_by(clst) %>% 
                  dplyr::arrange(delta) %>% 
                  dplyr::summarise(npts      = n(),
                                   ..npts    = round(npts*0.1),
                                   m.delta   = mean(delta[1:..npts]),
                                   sd.delta  = sd  (delta[1:..npts]),
                                   # m.Ang.pos = npts / nrow(ring.2),
                                   max.theta = max(theta),
                                   min.theta = min(theta),
                                   dist.l    = dist[which.min(theta)],
                                   dist.r    = dist[which.max(theta)]) %>%
                  dplyr::select(-..npts) %>% 
                  dplyr::filter(., m.delta < 5) %>% 
                  dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                  # segment.clean(., dist.thres = 2) %>% 
                  segment.clean(dist.thres = 3) %>% 
                  segment.clean(dist.thres = 2, 
                                segment.size = segment.size) %>% 
                  segment.overlap(., segment.size = segment.size) %>% 
                  as.data.frame()
                
              }
              
              # loop-control
              temp.section.sum
              .ringi_closeby <- 
                segment.select(ring.summary   = temp.section.sum,
                               neighbor.thres = 0.05,
                               rev = rev.control)
              
              # Only when there are something close-by, then make sense
              if(purrr::is_empty(.ringi_closeby$clst) != TRUE){
                
                # Only when segments are close enough, then make sense
                .range <- 
                  max(.ringi_closeby$max.theta) - 
                  min(.ringi_closeby$min.theta)
                # (1) Scenario-01
                if(S.Flip == TRUE){
                  
                  if(((.ringi_closeby$min.theta[1] - temp.sec.edge01) > 0.1 &
                      sum(.ringi_closeby$segment.l) < segment.size) |
                     (.range < 0.1 &
                      nrow(.ringi_closeby) < 3 &
                      sum(.ringi_closeby$segment.l) < 0.8*segment.size)){
                    
                    .ringi_closeby <- .ringi_closeby[vector(),]
                    
                  }
                  
                }
                # (2) Scenario-02
                if(S.Flip != TRUE){
                  
                  if(((temp.sec.edge02 - .ringi_closeby$max.theta[1]) > 0.1 &
                      sum(.ringi_closeby$segment.l) < segment.size) | 
                     (.range < 0.1 &
                      nrow(.ringi_closeby) < 3 &
                      sum(.ringi_closeby$segment.l) < 0.8*segment.size)){
                    
                    .ringi_closeby <- .ringi_closeby[vector(),]
                    
                  }
                  
                }
                
              }
              
              if(purrr::is_empty(.ringi_closeby$clst)){break}
              
              temp.small.list[[temp.small.i]] <- .ringi_closeby
              temp.small.i <- temp.small.i + 1
              
            }
            
            # temp.small.sum
            temp.small.sum <- 
              temp.small.list %>% 
              dplyr::bind_rows() %>% 
              dplyr::arrange(min.theta)
            
            # Only when segments are big enough, then make sense
            # 1000 = 0.1 theta range
            # old: temp.small.sum$npts %>% sum() < 1000
            .temp.range <- 
              max(temp.small.sum$max.theta) - min(temp.small.sum$min.theta)
            if(.temp.range < 0.1 |
               temp.small.sum$segment.l %>% sum() < segment.size){
              temp.small.sum <- temp.small.sum[vector(),]
            }
            
            # Small-segments
            # temp.small.clst <- temp.small.sum$clst
            # temp.small <-
            #   temp.ring.sc %>%
            #   dplyr::filter(., clst %in% temp.small.sum$clst)
            .small.tail <- temp.small.sum
            
            # Fix: Resume the temp.section.sum
            if(purrr::is_empty(.small.tail$clst)){
              
              # report
              print("# Small Search from Tail ... Failed")
              temp.section.sum <- .temp.section.sum
              
            }
            
            
          }else{
            
            # Nothing should pass through
            .small.tail <- temp.section.sum[vector(),]
            
          }
          
          
          # Special condition
          #'@note:
          #' When the derived "temp.section.sum" have other clusters left,
          #' We need to check whether there are other "grouped" segments.
          if(Small.Back == TRUE){ # Small.Back-control
            
            #' [02_01_Small-Segments Grouped close to Edge.01]
            # Already searched.
            
            # temp.small.sum  <- temp.section.sum
            # temp.small.clst <- temp.small.sum$clst
            # temp.small      <- temp.ring.sc %>% 
            #   dplyr::filter(., clst %in% temp.small.sum$clst)
            # Small.Back <- FALSE
            
          }else if(purrr::is_empty(temp.section.sum$clst) != TRUE){
            
            #' [02_03_Small-Segments_Backward-Search]
            #  This is to capture any small-segment-groups
            #  in-between Edge.01 and min(big.segments$theta) or temp.small
            
            # Save the previously searched results:
            .small.middle <- list()
            .ms <- 1
            # segment.select control
            if(S.Flip == TRUE){rev.control <- TRUE}else{rev.control <- FALSE}
            while (TRUE) {
              
              # loop-control - I
              if(purrr::is_empty(temp.section.sum$clst) == TRUE){
                
                if(.ms == 1){
                  .small.middle[[.ms]] <- temp.section.sum
                }
                break
                
              }
              
              temp.small.sum <- 
                temp.section.sum %>% 
                segment.select(neighbor.thres = 0.05, rev = rev.control)
              
              # Only when segments are big enough, then make sense
              # 1000 = 0.1 theta range
              # old: temp.small.sum$npts %>% sum() >= 1000
              .temp.range <- 
                max(temp.small.sum$max.theta) - min(temp.small.sum$min.theta)
              if(sum(temp.small.sum$segment.l) >= segment.size |
                 .temp.range >= 0.1){
                
                #'@note
                #' For specific purpose.
                .small.middle[[.ms]] <- temp.small.sum
                .ms <- .ms + 1
                
              }
              
              # loop-control
              temp.sec.ID <- dplyr::setdiff(temp.section.sum$clst, 
                                            temp.small.sum$clst)
              temp.section.sum <-
                temp.section.sum %>% 
                dplyr::filter(., clst %in% temp.sec.ID)
              
            }
            
          }else{
            
            .small.middle <- list()
            .small.middle[[1]] <- temp.section.sum[vector(),]
            
          }
          
          # Summarize the searched small segments
          if(small.mid == TRUE){
            
            # Signal
            small.mid <- FALSE
            
            # Small Segments already found ...
            temp.small.sum
            temp.small.clst
            temp.small
            
          }else if(Small.Back == TRUE){ # Small.Back-control
            
            #' [02_01_Small-Segments Grouped close to Edge.01]
            # Already searched.
            
            temp.small.sum  <- temp.section.sum
            temp.small.clst <- temp.small.sum$clst
            temp.small      <- temp.ring.sc %>% 
              dplyr::filter(., clst %in% temp.small.sum$clst)
            Small.Back <- FALSE
            
          }else{
            
            # Targets:
            .small.tail
            .small.middle
            
            # Preparation (Remove the empty list)
            # Fix
            if(S.Flip == TRUE){
              .temp.small <- append(list(.small.tail), .small.middle)
            }else{
              .temp.small <- append(.small.middle, list(.small.tail))
            }
            
            .ms <- 1
            while (TRUE) {
              
              if(purrr::is_empty(.temp.small[[.ms]]$clst)){
                .temp.small[[.ms]] <- NULL
                # loop-control
                .ms <- 1
              }else{
                # loop-control
                .ms <- .ms + 1
              }
              
              if(.ms > length(.temp.small)){break}
              
            }
            
            #' [Add-in] : S.Flip
            #'@special
            #' when big-segment is missing, the argument will be handled
            #' in "No Big Segments EXIST" section.
            #' However, when
            #' (1) .small.tail is found 
            #' (2) Big-segment is missing
            #' if there are anything left before .small.middle,
            #' this will be ignored by the search...
            #' To avoid this,
            #' this section offer another way to signal "S.Flip"
            #' as without "big-segment",
            #' such search from Edge.02 (.small.tail) just acts as "S.Flip".
            #'@details 
            #' (1) .gap.r < .gap.l & .gap.l > 1
            #' (2) To release more morph-search from .small.middle, 
            #'     distance between each group should < 0.15
            if(purrr::is_empty(.temp.small) != TRUE &
               purrr::is_empty(temp.sec.big.sum$clst)){
              
              if(S.Flip == TRUE){
                
                .pos <- length(.temp.small)
                .gap.l <- min(.temp.small[[.pos]]$min.theta) - Edge.01
                .gap.r <- Edge.02 - max(.temp.small[[1]]$max.theta)
                
              }else{
                
                .pos <- length(.temp.small)
                .gap.l <- min(.temp.small[[1]]$min.theta) - Edge.01
                .gap.r <- Edge.02 - max(.temp.small[[.pos]]$max.theta)
                
              }
              
              S.Flip.Adj <- FALSE
              if(.gap.r < .gap.l){
                
                if(S.Flip != TRUE){
                  # Report
                  message("# Search-Flip from Tail #")
                  S.Flip.Adj <- TRUE
                }
                S.Flip <- TRUE
                
                # Check
                print("# Section III")
                
              }else{
                
                if(S.Flip == TRUE){
                  # Report
                  message("# Search-Flip BACK ... #")
                  S.Flip.Adj <- TRUE
                }
                S.Flip <- FALSE
                
              }
              
              if(S.Flip.Adj == TRUE){
                
                # Fix the order of ".temp.small"
                ..temp.small <- .temp.small %>% dplyr::bind_rows()
                if(S.Flip == TRUE){rev.control <- TRUE
                }else{rev.control <- FALSE}
                
                .ms = 1
                ...temp.small <- list()
                while (TRUE) {
                  
                  # loop-control - I
                  if(purrr::is_empty(..temp.small$clst) == TRUE){
                    
                    if(.ms == 1){
                      ...temp.small[[.ms]] <- ..temp.small
                    }
                    break
                    
                  }
                  
                  temp.small.sum <- 
                    ..temp.small %>% 
                    segment.select(neighbor.thres = 0.05, rev = rev.control)
                  
                  ...temp.small[[.ms]] <- temp.small.sum
                  .ms <- .ms + 1
                  
                  # loop-control
                  temp.sec.ID <- dplyr::setdiff(..temp.small$clst, 
                                                temp.small.sum$clst)
                  ..temp.small <-
                    ..temp.small %>% 
                    dplyr::filter(., clst %in% temp.sec.ID)
                  
                }
                .temp.small <- ...temp.small
                
              }
              
              if(.pos > 1){
                
                # (1) Calculate the distance in between groups
                # (2) To decide which side to select -> No need
                #     as "S.Flip" activated, always get from right.
                
                if(S.Flip == TRUE){
                  
                  .s.save  <- .temp.small
                  .s.dist  <- vector()
                  .s.theta <- vector()
                  .sd <- 1
                  while (TRUE) {
                    
                    .temp.theta <- 
                      min(.temp.small[[ .sd ]]$min.theta) - 
                      max(.temp.small[[.sd+1]]$max.theta)
                    
                    .temp.dist <- 
                      min(.temp.small[[ .sd ]]$dist.l) - 
                      max(.temp.small[[.sd+1]]$dist.r)
                    
                    # loop control
                    if(.temp.theta < 0.15 &
                       .temp.dist  < 5){
                      
                      .s.theta <- append(.s.theta, .temp.theta)
                      .s.dist  <- append( .s.dist, .temp.dist)
                      .sd <- .sd - 1
                      
                    }else{
                      
                      .delete <- .pos
                      while (.delete >= (.sd+1)) {
                        .temp.small[[.delete]] <- NULL
                        # loop-control
                        .delete <- .delete - 1
                      }
                      break
                      
                    }
                    
                    # loop-control II
                    if(.sd <= 1){break}
                    
                  }
                  
                }else{
                  
                  .s.save  <- .temp.small
                  .s.dist  <- vector()
                  .s.theta <- vector()
                  .sd <- 1
                  while (TRUE) {
                    
                    .temp.theta <- 
                      min(.temp.small[[.sd+1]]$min.theta) -
                      max(.temp.small[[ .sd ]]$max.theta)
                    
                    .temp.dist <-
                      min(.temp.small[[.sd+1]]$dist.l) -
                      max(.temp.small[[ .sd ]]$dist.r)
                    
                    # loop control
                    if(.temp.theta < 0.15 &
                       .temp.dist  < 5){
                      
                      .s.theta <- append(.s.theta, .temp.theta)
                      .s.dist  <- append( .s.dist, .temp.dist)
                      .sd <- .sd + 1
                      
                    }else{
                      
                      .delete <- .pos
                      while (.delete >= (.sd+1)) {
                        .temp.small[[.delete]] <- NULL
                        # loop-control
                        .delete <- .delete - 1
                      }
                      break
                      
                    }
                    
                    # loop-control II
                    if(.sd == .pos){break}
                    
                  }
                  
                }
                
              }
              
            }
            
            
            if(purrr::is_empty(.temp.small)){
              
              # Nothing Should be passed out from this trunk
              temp.small.sum  <- .small.tail[vector(),]
              temp.small.clst <- temp.small.sum$clst
              temp.small      <- temp.ring.sc[vector(),]
              
            }else if(length(.temp.small) > 1){
              
              # More than 1 group...
              # To Pass through the morphing section,
              # It is better to merge these small groups with 
              # "insertion-supports"
              
              # Signal (No need)
              # The added segments will be filtered out.
              # Small.sup <- TRUE
              
              # Fix the order of ".temp.small"
              ..temp.small <- .temp.small %>% dplyr::bind_rows()
              if(S.Flip == TRUE){rev.control <- TRUE
              }else{rev.control <- FALSE}
              
              .ms = 1
              ...temp.small <- list()
              while (TRUE) {
                
                # loop-control - I
                if(purrr::is_empty(..temp.small$clst) == TRUE){
                  
                  if(.ms == 1){
                    ...temp.small[[.ms]] <- ..temp.small
                  }
                  break
                  
                }
                
                temp.small.sum <- 
                  ..temp.small %>% 
                  segment.select(neighbor.thres = 0.05, rev = rev.control)
                
                ...temp.small[[.ms]] <- temp.small.sum
                .ms <- .ms + 1
                
                # loop-control
                temp.sec.ID <- dplyr::setdiff(..temp.small$clst, 
                                              temp.small.sum$clst)
                ..temp.small <-
                  ..temp.small %>% 
                  dplyr::filter(., clst %in% temp.sec.ID)
                
              }
              .temp.small <- ...temp.small
              
              # How much are there?
              n.ss <- length(.temp.small) - 1
              
              # Make
              small.sup <- list()
              .ss <- 1
              while (TRUE) {
                
                if(S.Flip != TRUE){
                  
                  lb.clst <- .temp.small[[.ss]]    $clst
                  lb.clst <- lb.clst[length(lb.clst)]
                  rb.clst <- .temp.small[[(.ss+1)]]$clst[1]
                  
                }else{
                  
                  rb.clst <- .temp.small[[.ss]]    $clst
                  rb.clst <- rb.clst[length(rb.clst)]
                  lb.clst <- .temp.small[[(.ss+1)]]$clst[1]
                  
                }
                
                
                # Morph the "insertion-supports" section
                lb <- temp.ring.sc %>% dplyr::filter(clst %in% lb.clst)
                rb <- temp.ring.sc %>% dplyr::filter(clst %in% rb.clst)
                .ss.seg <- 
                  Section.Morph(Edge.01 = lb$theta %>% max(),
                                Edge.02 = rb$theta %>% min(),
                                dist.01 = lb$dist[which.max(lb$theta)],
                                dist.02 = rb$dist[which.min(rb$theta)],
                                ref.ring = ring.2)
                .ss.seg$clst <- "Small.sup"
                
                # loop-control
                small.sup[[.ss]] <- .ss.seg
                .ss <- .ss + 1
                if(.ss == length(.temp.small)){break}
                
              }
              small.sup <- small.sup %>% dplyr::bind_rows()
              
              # Output
              temp.small.sum  <- 
                .temp.small %>% dplyr::bind_rows()
              temp.small.clst <- temp.small.sum$clst
              temp.small      <- 
                temp.ring.sc %>% 
                dplyr::filter(clst %in% temp.small.clst) %>% 
                dplyr::select(c("clst", "theta", "dist")) %>% 
                rbind(., small.sup)
              
            }else{
              
              # Only 1 small segment group is found
              # Pass safely.
              temp.small.sum  <- 
                .temp.small %>% dplyr::bind_rows()
              temp.small.clst <- temp.small.sum$clst
              temp.small      <- 
                temp.ring.sc %>% 
                dplyr::filter(clst %in% temp.small.clst)
              
            }
            
          }
          
          
          
          ## (07)_Check "Gap" in-between ----
          #' [03_Check "Gap" in-between Big-Segment and found "small-Segments"]
          #' @description 
          #  As the small segments searching in "Big-Segment Exist"
          #  Only handles one groups of small segments,
          #  If the "Gap" is too big,
          #  There are opportunities to find other segments.
          if(purrr::is_empty(temp.sec.big.sum$clst) != TRUE &
             purrr::is_empty(temp.small.sum$clst) != TRUE){
            
            # Boundary from Big-Segments
            #' [Add-in]:S.Flip
            if(S.Flip == TRUE){
              
              # Correspond Boundary from Grouped Big Segments
              temp.big.rb <- temp.sec.big.sum$max.theta %>% max()
              
              # Correspond Boundary from Grouped small Segments
              temp.small.lb <- temp.small.sum$min.theta %>% min()
              
              # Gap_in_Between
              Gap_in_Between <- temp.small.lb - temp.big.rb
              
            }else{
              
              # Correspond Boundary from Grouped Big Segments
              temp.big.lb <- temp.sec.big.sum$min.theta %>% min()
              
              # Correspond Boundary from Grouped small Segments
              temp.small.rb <- temp.small.sum$max.theta %>% max()
              
              # Gap_in_Between
              Gap_in_Between <- temp.big.lb - temp.small.rb
              
            }
                       
            if(Gap_in_Between > 1){
              
              # Gap is to big
              # -> Clean Big-Segment Position to Pass the code to
              # -> "No-Big-Segment Exist"
              big.segment.pos  <- vector()
              temp.sec.big     <- ring.ck[vector(),]
              temp.sec.big.sum <- temp.section.sum[vector(),]
              temp.section.sum 
              
              # Report
              message("Huge Gap in between grouped small and big segments...keep small")
              
              # Big.Group Control
              if(Big.Group == TRUE |  Big.Group.Mid == TRUE){
                
                Big.Group <- FALSE
                Big.Group.Mid <- FALSE
                
              }
            }
            
          }
          
          
          
          ## (08)_Last search of small segments ----
          #'@special
          #' When "temp.sec.big.sum" & "temp.small.sum" are EMPTY
          #' However, in the "temp.ring.sc.sum" is not EMPTY,
          #' Any segments in such region with small m.delta will be taken.
          #' With this fact in mind,
          #' This search is to perform additional search (small segments)
          #' before the k-loop is closed by "[Nothing Remains]"
          if(purrr::is_empty(temp.sec.big.sum$clst) == TRUE &
             purrr::is_empty(temp.small.sum  $clst) == TRUE &
             purrr::is_empty(temp.ring.sc.sum$clst) != TRUE){
            
            message("# Final small-segment search #")
            print(paste("Between Edge.01 (", Edge.01,
                        ") & Edge.02 (", Edge.02,
                        "); dist = ", 
                        round((Edge.02 - Edge.01),digits = .digits), sep = ""))
            
            # Check
            if((Edge.02 - Edge.01) < 0.5){
              
              # Preparation
              if(purrr::is_empty(morph.list)!=TRUE){
                temp.final <- morph.list %>% dplyr::bind_rows()
              }else{
                temp.final <- ring.i
              }
              
              # Segments from Edge.01
              temp.final.l <- temp.final[which(temp.final$theta < Edge.02),]
              if(nrow(temp.final.l) > 10*ang.npts){ # 10 means 10 degrees
                temp.final.l <- temp.final.l[order(-temp.final.l$theta),]
                temp.final.l <- temp.final.l[c(1:(10*ang.npts)),]
              }
              # Segments from Edge.02
              temp.final.r <- temp.final[which(temp.final$theta > Edge.01),]
              if(nrow(temp.final.r) > 10*ang.npts){ # 10 means 10 degrees
                temp.final.r <- temp.final.r[order( temp.final.r$theta),]
                temp.final.r <- temp.final.r[c(1:(10*ang.npts)),]
              }
              # Segments from temp.ring.sc.sum
              temp.final.s <- 
                temp.ring.sc %>% 
                dplyr::filter(clst %in% temp.ring.sc.sum$clst) %>% 
                dplyr::select(c("clst", "theta", "dist"))
              
              # temp.smdf
              temp.smdf <- rbind(temp.final.l, temp.final.r)
              temp.smdf <- rbind(temp.smdf, temp.final.s)
              
              temp.sm.k  <- smooth.spline(x  = temp.smdf$theta,
                                          y  = temp.smdf$dist,
                                          df = 10)
              temp.rm <-
                predict(temp.sm.k, 
                        x = seq(Edge.01, Edge.02, by = ReSample.Size)) %>%
                as.data.frame()
              names(temp.rm) <- c("theta", "dist")
              temp.rm$clst   <- "Morph"
              temp.rm$theta  <- temp.rm$theta %>% round(digits = .digits)
              
              temp.section.r <-
                ring.ck %>% 
                segment.filter(Edge.01 = Edge.01,
                               Edge.02 = Edge.02,
                               ring.segments = .) %>% 
                ring.join(., ring.segment = temp.rm)
              
              head(temp.section.r)
              
              if(purrr::is_empty(temp.section.r$clst)){
                
                # Nothing in between
                temp.ring.sc.sum <- temp.ring.sc.sum %>% as.data.frame()
                temp.ring.sc.sum <- temp.section.sum[vector(),]
                
              }else{
                
                # Something in between
                temp.ring.sc.sum <-
                  temp.section.r %>% 
                  dplyr::mutate(., delta = abs(dist - dist.Morph)) %>% 
                  dplyr::group_by(clst) %>% 
                  dplyr::arrange(desc(theta)) %>% 
                  dplyr::summarise(npts      = n(),
                                   ..npts    = round(0.2*npts),
                                   m.delta   = mean(delta[1:..npts]),
                                   sd.delta  = sd  (delta[1:..npts]),
                                   max.theta = max(theta),
                                   min.theta = min(theta),
                                   dist.l    = dist[which.min(theta)],
                                   dist.r    = dist[which.max(theta)]) %>%
                  dplyr::select(-..npts) %>%
                  dplyr::filter(., m.delta < 5) %>% 
                  dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                  # segment.clean(., dist.thres = 2) %>% 
                  segment.clean(dist.thres = 3) %>% 
                  segment.clean(dist.thres = 2, 
                                segment.size = segment.size) %>% 
                  segment.overlap(., segment.size = segment.size) %>% 
                  as.data.frame()
                
              }
              
              # Check
              temp.ring.sc.sum
              
              # Plot Check
              # test <- rbind(temp.final.l, temp.final.r)
              # plot(test$theta, test$dist)
              # points(temp.final.s$theta, 
              #        temp.final.s$dist, col = "lightblue")
              # lines(temp.sm.k, col = "darkgreen")
              
              #' [Fix:Tail-check]
              temp.ring.sc.sum <- 
                temp.ring.sc.sum %>% dplyr::arrange(min.theta)
              if(purrr::is_empty(temp.ring.sc.sum$clst) != TRUE &
                 nrow(temp.ring.sc.sum) != 1){
                
                # Left-Check
                .gap.l  <- temp.ring.sc.sum$min.theta[1] - Edge.01
                .dist.l <- abs(temp.ring.sc.sum$dist.l[1] - dist.01)
                if(.gap.l < 0.1 & .dist.l > 5){
                  temp.ring.sc.sum <- temp.ring.sc.sum[-1,]
                }
                
                # Right-Check
                .gap.r  <- Edge.02 - tail(temp.ring.sc.sum$max.theta, n = 1)
                .dist.r <- abs(dist.02 - tail(temp.ring.sc.sum$dist.r, n = 1))
                if(.gap.r < 0.1 & .dist.r > 5){
                  temp.ring.sc.sum <- temp.ring.sc.sum[-nrow(temp.ring.sc.sum),]
                }
                
              }
              
              if(purrr::is_empty(temp.ring.sc.sum$clst) != TRUE &
                 nrow(temp.ring.sc.sum) != 1){
                
                # Group Small-segments
                temp.sum.e1 <- 
                  temp.ring.sc.sum %>% 
                  segment.select(neighbor.thres = 0.05)
                .gap.l <- min(temp.sum.e1$min.theta) - Edge.01
                
                temp.sum.e2 <- 
                  temp.ring.sc.sum %>% 
                  segment.select(neighbor.thres = 0.05, rev = TRUE)
                .gap.r <- Edge.02 - max(temp.sum.e2$max.theta)
                
                # Output
                temp.small.sum  <- temp.sum.e1
                temp.small.clst <- temp.small.sum$clst
                temp.small      <- 
                  temp.ring.sc %>% 
                  dplyr::filter(clst %in% temp.small.clst)
                
                if(.gap.r < .gap.l){
                  
                  # Output
                  temp.small.sum  <- temp.sum.e2
                  temp.small.clst <- temp.small.sum$clst
                  temp.small      <- 
                    temp.ring.sc %>% 
                    dplyr::filter(clst %in% temp.small.clst) %>% 
                    dplyr::select(c("clst", "theta", "dist"))
                  
                  # Report
                  message("# Search-Flip from Tail #")
                  S.Flip <- TRUE
                  
                  # Check
                  print("# Section IV")
                  
                }
                
                
              }else{
                
                message("# Final Search Ends ")
                
                # Ensure nothing pass through
                temp.small.sum  <- temp.ring.sc.sum[vector(),]
                temp.small.clst <- temp.small.sum$clst
                temp.small      <- 
                  ring.ck %>% 
                  dplyr::filter(clst %in% temp.small.clst) %>% 
                  dplyr::select(c("clst", "theta", "dist"))
                
              }
              
            }else{
              
              # stop("Gap Check")
              
              #'@description 
              #' When the gap is large, the final check is to 
              #' use the closest small segments to the edge.01/ edge.02,
              #' Mimic a segment to parts of edge, and perform a 
              #' "Middle Check"
              #'@note 
              #' The section block "Middle Check" "Tail Check"
              #' Could have been built as function in the future.
              #'@details:
              #' To decide which side of small segment is used to 
              #' fake as a "big-middle", 
              #' the previously searched segments by "close-by" 
              #' should be removed.
              
              ..e1 <- segment.select(ring.summary   = temp.ring.sc.sum,
                                     neighbor.thres = 0.05,
                                     dist.check = TRUE,
                                     dist.thres = 5)
              if((min(..e1$min.theta) - Edge.01)<0.1){
                
                .pos <- c(1:nrow(..e1))
                ..e1 <- 
                  segment.select(ring.summary   = temp.ring.sc.sum[-.pos,],
                                 neighbor.thres = 0.05,
                                 dist.check = TRUE,
                                 dist.thres = 5)
                
              }
              if(purrr::is_empty(..e1$clst)){
                .gap.l <- Inf
              }else{
                .gap.l <- min(..e1$min.theta) - Edge.01
              }
              ..e2 <- segment.select(ring.summary   = temp.ring.sc.sum,
                                     neighbor.thres = 0.05,
                                     rev = TRUE,
                                     dist.check = TRUE,
                                     dist.thres = 5)
              if((Edge.02 - max(..e2$max.theta))<0.1){
                
                .pos <- seq(nrow(temp.ring.sc.sum), 
                            by = -1, length.out = nrow(..e2))
                ..e2 <- 
                  segment.select(ring.summary   = temp.ring.sc.sum[-.pos,],
                                 neighbor.thres = 0.05,
                                 rev = TRUE,
                                 dist.check = TRUE,
                                 dist.thres = 5)
                
              }
              if(purrr::is_empty(..e2$clst)){
                .gap.r <- Inf
              }else{
                .gap.r <- Edge.02 - max(..e2$max.theta)
              }
              
              
              # Selection
              #'@note
              #' When ..e1 & ..e2 are Empty,
              #' It will not matter whether .mid.big.sum = ..e1 or ..e2
              #' Additionally, 
              #' always let the search from edge.01 goes first.
              if(.gap.l <= .gap.r){
                .mid.big.sum <- ..e1[1,]
              }
              if(.gap.r < .gap.l){
                .mid.big.sum <- ..e2[1,]
                
                # Report
                message("# Search-Flip from Tail #")
                S.Flip <- TRUE
                print("# Section IV")
                
              }
              
              # Search within the final gap
              # (Middle Check Structure)
              if(purrr::is_empty(.mid.big.sum$clst) != TRUE){
                
                # While control ..e1 and ..e2
                while (TRUE) {
                  
                  # While EMPTY-control
                  if(purrr::is_empty(..e1$clst) &
                     purrr::is_empty(..e2$clst)){
                    
                    if(S.Flip == TRUE){S.Flip <- FALSE}
                    message("# Final Atempt Failed ... Nothing left")
                    break
                    
                  }
                  
                  # Report
                  # message("# Middle Check ...")
                  
                  # Target .mid.big assigned... ABOVE
                  .mid.big.sum 
                  
                  ####DEV----
                  #'[Forced:Small-Middle.Search]
                  mid.big.id   <- .mid.big.sum$clst
                  mid.big.df   <- ring.ck %>% filter(clst %in% mid.big.id)
                  npts.range   <- 0.2
                  middle.check <- 
                    segment.MiddleCheck(target.segment.df   = mid.big.df,
                                        considered.segments = ring.ck,
                                        ref.ring            = ring.2,
                                        npts.range          = npts.range,
                                        Force.small         = TRUE,
                                        segment.SUM         = lookup.SUM,
                                        segment.size        = segment.size,
                                        ReSample.Size       = ReSample.Size)
                  
                  # Use Results directly from "middle.check"
                  .temp.section.sum <- middle.check$temp.section.sum
                  .mid.big.sum      <- middle.check$mid.big.sum
                  morph.id          <- middle.check$morph.id
                  .mid.big          <- middle.check$mid.big
                  mid.ring.2        <- middle.check$mid.ring.2
                  
                  # Check
                  if(.plotcheck == TRUE){
                    
                    points(.mid.big$theta, 
                           .mid.big$dist, cex = 0.3, col = "darkred")
                    
                  }
                  
                  # stop("Check Section 04")
                  
                  #'[GapDist_Check:Forced-Search close to Both Edges]
                  #' Middle-Searched Reference
                  #'@note Use mid.ring.2 instead of .mid.big to 
                  .temp.mid.ring2 <- 
                    mid.ring.2 %>% 
                    dplyr::filter(theta > (Edge.01 + .M.Step)) %>% 
                    dplyr::filter(theta < (Edge.02 - .M.Step)) %>% 
                    dplyr::arrange(theta)
                  
                  # Boundaries & dists close to both edges
                  ..mid.lb <- .temp.mid.ring2$theta %>% min()
                  ..mid.rb <- .temp.mid.ring2$theta %>% max()
                  ..mid.ld <- .temp.mid.ring2$dist  %>% head(1)
                  ..mid.rd <- .temp.mid.ring2$dist  %>% tail(1)
                  
                  .d.theta <- abs(c((Edge.01 - ..mid.lb),(Edge.02 - ..mid.rb)))
                  .d.dist  <- abs(c((dist.01 - ..mid.ld),(dist.02 - ..mid.rd)))
                  if(any(.d.theta < 0.05 & .d.dist > 5)){
                    
                    # Report
                    message("# Final Search with Tail Adjustment ...")
                    
                    # Check
                    # points(mid.ring.2$theta,
                    #        mid.ring.2$dist,
                    #        cex = 0.3, col = "lightblue")
                    
                    # .ring.rm Adjustments
                    .mid.big.lb <- .mid.big$theta %>% min()
                    .mid.big.rb <- .mid.big$theta %>% max()
                    .adj.ln  <- .temp.mid.ring2$theta %>% length()
                    .adj.dl  <- dist.01 - ..mid.ld
                    .adj.dr  <- dist.02 - ..mid.rd
                    dist.adj <- seq(.adj.dl, .adj.dr, length.out = .adj.ln)
                    .ring.rm <-
                      .temp.mid.ring2 %>% 
                      dplyr::mutate(dist = dist + dist.adj) %>% 
                      dplyr::filter(theta > (.mid.big.lb - .M.Step)) %>% 
                      dplyr::filter(theta < (.mid.big.rb + .M.Step)) %>% 
                      dplyr::arrange(theta)
                    
                    # Morph-searched & ring.i
                    .ring.morph <- 
                      morph.list %>% 
                      dplyr::bind_rows() %>% segment.GroupRefine()
                    .ring.i  <- ring.ck %>% dplyr::filter(clst %in% Ring.ID.i)
                    
                    # In-between
                    .ring.M.L <- 
                      Section.Morph(Edge.01 = Edge.01,
                                    Edge.02 = .ring.rm$theta %>% head(1),
                                    dist.01 = dist.01,
                                    dist.02 = .ring.rm$dist  %>% head(1),
                                    ref.ring = ring.2)
                    .ring.M.R <- 
                      Section.Morph(Edge.01 = .ring.rm$theta %>% tail(1),
                                    Edge.02 = Edge.02,
                                    dist.01 = .ring.rm$dist  %>% tail(1),
                                    dist.02 = dist.02,
                                    ref.ring = ring.2)
                    .ring.M <- rbind(.ring.M.L, .ring.M.R)
                    
                    # Compilation
                    ring.rm <- 
                      .ring.rm %>% 
                      rbind(., .ring.morph) %>% 
                      rbind(., .ring.i) %>% 
                      rbind(., .ring.M)
                    ring.rm$theta <- ring.rm$theta %>% round(digits = .digits)
                    ring.rm$clst  <- "ref"
                    
                    # Check
                    if(nrow(ring.rm) != nrow(ring.2)){
                      stop("Check 04_Last Search ... ring.rm")}
                    
                    # Check
                    # points(.ring.morph$theta, 
                    #        .ring.morph$dist, cex = 0.5, col = "lightblue")
                    # points(.ring.rm$theta, 
                    #        .ring.rm$dist, cex = 0.5, col = "lightblue")
                    # points(.ring.i$theta, 
                    #        .ring.i$dist, cex = 0.5, col = "darkred")
                    # points(.ring.M$theta, 
                    #        .ring.M$dist, cex = 0.5, col = "darkred")
                    
                    .sum <-
                      temp.ring.sc %>% 
                      ring.join(., ring.segment = ring.rm)%>% 
                      dplyr::group_by(clst) %>% 
                      dplyr::mutate(., delta = abs(dist - dist.ref)) %>% 
                      dplyr::summarise(npts      = n(),
                                       ..npts    = round(0.2*npts),
                                       m.delta   = mean(delta[1:..npts]),
                                       sd.delta  = sd  (delta[1:..npts]),
                                       max.theta = max(theta),
                                       min.theta = min(theta),
                                       dist.l    = dist[which.min(theta)],
                                       dist.r    = dist[which.max(theta)]) %>%
                      dplyr::select(-..npts) %>%  
                      dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                      dplyr::filter(m.delta < 2) %>% 
                      dplyr::arrange(min.theta) %>% 
                      segment.overlap(., overlap.check = "All") %>% 
                      as.data.frame()
                    
                    # S.Flip Control
                    if(purrr::is_empty(.sum$clst) != TRUE){
                      
                      .gap.l <- min(.sum$min.theta) - Edge.01
                      .gap.r <- Edge.02 - max(.sum$max.theta)
                      
                      if(.gap.r < .gap.l){
                        
                        if(S.Flip == FALSE){
                          # Report
                          message("# Search-Flip from Tail #")
                          S.Flip <- TRUE
                          
                          # Check
                          print("# Section IV")
                        }
                        
                        # rev.control
                        rev.control <- TRUE
                        
                      }else{
                        
                        if(S.Flip == TRUE){
                          # Report
                          message("# Search-Flip BACK ... #")
                          S.Flip <- FALSE
                        }
                        
                        # rev.control
                        rev.control <- FALSE
                        
                      }
                      
                    }
                    
                    if(S.Flip == TRUE){
                      .sum <- .sum[order(-.sum$min.theta),]
                    }else{
                      .sum <- .sum[order( .sum$min.theta),]
                    }
                    
                    # Out-put Control
                    if(any(.sum$segment.l > segment.size)){
                      
                      .pos <- which(.sum$segment.l > segment.size)[1]
                      
                      # Big-Segment exist
                      temp.sec.big     <- 
                        ring.ck %>% dplyr::filter(clst %in% .sum$clst[.pos])
                      temp.sec.big.sum <- .sum[.pos,]
                      
                      # Small Segment exist
                      temp.small.sum  <- 
                        .sum[c(1:(.pos-1)),] %>% 
                        segment.select(neighbor.thres = 0.05, 
                                       rev = rev.control)
                      temp.small.clst <- temp.small.sum$clst
                      temp.small      <- 
                        ring.ck %>% dplyr::filter(clst %in% temp.small.clst)
                      
                    }else{
                      
                      # Big-Segment No-existance
                      temp.sec.big     <- ring.ck[vector(),]
                      temp.sec.big.sum <- .sum[vector(),]
                      
                      # Small Segment exist
                      temp.small.sum  <- 
                        .sum %>% 
                        segment.select(neighbor.thres = 0.05, 
                                       rev = rev.control)
                      temp.small.clst <- temp.small.sum$clst
                      temp.small      <- 
                        ring.ck %>% dplyr::filter(clst %in% temp.small.clst)
                      
                    }
                    
                    # Update ring.2
                    message("# [Special] Re-Structure ring.2 ...Last Search")
                    ring.2 <- ring.rm
                    
                    # Loop-control
                    break
                    
                  }
                  
                  #'[Get_Segments:based on Middle-Check]
                  #'@note
                  #' The section ONLY works when there is no huge distance gap
                  #' in between .mid.big & morph-searched results.
                  #'@seealso 
                  #' Section: "GapDist_Check"
                  
                  # Segments found by middle-search
                  morph.id
                  
                  # Code Option
                  # (1) Check when there is any big-segment in morph.ID
                  # (2) Segments in morph.id might not exist in *sum data frame
                  
                  # Structure the cluster-summary
                  .sum <-
                    temp.ring.sc %>% 
                    dplyr::filter(clst %in% morph.id) %>% 
                    dplyr::select(names(ring.2)) %>% 
                    ring.join(., ring.segment = ring.2)%>% 
                    dplyr::group_by(clst) %>% 
                    dplyr::mutate(., delta = abs(dist - dist.ref)) %>% 
                    dplyr::summarise(npts      = n(),
                                     ..npts    = round(0.2*npts),
                                     m.delta   = mean(delta[1:..npts]),
                                     sd.delta  = sd  (delta[1:..npts]),
                                     max.theta = max(theta),
                                     min.theta = min(theta),
                                     dist.l    = dist[which.min(theta)],
                                     dist.r    = dist[which.max(theta)]) %>%
                    dplyr::select(-..npts) %>%
                    dplyr::left_join(., lookup.SUM, by = "clst") %>% 
                    segment.overlap(segment.size = segment.size,
                                    overlap.check = "All")%>% 
                    as.data.frame()
                  
                  #'@special
                  #' this limits the search for too small segments 
                  #' dominates this search procedure.
                  if(all(.sum$segment.l < segment.size) & (nrow(.sum) == 1)){
                    .sum <- .sum[vector(),]
                  }
                  
                  # While control
                  if(purrr::is_empty(.sum$clst) != TRUE){
                    
                    #' Fix: Adjust .sum Based on how close it is to Edges
                    .gap.l <- min(.sum$min.theta) - Edge.01
                    .gap.r <- Edge.02 - max(.sum$max.theta)
                    if(.gap.r < .gap.l){S.Flip <- TRUE}else{FALSE}
                    
                    if(S.Flip == TRUE){
                      .sum <- .sum[order(-.sum$min.theta),]
                    }else{
                      .sum <- .sum[order( .sum$min.theta),]
                    }
                    
                    if(any(.sum$segment.l > segment.size)){
                      
                      .pos <- which(.sum$segment.l > segment.size)[1]
                      
                      # Big-Segment exist
                      temp.sec.big     <- 
                        ring.ck %>% dplyr::filter(clst %in% .sum$clst[.pos])
                      temp.sec.big.sum <- .sum[.pos,]
                      
                      # Small Segment exist
                      # Fix
                      if(.pos != 1){
                        temp.small.sum <- .sum[c(1:(.pos-1)),]
                      }else{
                        temp.small.sum <- .sum[vector(),]
                      }
                      temp.small.clst <- temp.small.sum$clst
                      temp.small      <- 
                        ring.ck %>% dplyr::filter(clst %in% temp.small.clst)
                      
                    }else{
                      
                      # Small Segment exist
                      temp.small.sum  <- .sum
                      temp.small.clst <- temp.small.sum$clst
                      temp.small      <- 
                        ring.ck %>% dplyr::filter(clst %in% temp.small.clst)
                      
                    }
                    
                    message("# Segments acquired...")
                    break
                    
                  }
                  
                  # Switch-Control
                  if(purrr::is_empty(..e1$clst) != TRUE){
                    
                    if(.mid.big.sum$clst %in% ..e1$clst){
                      
                      .mid.big.sum <- ..e2[1,]
                      ..e1 <- ..e1[vector(),]
                      
                      # Report
                      message("# Search-Flip from Tail #")
                      S.Flip <- TRUE
                      
                      # Check
                      print("# Section III")
                      
                    }else{ #.mid.big.sum$clst %in% ..e2$clst
                      
                      .mid.big.sum <- ..e1[1,]
                      ..e2 <- ..e2[vector(),]
                      
                      # Report
                      message("# Search-Flip Failed... loop back #")
                      S.Flip <- FALSE
                      
                    }
                    
                  }else if(purrr::is_empty(..e2$clst) != TRUE){
                    
                    if(.mid.big.sum$clst %in% ..e2$clst){
                      
                      .mid.big.sum <- ..e1[1,]
                      ..e2 <- ..e2[vector(),]
                      
                      # Report
                      message("# Search-Flip Failed... Search Stop #")
                      S.Flip <- FALSE
                      
                    }
                    
                  }
                  
                }
                # While control ..e1 and ..e2
                
              }else{
                
                # Check the Situation
                # If it is ok, then PASS
                message("# No segments found with in large GAP : END")
                
              }
              
            }
            
          }
          
          # ------------------------------------------------------------------ #
          # 05_temp.ring.sc.sum Update?
          #' @note 
          #' No need to update "temp.ring.sc.sum"
          #' (1) [Nothing Remains] -> No Need
          #' (2) [No Big Segments EXIST] -> No Need (Only Temp.small*)
          #' (3) [Big Segments EXIST] -> Need Original One
          # ------------------------------------------------------------------ #
          
          
          
          # 2.2_Morph Segments ----
          ## (01)_Nothing Remains ----
          #' [Nothing Remains]
          if(purrr::is_empty(temp.sec.big.sum$clst) == TRUE &
             purrr::is_empty(temp.small.sum$clst)   == TRUE){
            
            ring.keep[[k]] <- data.frame(clst  = character(),
                                         theta = numeric(),
                                         dist  = numeric())
            break
            
          }else{
            
            # Signal for Final Morph Control
            #'@note
            #' "middle.II" is ONLY used when nothing is found
            #' after the specific "Middle-Check" section.
            #' If this is the case, 
            #' the rest gap will be filled by "ring.rm" from 
            #' such section.
            middle.II <- FALSE
            
          }
          
          
          
          ## (02)_No Big Segments EXIST ----
          #' [No Big Segments EXIST]
          if(purrr::is_empty(temp.sec.big.sum$clst) == TRUE){
            
            # Preparation
            #' [Add-in]: Better Morphing section
            #' If the small segment morphing does not consider the close edges,
            #' The connection in-between will not be smooth.
            #' reference ring: ring.i
            if(S.Flip == FALSE){
              
              temp.smdf <- rbind(ring.i, morph.list %>% dplyr::bind_rows())
              temp.smdf <- temp.smdf[which(temp.smdf$theta < Edge.02),]
              if(nrow(temp.smdf) > 10*ang.npts){ # 10 means 10 degrees
                temp.smdf <- temp.smdf[order(-temp.smdf$theta),]
                temp.smdf <- temp.smdf[c(1:(10*ang.npts)),]
              }
              
            }else{
              
              temp.smdf <- ring.i[which(ring.i$theta > Edge.01),]
              if(nrow(temp.smdf) > 10*ang.npts){ # 10 means 10 degrees
                temp.smdf <- temp.smdf[order( temp.smdf$theta),]
                temp.smdf <- temp.smdf[c(1:(10*ang.npts)),,]
              }
              
            }
            .temp.small <- rbind(temp.small[names(temp.smdf)], temp.smdf)
            
            #' [Step01_Morph within small ring segments]
            if(purrr::is_empty(temp.small.clst) == TRUE){
              
              # Normally This should NOT Happen, the loop should already break
              message("[Warning]: The loop still continues")
              temp.ring.m <- data.frame(clst  = character(),
                                        theta = numeric(),
                                        dist  = numeric())
              
            }else if(length(temp.small.clst) != 1){
              
              temp.theta.r <- temp.small$theta %>% range() %>% dist()
              temp.dist.r  <- temp.small$dist  %>% range() %>% dist()
              if(temp.theta.r > 1){
                
                temp.smsp <- smooth.spline(x  = .temp.small$theta,
                                           y  = .temp.small$dist,
                                           df = 30)
                
              }else if(temp.dist.r  > 8 | temp.theta.r > 0.5){
                
                temp.smsp <- smooth.spline(x  = .temp.small$theta,
                                           y  = .temp.small$dist,
                                           df = 10)
                
              }else if(temp.theta.r < (5/180*pi)){
                
                temp.smsp <- smooth.spline(x  = .temp.small$theta,
                                           y  = .temp.small$dist,
                                           df = 3)
                
              }else{
                
                temp.smsp <- smooth.spline(x  = .temp.small$theta,
                                           y  = .temp.small$dist,
                                           df = 4)
                
              }
              
              temp.x      <- seq(min(temp.small$theta),
                                 max(temp.small$theta),
                                 by = 0.0001)
              temp.sm     <- predict(temp.smsp, x = temp.x) %>% as.data.frame()
              temp.ring.m <- data.frame(clst  = "clst.M",
                                        theta = temp.sm$x,
                                        dist  = temp.sm$y)
              
              # Plot Check
              if(.plotcheck == TRUE){
                
                points(temp.ring.m$theta, 
                       temp.ring.m$dist, 
                       col = "pink", cex = 0.3)
                
              }
              
              
            }else{ # ONLY ONE small segment
              
              temp.ring.m <- data.frame(clst  = temp.small$clst,
                                        theta = temp.small$theta,
                                        dist  = temp.small$dist)
              
            }
            
            # Save the selected ring segments
            ring.keep[[k]] <- 
              temp.small %>% 
              dplyr::mutate(clst = clst %>% as.character())%>% 
              dplyr::select(c("clst", "theta", "dist")) %>% 
              dplyr::filter(clst != "Small.sup")
            
            
            
            #' [Step02_Morph within large ring segments]
            #' [ONLY Left Space]
            if(purrr::is_empty(temp.ring.m$clst) != TRUE){
              
              #' [Left of small ring segment group]
              
              #' [Add-in] : S.Flip
              if(S.Flip == FALSE){
                
                #'@note
                # clst.L is the target ring segments 
                # clst.R is temp.ring.m
                
                temp.edge01 <- Edge.01
                temp.edge02 <- temp.ring.m$theta %>% min()
                
                temp.dist01 <- dist.01
                temp.dist02 <- temp.ring.m$dist[which.min(temp.ring.m$theta)]
                
              }else{
                
                # S.Flip
                
                temp.edge01 <- temp.ring.m$theta %>% max()
                temp.edge02 <- Edge.02
                
                temp.dist01 <- temp.ring.m$dist[which.max(temp.ring.m$theta)]
                temp.dist02 <- dist.02
                
              }
              
              temp.sec.ls <- Section.Morph(Edge.01  = temp.edge01,
                                           Edge.02  = temp.edge02,
                                           dist.01  = temp.dist01,
                                           dist.02  = temp.dist02,
                                           ref.ring = ring.2)
              
              # Plot Check
              if(.plotcheck == TRUE){
                
                points(temp.sec.ls$theta, 
                       temp.sec.ls$dist, 
                       col = "red3", cex = 0.3)
                
              }
              
              
              temp.sec.k <- temp.sec.ls
              
            }else{
              
              # Normally This should NOT Happen, the loop should already break
              message("[Warning]: The loop still continues")
              temp.sec.k <- data.frame(clst  = character(),
                                       theta = numeric(),
                                       dist  = numeric())
              
            }
            
            # Continue with:
            #' [Step03_Compile & Update temp.ring.sc to process the loop]
            
          }
          
          
          
          ## (03)_Big Segments EXIST ----
          #' [Big Segments EXIST]
          if(purrr::is_empty(temp.sec.big.sum$clst) != TRUE){
            
            #' [Step01_Morph within small ring segments]
            if(purrr::is_empty(temp.small.sum$clst) != TRUE){
              
              # Preparation
              temp.small.clst
              temp.small
              
              # stop("BIG temp.check")
              
              #' [Add-in]: Better Morphing section
              #' If the small segment morphing does not consider close edges,
              #' The connection in-between will not be smooth.
              #' reference ring: ring.i
              if(S.Flip == FALSE){
                
                temp.smdf.l <- rbind(ring.i, morph.list %>% dplyr::bind_rows())
                temp.smdf.l <- temp.smdf.l[which(temp.smdf.l$theta < Edge.02),]
                if(nrow(temp.smdf.l) > 5*ang.npts){ # 10 means 10 degrees
                  temp.smdf.l <- temp.smdf.l[order(-temp.smdf.l$theta),]
                  temp.smdf.l <- temp.smdf.l[c(1:(5*ang.npts)),]
                }
                
                temp.smdf.r <- temp.sec.big
                if(nrow(temp.smdf.r) > 5*ang.npts){ # 10 means 10 degrees
                  temp.smdf.r <- temp.smdf.r[order( temp.smdf.r$theta),]
                  temp.smdf.r <- temp.smdf.r[c(1:(5*ang.npts)),]
                }
                
              }else{
                
                temp.smdf.l <- temp.sec.big
                if(nrow(temp.smdf.l) > 5*ang.npts){ # 10 means 10 degrees
                  temp.smdf.l <- temp.smdf.l[order(-temp.smdf.l$theta),]
                  temp.smdf.l <- temp.smdf.l[c(1:(5*ang.npts)),]
                }
                
                temp.smdf.r <- ring.i[which(ring.i$theta > Edge.01),]
                if(nrow(temp.smdf.r) > 5*ang.npts){ # 10 means 10 degrees
                  temp.smdf.r <- temp.smdf.r[order( temp.smdf.r$theta),]
                  temp.smdf.r <- temp.smdf.r[c(1:(5*ang.npts)),]
                }
                
              }
              temp.smdf   <- rbind(temp.smdf.l, temp.smdf.r)
              .temp.small <- rbind(temp.small[names(temp.smdf)], temp.smdf)
              
              if(length(temp.small.clst) != 1){
                
                temp.theta.r <- temp.small$theta %>% range() %>% dist()
                temp.dist.r  <- temp.small$dist  %>% range() %>% dist()
                if(temp.theta.r > 1){
                  
                  temp.smsp <- smooth.spline(x  = .temp.small$theta,
                                             y  = .temp.small$dist,
                                             df = 30)
                  
                }else if(temp.dist.r  > 8 | temp.theta.r > 0.5){
                  
                  temp.smsp <- smooth.spline(x  = .temp.small$theta,
                                             y  = .temp.small$dist,
                                             df = 10)
                  
                }else if(temp.theta.r < (5/180*pi)){
                  
                  temp.smsp <- smooth.spline(x  = .temp.small$theta,
                                             y  = .temp.small$dist,
                                             df = 3)
                  
                }else{
                  
                  temp.smsp <- smooth.spline(x  = .temp.small$theta,
                                             y  = .temp.small$dist,
                                             df = 4)
                  
                }
                
                temp.x    <- seq(min(temp.small$theta),
                                 max(temp.small$theta),
                                 by = ReSample.Size)
                temp.sm   <- predict(temp.smsp, x = temp.x) %>% as.data.frame()
                temp.ring.m <- data.frame(clst  = "clst.M",
                                          theta = temp.sm$x,
                                          dist  = temp.sm$y)
                
                # Plot Check
                if(.plotcheck == TRUE){
                  
                  points(x   = temp.ring.m$theta, 
                         y   = temp.ring.m$dist, 
                         col = "pink", cex = 0.3)
                  
                }
                
                
              }else if(temp.small.sum$segment.l > segment.size){
                
                temp.ring.m <- temp.small[,c("clst", "theta", "dist")]
                
              }else{
                
                temp.ring.m <- data.frame(clst  = temp.small$clst,
                                          theta = temp.small$theta,
                                          dist  = temp.small$dist)
                
              }
              
            }else{
              
              temp.small  <- temp.ring.sc[vector(),]
              temp.ring.m <- data.frame(clst  = character(),
                                        theta = numeric(),
                                        dist  = numeric())
              
            }
            
            temp.small <- 
              temp.small %>% 
              dplyr::mutate(clst = clst %>% as.character())%>% 
              dplyr::select(c("clst", "theta", "dist")) %>% 
              dplyr::filter(clst != "Small.sup")
            
            
            
            #' [Step02_Morph within large ring segments]
            #' [ONLY Left Space]
            
            #' [Preparation]
            temp.sec.big
            
            #' [Case01_small ring segment group EXISTS]
            if(purrr::is_empty(temp.ring.m$clst) != TRUE){
              
              #' [Left of small ring segment group]
              
              #' [Add-in] : S.Flip
              if(S.Flip == FALSE){
                
                #'@note
                # clst.L is the target ring segments 
                # clst.R is temp.ring.m
                
                temp.edge01 <- Edge.01
                temp.edge02 <- temp.ring.m$theta %>% min()
                
                temp.dist01 <- dist.01
                temp.dist02 <- temp.ring.m$dist[which.min(temp.ring.m$theta)]
                
              }else{
                
                # S.Flip
                
                temp.edge01 <- temp.sec.big$theta %>% max()
                temp.edge02 <- temp.ring.m $theta %>% min()
                
                temp.dist01 <- temp.sec.big$dist[which.max(temp.sec.big$theta)]
                temp.dist02 <- temp.ring.m $dist[which.min(temp.ring.m $theta)]
                
              }
              
              temp.sec.ls <- Section.Morph(Edge.01  = temp.edge01,
                                           Edge.02  = temp.edge02,
                                           dist.01  = temp.dist01,
                                           dist.02  = temp.dist02,
                                           ref.ring = ring.2)
              
              # Plot Check
              if(.plotcheck == TRUE){
                
                points(x   = temp.sec.ls$theta, 
                       y   = temp.sec.ls$dist, 
                       col = "red3", cex = 0.3)
                
              }
              
              #' [Right of small ring segment group]
              
              #' [Add-in] : S.Flip
              if(S.Flip == FALSE){
                
                #'@note 
                # clst.L is temp.ring.m
                # clst.R is temp.sec.big at big.segment.pos[k]
                
                temp.edge01 <- temp.ring.m $theta %>% max()
                temp.edge02 <- temp.sec.big$theta %>% min()
                
                temp.dist01 <- temp.ring.m $dist[which.max(temp.ring.m $theta)]
                temp.dist02 <- temp.sec.big$dist[which.min(temp.sec.big$theta)]
                
              }else{
                
                # S.Flip
                
                temp.edge01 <- temp.ring.m$theta %>% max()
                temp.edge02 <- Edge.02
                
                temp.dist01 <- temp.ring.m$dist[which.max(temp.ring.m$theta)]
                temp.dist02 <- dist.02
                
              }
              
              temp.sec.rs <- Section.Morph(Edge.01  = temp.edge01,
                                           Edge.02  = temp.edge02,
                                           dist.01  = temp.dist01,
                                           dist.02  = temp.dist02,
                                           ref.ring = ring.2)
              
              # Plot Check
              if(.plotcheck == TRUE){
                
                points(x   = temp.sec.rs$theta, 
                       y   = temp.sec.rs$dist, 
                       col = "red3", cex = 0.3)
                
              }
              
              #' [Compile the sectional morph products]
              temp.sec.k <- rbind(temp.sec.ls, temp.sec.rs)
              
            }
            
            # Not Using else just to emphasize
            #' [Case02_small ring segment group DO NOT EXISTS]
            if(purrr::is_empty(temp.ring.m$clst) == TRUE){
              
              #' [Add-in] : S.Flip
              if(S.Flip == FALSE){
                
                #'@note 
                # clst.L is the target ring segments 
                # clst.R is temp.sec.big at big.segment.pos[k]
                
                temp.edge01 <- Edge.01
                temp.edge02 <- temp.sec.big$theta %>% min()
                
                temp.dist01 <- dist.01
                temp.dist02 <- temp.sec.big$dist[which.min(temp.sec.big$theta)]
                
              }else{
                
                # S.Flip
                temp.edge01 <- temp.sec.big$theta%>% max()
                temp.edge02 <- Edge.02
                
                temp.dist01 <- temp.sec.big$dist[which.max(temp.sec.big$theta)]
                temp.dist02 <- dist.02
                
              }
              
              
              temp.sec.bs <- Section.Morph(Edge.01  = temp.edge01,
                                           Edge.02  = temp.edge02,
                                           dist.01  = temp.dist01,
                                           dist.02  = temp.dist02,
                                           ref.ring = ring.2)
              
              # Plot Check
              if(.plotcheck == TRUE){
                
                points(x   = temp.sec.bs$theta, 
                       y   = temp.sec.bs$dist, 
                       col = "red3", cex = 0.3)
                
              }
              
              
              #' [Compile the sectional morph products]
              temp.sec.k <- temp.sec.bs
              
            }
            
            #' [Compile the big ring Segment]
            if(Big.Group.Mid == TRUE){
              
              temp.sec.big <- temp.sec.big[, c("clst", "theta", "dist")]
              temp.sec.k   <- rbind(temp.sec.k, ..mid.big)
              
              # Disable Signal
              Big.Group.Mid <- FALSE
              
            }else if(Big.Group == TRUE){
              
              temp.sec.big <- .temp.big[, c("clst", "theta", "dist")]
              temp.sec.k   <- rbind(temp.sec.k, temp.sec.big)
              temp.sec.k   <- rbind(temp.sec.k, temp.group.morph)
              
              # Disable Signal
              Big.Group <- FALSE
              
            }else{
              
              temp.sec.big <- temp.sec.big[, c("clst", "theta", "dist")]
              temp.sec.k   <- rbind(temp.sec.k, temp.sec.big)
              
            }
            
            temp.ring      <- rbind(temp.sec.big, temp.small)
            temp.ring$clst <- temp.ring$clst %>% as.character()
            ring.keep[[k]] <- temp.ring
            
            # Continue with:
            #' [Step03_Compile & Update temp.ring.sc to process the loop]
            
          }
          
          
          
          ## (04)_Compile & Update ----
          #' [Step03_Compile & Update temp.ring.sc to process the loop]
          # Compile
          temp.sec.k <- rbind(temp.ring.m, temp.sec.k)
          temp.sec.k$clst <- temp.sec.k$clst %>% as.character()
          morph.list[[k]] <- temp.sec.k
          
          
          # Simplified the side pattern of section.k for Edge.01 update
          #' [Add-in] : S.Flip
          if(S.Flip == FALSE){
            
            # Edge.01 update
            #'@description 
            #' re-write into the logic as "S.Flip"
            .r.morph.k <- 
              morph.list %>% 
              dplyr::bind_rows() %>% 
              rbind(ring.i)
            
            temp.smdf  <- 
              .r.morph.k[which(.r.morph.k$theta < (Edge.02-.M.Step)),]
            
            if(nrow(temp.smdf) > 50*ang.npts){ # 10 means 10 degrees
              temp.smdf <- temp.smdf[order(-temp.smdf$theta),]
              temp.smdf <- temp.smdf[c(1:(50*ang.npts)),]
            }
            
            .range <- temp.smdf$theta %>% range() %>% diff()
            if(.range < 0.5){sp.df <- 5}else{sp.df <- 10}
            
            temp.sm.k  <- smooth.spline(x  = temp.smdf$theta,
                                        y  = temp.smdf$dist,
                                        df = sp.df)
            temp.sp <-
              predict(temp.sm.k,
                      x = temp.smdf$theta) %>%
              as.data.frame()
            names(temp.sp) <- c("theta", "dist")
            
            
            # Update Edge.01
            Edge.01 <- temp.sp$theta %>% max() %>% round(digits = .digits)
            dist.01 <- temp.sp$dist[which.max(temp.sp$theta)]
            
            
          }else{
            
            # Edge.02 update
            #'@description 
            #' The strategy is to add the Flip-searched big-segments
            #' in "ring.i" to allow change the "edge.02" (with Tag).
            #' Later, it is much easier to filter this out.
            .s.flip.k      <-  temp.sec.k
            .s.flip.k$clst <- "S.Flip"
            ring.i         <- rbind(ring.i, .s.flip.k)
            
            temp.smdf <- ring.i[which(ring.i$theta > Edge.01),]
            if(nrow(temp.smdf) > 50*ang.npts){ # 10 means 10 degrees
              temp.smdf <- temp.smdf[order( temp.smdf$theta),]
              temp.smdf <- 
                temp.smdf[c(1:(50*ang.npts)),]
            }
            .range <- temp.smdf$theta %>% range() %>% diff()
            if(.range < 0.5){sp.df <- 5}else{sp.df <- 10}
            temp.sm.k  <- smooth.spline(x  = temp.smdf$theta,
                                        y  = temp.smdf$dist,
                                        df = sp.df)
            temp.sp <-
              predict(temp.sm.k,
                      x = temp.smdf$theta) %>%
              as.data.frame()
            names(temp.sp) <- c("theta", "dist")
            
            
            # Update Edge.02
            Edge.02 <- temp.sp$theta %>% min() %>% round(digits = .digits)
            dist.02 <- temp.sp$dist[which.min(temp.sp$theta)]
            
            # disable Flag "S.Flip"
            S.Flip <- FALSE
            
          }
          
          # Update the loop
          ring.m <- Section.Morph(Edge.01  = Edge.01,
                                  Edge.02  = Edge.02,
                                  dist.01  = dist.01,
                                  dist.02  = dist.02,
                                  ref.ring = ring.2)
          
          temp.ring.sc <-
            ring.ck %>% 
            segment.filter(Edge.01 = Edge.01,
                           Edge.02 = Edge.02,
                           ring.segments = .) %>% 
            ring.join(., ring.segment = ring.m)
          
          head(temp.ring.sc)
          
          if(purrr::is_empty(temp.ring.sc$clst)){
            
            temp.ring.sc.sum <- temp.ring.sc.sum %>% as.data.frame()
            temp.ring.sc.sum <- temp.ring.sc.sum[vector(),]
            
          }else{
            
            temp.ring.sc.sum <-
              temp.ring.sc %>% 
              dplyr::mutate(delta = abs(dist - dist.Morph)) %>% 
              dplyr::group_by(clst) %>% 
              dplyr::arrange(delta) %>% 
              dplyr::summarise(npts      = n(),
                               ..npts    = round(npts*0.1),
                               m.delta   = mean(delta[1:..npts]),
                               sd.delta  = sd  (delta[1:..npts]),
                               # m.Ang.pos = npts / nrow(ring.2),
                               max.theta = max(theta),
                               min.theta = min(theta),
                               dist.l    = dist[which.min(theta)],
                               dist.r    = dist[which.max(theta)]) %>% 
              dplyr::select(-..npts) %>%
              dplyr::filter(., m.delta < 5) %>% 
              dplyr::left_join(., lookup.SUM, by = "clst") %>% 
              segment.clean(., dist.thres = 3) %>% 
              segment.clean(., dist.thres = 2, segment.size = segment.size) %>% 
              segment.overlap(., segment.size = segment.size) %>% 
              as.data.frame()
            
          }
          
          temp.ring.sc.sum
          
          # Update k for saving data
          k = k + 1
          
        }
        
        # -------------------------------------------------------------------- #
        # ___K_END ----------------------------------------------------------- 
        # -------------------------------------------------------------------- #
        
        ## 2.3 Compile Found Segments ----
        #'[Compile tree ring structure based on segment i (Ring.i)]
        ring.keep  <- dplyr::bind_rows(ring.keep)
        ring.morph <- dplyr::bind_rows(morph.list)
        
        if(purrr::is_empty(ring.keep$clst)){
          
          #' [No Segment Found]
          #' Segments are too small to be accepted
          #' The desired temp.rest == ring.m
          ring.keep  <- ring.ck[vector(),]
          ring.morph <- ring.ck[vector(),]
          print(paste("Segments are too small to be accepted"),)
          print(paste("ID: ", temp.1st.sum$clst, sep = ""),)
          
          # Morph
          temp.rest <- ring.m
          
        }else{
          
          # Last position in the morph list 
          # This is the one has a gap to Edge.02 (Not ANYMORE...)
          # .pos      <- length(morph.list)
          # .last.sec <- morph.list[[.pos]]
          
          temp.edge01 <- Edge.01
          temp.edge02 <- Edge.02
          
          temp.dist01 <- dist.01
          temp.dist02 <- dist.02
          
          # Morph
          #'@concept
          #' Special TRY
          #' Directly Modify "ring.2" when "middle.II == TRUE"
          temp.rest <- Section.Morph(Edge.01  = temp.edge01,
                                     Edge.02  = temp.edge02,
                                     dist.01  = temp.dist01,
                                     dist.02  = temp.dist02,
                                     ref.ring = ring.2)
          
        }
        
        # Plot Check
        if(.plotcheck == TRUE){
          
          points(temp.rest$theta, temp.rest$dist, col = "red3", cex = 0.3)
          
        }
        
      }else{
        
        #' [No Segment Found]
        ring.keep  <- ring.ck[vector(),]
        ring.morph <- ring.ck[vector(),]
        
        #' Jump to:
        #'[Refine tree ring structure based on segment i (Ring.i)]
        
        # Morph from Original Edges and Dist
        # temp.edge01 <- Edge.01
        # temp.edge02 <- Edge.02
        # 
        # temp.dist01 <- dist.01
        # temp.dist02 <- dist.02
        # 
        # # Morph
        # temp.rest <- Section.Morph(Edge.01  = temp.edge01,
        #                            Edge.02  = temp.edge02,
        #                            dist.01  = temp.dist01,
        #                            dist.02  = temp.dist02,
        #                            ref.ring = ring.2)
        
        #' [Thus temp.rest == ring.m]
        # Morph
        temp.rest <- ring.m
        
        # Plot Check
        if(.plotcheck == TRUE){
          
          points(temp.rest$theta, temp.rest$dist, col = "red3", cex = 0.3)
          
        }
        
      }
      
      #'[Add-in] : S.Flip
      #' Filter the Tagged segments in ring.i ("S.Flip")
      ring.i <- 
        ring.i %>% 
        dplyr::mutate(clst = clst %>% as.character())%>% 
        dplyr::select(c("clst", "theta", "dist")) %>% 
        dplyr::filter(clst != "S.Flip")
      
      # Compile the morph data frame
      ring.keep  <- rbind(ring.keep[,c("clst", "theta", "dist")], 
                          ring.i   [,c("clst", "theta", "dist")])
      ring.morph <- rbind(ring.morph, temp.rest)
      ring.morph <- rbind(ring.morph, ring.i[,c("clst", "theta", "dist")])
      
      # Check 
      #'@note
      #' In case, in the future, 
      #' there's time to fix this issue, it happens in the ring.morph that 
      #' there are segments with ID overlapped with the morphed segments.
      if(nrow(ring.morph) != nrow(ring.2)){
        
        message("PHASE I - Unequal Length: segment.GroupRefine()")
        # stop("The morph ring has unequal length to the morph base")
        
        # Quick Solution:
        ring.morph <- ring.morph %>% segment.GroupRefine()
        if(nrow(ring.morph) != nrow(ring.2)){stop("...Failled")}
        
      }
      
      # Convert back to 0-2pi region
      if(Padding == TRUE){
        
        padpos.keep  <- which(ring.keep $theta > 2*.pi)
        padpos.morph <- which(ring.morph$theta > 2*.pi)
        
        ring.keep $theta[padpos.keep]  <- ring.keep $theta[padpos.keep]  -2*.pi
        ring.morph$theta[padpos.morph] <- ring.morph$theta[padpos.morph] -2*.pi
        
        # Restore
        Padding <- vector()
        
      }
      
      # Restore
      ring.2  <- .ring.2 
      ring.ck <- .ring.ck
      
      # i-loop control
      i = i + 1
      
      ## 2.4 Refine Compiled Structure ----
      #'[Refine tree ring structure based on segment i (Ring.i)]
      #' @description
      #' Update ring.keep & ring.morph to include tiny segments
      
      # Check
      # stop("check pick.sum")
      #' [Hold Ring Segments Control]
      # l-loop control
      hold.check.morph <- 
        ring.morph %>% dplyr::filter(clst != "Morph") %>% nrow()
      
      # Use "ring.morph" to extract the small segments
      #'@note choose a "morph ring"
      if(hold.check.morph / nrow(ring.2) > 0.5 &
         abs(temp.edge01 - temp.edge02) < 1 ){
        
        temp.morph <- 
          ring.keep %>% 
          segment.smooth(ring.ref = ring.2,
                         coverage = 0.7,
                         name.write = "Morph",
                         Morph.step = ReSample.Size)
        
      }else{
        
        temp.morph <- 
          ring.morph[, c("theta", "dist")]
        
      }
      temp.morph$clst <- c("Morph")
      
      temp.sum   <- 
        ring.ck %>% 
        ring.join(., ring.segment = temp.morph) %>% 
        dplyr::group_by(clst) %>% 
        dplyr::mutate(delta = abs(dist - dist.Morph)) %>% 
        dplyr::arrange(delta) %>% 
        dplyr::summarise(npts      = n(),
                         m.delta   = mean(delta[1:round(0.5*npts)]),
                         sd        = sd  (delta[1:round(0.5*npts)]),
                         max.theta = max(theta),
                         min.theta = min(theta)) %>% 
        dplyr::filter(., m.delta < 2)
      
      pick.clst <- setdiff(temp.sum$clst, unique(ring.keep$clst))
      pick.sum  <- 
        temp.sum %>% 
        dplyr::filter(., clst %in% pick.clst) %>% 
        segment.overlap(., overlap.check = "All") %>% 
        dplyr::mutate(det = 
                        min(ring.i$theta) > min.theta & 
                        min(ring.i$theta) < max.theta) %>% 
        dplyr::filter(., det == FALSE)
      ring.pick <- 
        ring.ck %>% filter(., clst %in% pick.sum$clst)
      
      # Update "ring.keep" and "ring.morph"
      ring.keep  <- rbind(ring.keep, ring.pick[,names(ring.keep)])
      ring.morph <- ring.substitute(ring.1 = ring.morph,
                                    ring.segment = ring.pick)
      
      # Check
      # stop("pick.sum Check")
      
      # Check 
      if(nrow(ring.morph) != nrow(ring.2)){
        
        message("PHASE II - Unequal Length: segment.GroupRefine()")
        # stop("The morph ring has unequal length to the morph base")
        
        # Quick Solution:
        ring.morph <- ring.morph %>% segment.GroupRefine()
        if(nrow(ring.morph) != nrow(ring.2)){stop("...Failled")}
        
      }
      
      
      
      # 3_Hold Ring Segments Control ----
      #'[Hold Ring Segments Control]
      # l-loop control
      hold.check.morph <- 
        ring.morph %>% dplyr::filter(clst != "Morph") %>% nrow()
      
      #'[Completeness Control of the segments found]
      if(hold.check.morph / nrow(ring.2) > 0.5 &
         nrow(ring.keep) / nrow(ring.2) > sapw.control[ii] ){
        
        # & nrow(ring.keep) < nrow(ring.2)*0.7 ## Not use
        
        message("# Segment Holding #")
        
        #'[01:Preparation]
        #'@note
        #' (1) .hold.morph: smooth.splined ring.keep for further comparison
        #' (2) .hold.dist: distance to ring.2
        
        #'(1) .hold.morph
        #'@note Fill the huge gaps within ring.keep
        .hold.morph <- 
          ring.keep %>%
          segment.smooth(ring.ref = ring.2,
                         coverage = 0.7,
                         name.write = "ref",
                         Morph.step = ReSample.Size)
        
        #'(2) .hold.dist
        ring.2      <- ring.2 %>% dplyr::arrange(theta)
        .hold.dist  <- (.hold.morph$dist - ring.2$dist) %>% mean()
        
        # Ensure Data Structure
        ring.keep$theta <- ring.keep$theta %>% round(digits = .digits)
        
        # Check
        # stop("Hold Ring Segments Control...")
        
        #'[02:Compare & Consider the hold Ring Segments]
        if(l > 1){
          
          #'[Det_Check: How many unique cluster in ring.keep]
          #'@note 
          #' [Situation:]
          #' (1) hold.det < 0.1
          #' (2) hold.det > 0.1
          #' [Meaning:]
          #' [hold.det < 0.1]
          #' This means the difference in cluster "Number" is minimum.
          #' Function segment.DetHold considers segment SIZE & Number of clsts.
          #' [Decision:]
          #' (1) Huge number of disagreements -> different ring structure.
          #' (2) All segments will contribute voting process to decide which 
          #'     structure is the "Most" robust structure agreed by "Most"
          #'     considered tree ring segments.
          #' - Special Check-III when .hold.sum has segment > 0.1*nrow(ring.2)
          #'   (This is to avoid situation that a big segment is acquired...)
          
          Hold.det.list <- segment.DetHold(target.ring    = ring.keep,
                                           target.morph   = .hold.morph,
                                           hold.id        = hold.id, 
                                           hold.morph     = hold.morph,
                                           check.segments = ring.ck,
                                           ReSample.Size  = ReSample.Size)
          hold.det      <- Hold.det.list$hold.det
          hold.sum.list <- Hold.det.list$hold.sum.list
          
          #'[Check-III]
          #' if(any(hold.det < 0.1 &
          #'        purrr::map_lgl(hold.sum.list, 
          #'                       ~any(.x$npts > 0.1*nrow(ring.2))))){
          #'   
          #'   #'(1) ck.pos
          #'   ck.pos <- 
          #'     which((hold.det < 0.1 & 
          #'              purrr::map_lgl(hold.sum.list, 
          #'                             ~any(.x$npts > 0.1*nrow(ring.2)))) == TRUE)
          #'   #'(2) .replace
          #'   .replace <- which(hold.det <= 0.1)
          #'   if(length(.replace) > 1){.replace <- which.min(hold.det)}
          #'   
          #'   # ck.pos vs. .replace
          #'   if(identical(ck.pos, .replace)){
          #'     
          #'     stop("Huge Segment disagreements ... Check")
          #'     
          #'     # Check specific hold tree ring structure at ck.pos
          #'     ck.hold <- hold.keep[[ck.pos]]
          #'     ck.hold$clst <- "hold.ck"
          #'     
          #'     ck.keep <- ring.keep
          #'     ck.keep$clst <- "keep.ck"
          #'     
          #'     plot(ck.keep$theta, ck.keep$dist, cex = 0.3)
          #'     points(ck.hold$theta, ck.hold$dist, cex = 0.3, col = "darkred")
          #'     
          #'   }
          #'   
          #' }
          
          #'[03:Determine Current Found Ring]
          if(all(hold.det > 0.1) == TRUE){
            
            # Signals For Existence of Another Similar Options
            .replace <- which(hold.det < 0.8)
            if(purrr::is_empty(.replace) != TRUE){
              
              if(length(.replace) > 1){.replace <- which.min(hold.det)}
              
              message("0.1 < hold.det < 0.8")
              print(paste("Most Similar to Hold.ID =", .replace))
              print(paste("hold.det =", 
                          round(hold.det[.replace], digits = 2)))
              
            }
            
            # Holding Signal
            message("Segment Hold")
            
            # Check
            # if(hold.det != 1){stop("Check... 0.8 < hold.det < 1")}
            print(hold.det)
            
            hold.id   [[l]] <- ring.keep$clst %>% unique()
            hold.keep [[l]] <- ring.keep
            hold.morph[[l]] <- .hold.morph
            hold.value      <- append(hold.value, nrow(ring.keep)/nrow(ring.2))
            hold.dist       <- append(hold.dist, .hold.dist)
            hold.clst.ck    <- append(hold.clst.ck, l)
            
            # l-loop control
            l = l + 1
            
          }else{
            
            #'@note
            #' Thanks to the previous condition,
            #' only really close segments (det < 0.1) will be handled here.
            #'@details 
            #' (1) Find overlapped hold segments 
            #' (2) Compare / Replace the one the overlapped the most.
            
            #'[_Similar Position in Hold*]
            #' Hold Positions where ring.keep is similar to
            .replace <- which(hold.det <= 0.1)
            if(length(.replace) > 1){.replace <- which.min(hold.det)}
            
            # Check-III
            # ck.hold      <- hold.keep[[.replace]]
            # ck.hold$clst <- "hold.ck"
            # ck.keep      <- ring.keep
            # ck.keep$clst <- "keep.ck"
            # plot(ck.keep$theta, ck.keep$dist, cex = 0.3)
            # points(ck.hold$theta, ck.hold$dist, cex = 0.3, col = "darkred")
            
            #'[_Considered Segments]
            #' ring.keep segments
            #' unique segments within ring.keep
            .keep.id  <- ring.keep$clst %>% unique()
            diff.clst <- dplyr::setdiff(.keep.id, hold.id[[.replace]])
            
            #'[03-01 Hold Scenario <= 0.1 Determination]
            if(purrr::is_empty(hold.sum.list[[.replace]]$clst) != TRUE){
              
              #'[Preparation:]
              .keep.id
              ..id <- append(.keep.id, hold.id [[.replace]]) %>% unique()
              ..id.morph <-
                ring.ck %>% 
                dplyr::filter(clst %in% ..id) %>% 
                dplyr::mutate(theta = round(theta, digits = .digits)) %>% 
                segment.smooth(ring.ref = ring.2,
                               coverage = 0.7,
                               name.write = "ref",
                               Morph.step = ReSample.Size)
              
              #'[For_ring.keep:]
              # Found ID
              .keep.id 
              
              # Difference in Clusters
              diff.clst <- dplyr::setdiff(.keep.id, hold.id[[.replace]])
              
              # Ensure Data Structure
              hold.diff       <- ring.ck %>% dplyr::filter(clst %in% diff.clst)
              hold.diff$theta <- hold.diff$theta %>% round(digits = .digits)
              .hold.ref       <- ..id.morph 
              .hold.ref$theta <- .hold.ref$theta %>% round(digits = .digits)
              
              #'Overlap-Check Target: hold.keep[[.replace]]
              .hold.keep <- hold.keep[[.replace]]
              .hold.keep$theta <- .hold.keep$theta %>% round(digits = .digits)
              
              #'[keep.sum] + [keep.overlap]
              .keep <-
                hold.diff %>% 
                ring.join(., .hold.ref) %>% 
                dplyr::mutate(delta = dist - dist.ref)
              keep.sum <- 
                .keep %>% 
                dplyr::group_by(clst) %>% 
                dplyr::mutate(theta   = round(theta, digits = .digits),
                              overlap = theta %in% .hold.keep$theta) %>% 
                dplyr::summarise(npts = n(),
                                 det  = mean(delta),
                                 sd   = sd(delta),
                                 overlap = sum(overlap))
              #'[keep.det]
              keep.tol <- .keep$delta %>% abs() %>% mean()
              
              #'[For_hold.keep:]
              # Hold ID
              hold.id [[.replace]]
              
              # Difference in Clusters
              diff.clst <- dplyr::setdiff(hold.id[[.replace]], .keep.id)
              
              # Ensure Data Structure
              hold.diff       <- ring.ck %>% dplyr::filter(clst %in% diff.clst)
              hold.diff$theta <- hold.diff$theta %>% round(digits = .digits)
              .hold.ref       <- ..id.morph 
              .hold.ref$theta <- .hold.ref$theta %>% round(digits = .digits)
              
              #'Overlap-Check Target: ring.keep
              ring.keep$theta <- ring.keep$theta %>% round(digits = .digits)
              
              #'[hold.sum] + [hold.overlap]
              .hold <- 
                hold.diff %>% 
                ring.join(., .hold.ref) %>% 
                dplyr::mutate(delta = dist - dist.ref)
              hold.sum <- 
                .hold %>% 
                dplyr::group_by(clst) %>% 
                dplyr::mutate(theta   = round(theta, digits = .digits),
                              overlap = theta %in% ring.keep$theta) %>% 
                dplyr::summarise(npts = n(),
                                 det  = mean(delta),
                                 sd   = sd(delta),
                                 overlap = sum(overlap))
              #'[hold.det]
              #'@note
              #' If .hold EMPTY, means all segments of hold.id.replace are
              #' ALSO in ring.keep. Hence, report: hold.tol = 0;
              #' Better performed ring.keep will be handled in "Justification"
              if(purrr::is_empty(.hold$clst)){hold.tol <- 0
              }else{hold.tol <- .hold$delta %>% abs() %>% mean()}
              
              #'[Justification:]
              .keep.sum <- keep.sum %>% dplyr::filter(overlap == 0)
              .hold.sum <- hold.sum %>% dplyr::filter(overlap == 0)
              if(purrr::is_empty(.keep.sum$clst) != TRUE &
                 min(hold.tol, keep.tol) < 5){
                
                # Report
                message("Similarity Justified ... Merge")
                
                # replacement aiming at: .keep.sum
                # (Keep if delta < 3) | No-overlap
                #'@note 
                #' the segment difference are really small, Hence,
                #' if there is no overlap, it will be support by the 
                #' the neighboring segments.(*)
                
                # ------------------------------------------------------------ #
                #'@special Check
                #' If there are any segments extremely large,
                #' stop and check...
                # if(any(.keep.sum$npts > 2000)){
                #   stop("Check Replace Justification Argument ...keep")}
                # ------------------------------------------------------------ #
                
                .keep.sum <- 
                  keep.sum %>% 
                  dplyr::filter((abs(det) + sd) < 3 & 
                                  (overlap == 0 | npts < 1000))
                
                diff.clst <- .keep.sum$clst
                hold.value[.replace] <- Inf
                
              }else if(hold.value[.replace] < nrow(ring.keep)/nrow(ring.2) &
                       min(hold.tol, keep.tol) < 5){
                
                # Report
                message("Similarity Justified ... Replace")
                
                # replacement aiming at: .hold.sum
                # (Keep if delta < 3) | No-overlap
                #'@note 
                #' the segment difference are really small, Hence,
                #' if there is no overlap, it will be support by the 
                #' the neighboring segments.(*)
                
                # ------------------------------------------------------------ #
                #'@special Check
                #' If there are any segments extremely large,
                #' stop and check...
                # if(any(.hold.sum$npts > 2000)){
                #   stop("Check Replace Justification Argument ... hold")}
                # ------------------------------------------------------------ #
                
                .hold.sum <-
                  hold.sum %>% 
                  dplyr::filter((abs(det) + sd) < 3 & 
                                  (overlap == 0 | npts < 1000))
                
                # Signal PASS
                hold.id[[.replace]] <- .hold.sum$clst
                
              }else{
                
                print(keep.sum)
                
                # No Need to Update
                hold.det[.replace] <- 0
                
              }
              
            }
            
            #'[03-02 Hold Scenario <= 0.1 Adjustments]
            # Processing ...
            if(hold.det[.replace] != 0){
              
              if(hold.value[.replace] < nrow(ring.keep)/nrow(ring.2)){
                
                # Report
                message("Replace Hold-Segment No.", .replace)
                
                #'[Update_hold_groups:]
                #'[.hold.id]
                .hold.id <- append(.keep.id, hold.id[[.replace]]) %>% unique()
                
                #'[.hold.morph]
                .ring.keep  <- ring.ck %>% dplyr::filter(clst %in% .hold.id)
                .hold.morph <- 
                  .ring.keep %>% 
                  segment.smooth(ring.ref = ring.2,
                                 coverage = 0.7,
                                 name.write = "ref",
                                 Morph.step = ReSample.Size)
                
                #'[.hold.dist]
                ring.2      <- ring.2 %>% dplyr::arrange(theta)
                .hold.dist  <- (.hold.morph$dist - ring.2$dist) %>% mean()
                
              }else{
                
                # Report
                message(paste(c("Add New Segment in Hold No.", .replace,
                                "with ID.", diff.clst), collapse = " "))
                
                #'[Update_hold_groups:]
                #'[.hold.id]
                .hold.id <- append(hold.id[[.replace]], diff.clst)
                
                #'[.hold.morph]
                .ring.keep  <- ring.ck %>% dplyr::filter(clst %in% .hold.id)
                .hold.morph <- 
                  .ring.keep %>% 
                  segment.smooth(ring.ref = ring.2,
                                 coverage = 0.7,
                                 name.write = "ref",
                                 Morph.step = ReSample.Size)
                
                #'[.hold.dist]
                ring.2      <- ring.2 %>% dplyr::arrange(theta)
                .hold.dist  <- (.hold.morph$dist - ring.2$dist) %>% mean()
                
              }
              
              # Update
              hold.id   [[.replace]] <- .hold.id
              hold.keep [[.replace]] <- .ring.keep
              hold.morph[[.replace]] <- .hold.morph
              hold.value [.replace]  <- nrow(.ring.keep)/nrow(ring.2)
              hold.dist  [.replace]  <- .hold.dist
              
            }else{ # hold.det[.replace] == 0
              
              # Report
              message("Segment Discarded: Not Better in Completeness")
              
            }
            
            # Update Vote (l-loop) Results
            hold.clst.ck <- append(hold.clst.ck, .replace)
            
          }
          
        }else{
          
          hold.id   [[l]] <- ring.keep$clst %>% unique()
          hold.keep [[l]] <- ring.keep
          hold.morph[[l]] <- .hold.morph
          hold.value      <- append(hold.value, nrow(ring.keep)/nrow(ring.2))
          hold.dist       <- append(hold.dist, .hold.dist)
          hold.clst.ck    <- append(hold.clst.ck, l)
          
          # l-loop control
          l = l + 1
          
        }
        
        
      }else{
        
        # Give the fail ones a record
        hold.clst.ck <- append(hold.clst.ck, 0)
        
      }
      
      #'[Full.segment Control]
      #'@description 
      #' Get ring segments in the 1st hold tree ring structure to 
      #' "contribute & vote" robust tree ring at similar position 
      #' based on their morph-searched results.
      if(hold.check.morph / nrow(ring.2) > 0.5 &
         nrow(ring.keep) / nrow(ring.2) > sapw.control[ii] ){
        
        if(Full.segment == FALSE){
          
          message("# Determine Ring Segment ...")
          
          # Save the current one, and add the additional search list
          .ring.SUM <- ring.SUM
          
          # Create the corresponding "ring.SUM"
          ring.SUM <-
            ring.ck %>% 
            dplyr::filter(clst %in% hold.id[[1]]) %>% 
            ring.join(., ring.segment = ring.2) %>% 
            dplyr::group_by(clst) %>%
            dplyr::summarise(delta = (dist - dist.ref) %>% mean()) %>% 
            dplyr::left_join(., lookup.SUM, by = "clst") %>% 
            dplyr::filter(., segment.l > 0.5*segment.size) %>%
            dplyr::arrange(delta)
          
          .pos <- which(ring.SUM$clst == Ring.ID.i)
          ring.SUM <- 
            rbind(.ring.SUM[c(1:(i-1)),], ring.SUM[-.pos,])
          
          # update ni
          ni = nrow(ring.SUM)
          
        }
        
        print(Ring.ID.i)
        
        # Report Segment Found
        Full.segment   <- TRUE
        
      }
      
      # Check "hold.clst.ck"
      if(length(hold.clst.ck) != (i-1)){stop("Check hold.clst.ck")}
      
      # Final-below_dist-Range
      #'@special
      #' As ring.SUM is depending on the 1st Full.segment
      #' And after ni == nrow(ring.SUM),
      #' The "vote" is completed.
      #' However, there might be segments large enough and in-between 
      #' current "voted" ring and ring.2,
      #' Hence, below offers this opportunity to capture them.
      if(i > ni & purrr::is_empty(hold.morph)){
        
        # Check
        # stop("ni check")
        
        ring.SUM <-
          ring.ck %>% 
          ring.join(., ring.segment = ring.2) %>% 
          dplyr::group_by(clst) %>%
          dplyr::summarise(delta = (dist - dist.ref) %>% mean()) %>% 
          dplyr::filter(., delta > 1) %>% 
          dplyr::left_join(., lookup.SUM, by = "clst") %>% 
          dplyr::filter(., segment.l > segment.size) 
        
        if(purrr::is_empty(ring.SUM$clst) != TRUE){
          
          .ring.sum <-
            ring.SUM %>% 
            dplyr::arrange(delta)
          
          ring.SUM <- 
            .ring.sum[c(1:(ni + 1)),] %>% 
            dplyr::filter(is.na(clst)!=TRUE)
          
          if(ni == nrow(ring.SUM)){
            
            # Signal
            ring.SUM <- ring.SUM %>% dplyr::slice(0)
            message("Nothing Left for Morph-Search ...")
            
          }else{
            
            message("Enlarge Search range from ring.2...")
            
          }
          
        }else{message("Nothing Left for Morph-Search ... Size")}
        
        # Update
        ni <- nrow(ring.SUM)
        
        
      }else if(i > ni){
        
        .pos <- which(hold.clst.ck != 0)
        vote.pos <- 
          sort(table(hold.clst.ck[.pos]),decreasing=TRUE)[1] %>% 
          names() %>% as.numeric()
        vote.hold <- 
          sort(table(hold.clst.ck[.pos]),decreasing=TRUE)[1] %>% 
          as.numeric()
        vote.ring <- hold.morph[[vote.pos]]
        vote.ring$clst <- ring.2$clst
        
        # The REST possibilities:
        temp.ring.sum <- 
          ring.ck %>% 
          ring.join(., ring.segment = vote.ring) %>% 
          dplyr::group_by(clst) %>%
          dplyr::summarise(delta = (dist - dist.ref) %>% mean()) %>% 
          dplyr::left_join(lookup.SUM, by = "clst") %>% 
          dplyr::filter(segment.l > 0.5*segment.size) %>% 
          dplyr::filter(delta < 1)
        
        .clst <- dplyr::setdiff(temp.ring.sum$clst, ring.SUM$clst)
        temp.ring.sum <- temp.ring.sum %>% dplyr::filter(clst %in% .clst)
        
        if(purrr::is_empty(temp.ring.sum$clst) != TRUE){
          
          # Check
          # stop("i check")
          
          # Update ring.SUM
          ring.SUM <- rbind(ring.SUM, temp.ring.sum)
          
          # Update ni
          ni = nrow(ring.SUM)
          
        }else{
          
          # message
          vote.pos
          vote.hold
          message("Hold-segment is voted for pos: No.",
                  vote.pos, "/", length(hold.morph),
                  " (", vote.hold, "/", nrow(ring.SUM), ")")
          
          # Check
          vote.ck <- ring.SUM
          vote.ck$hold.vote <- hold.clst.ck
          vote.sum <-
            vote.ck %>% 
            dplyr::group_by(hold.vote) %>% 
            summarise(count    = n(),
                      seg.size = mean(segment.l),
                      delta    = mean(delta))
          
          print(vote.sum, n = nrow(vote.sum))
          print(vote.ck[, c("clst", "segment.l", "hold.vote")],
                n = nrow(vote.ck))
          
          # Check
          # stop("Check")
          
          # ------------------------------------------------------------------ #
          # Determine the voted structure
          #'@description
          #' [Keys:]
          #' (1) vot.pos: Which structure gets the maximum votes;
          #' (2) vot.det: How different are other structure to the vote.pos
          #' (3) hold.dist: The distance between considered structure and  
          #'     the reference tree ring structure (ring.2).
          #' The scenarios triggered by (1) means there is a unique max vote, 
          #' when disagreed, this means might be overlapping scenario. Use 
          #' vot.det to clarify the independence of voted structure.
          #' [Scenarios:]
          #' (1) max vote  + non-independent other structures -> unique
          #' (2) max vote  + independent other structures     -> 1-2 unique
          #' (3) even vote + non-independent other structures -> overlap
          #' (2) even vote + independent other structures     -> 1-2 unique
          #' [Respective_Results:]
          #' (1) Remain vote.pos -> .update == TRUE
          #' (2) independent flowchart:
          #'     check hold.dist -> vote.pos.min   -> 1 unique, .update == TRUE
          #'     else -> check vote.ratio > 0.5    -> 2 unique, 
          #'     .update & .update.2 == TRUE; else -> 1 unique,
          #'     pick min.hold.dist, .update == TRUE
          #' (3) Overlap -> pick min.hold.dist
          #' (4) independent flowchart:
          #'     check hold.dist -> vote.pos.min   -> 1 unique, .update == TRUE
          #'     else -> check vote.ratio > 0.5    -> 2 unique, 
          #'     .update & .update.2 == TRUE; else -> 1 unique, 
          #'     pick min.hold.dist, .update == TRUE
          #' [Decision:]
          #' vote.pos -> 
          #' [vote.pos is min dist]
          #' :<Y> -> vote.pos; .update = TRUE
          #' :<N>
          #'   |
          #' [independent]
          #' :<Y> -> min.distance + vote.ratio
          #' :<N>
          #'   |
          #' [even vote]
          #' :<Y> -> overlap: min.distance
          #' :<N> -> No independent others + No even vote -> .update = TRUE
          #' 
          #'@detail
          #'[Independent_Flowchart]
          #'[01:]
          #' [vot.pos]
          #'  |
          #' [Any vot.det >= 0.9] :<N> -> [vot.pos]
          #' :<Y> 
          #'  |
          #' [which.min hold.dist] :<min = vot.pos> -> [vot.pos]
          #' :<min = pos.vote.det>=0.9>
          #'  |
          #' [vote.ratio > 0.5] :<Y> -> [vot.pos] + 
          #' :<N>                       [pos.vote.det>=0.9&min.hold.dist]
          #'  |
          #' [pos.vote.det>=0.9 & min.hold.dist]
          #'
          #'[02:]
          #' [Picked Structure]: [pos.vote.det>=0.9&min.hold.dist]
          #'  |
          #' [Any vot.det <= 0.1: n()>1] :<N> -> [Picked Structure]
          #' :<Y>
          #'  |
          #' [ANY [pos.vote.det>=0.9&min.hold.dist] %in% [vot.pos]] :<N>
          #' :<Y>                                                   |
          #'  |                                       |->{AP:Combine structures}
          #' {AP:Unique Structure} -------------------|             |
          #'                                                [Picked Structure]
          
          #'[Preparation:]
          #'@import 
          #'(1) vote.pos : which "hold structure" is voted;
          #'(2) vote.sum : vote results von Holding decision sections;
          #'(3) vote.hold : how many votes does the voted structure have;
          #'(3) hold.dist : distance between morphed structures & ref. ring;
          #'(4) hold.keep : kept ring segments;
          #'(5) hold.morph :morphed ring structure based on hold.keep.
          #'@eval 
          #'(1) vote.det : hold.id-comparison between vote.pos and others;
          #'(2) hold.dist
          #'(3) Overall Vote-Acquired Percentage
          #'@details 
          #'(1) 2 or more ring structures were detected in hold;
          #'(2) min hold.dist (dist to ref tree ring);
          #'(3) structures with same vote.
          
          #'[Det_Check]
          Vote.det.list <- 
            segment.DetHold(target.ring  = hold.keep [[vote.pos]],
                            target.morph = hold.morph[[vote.pos]],
                            hold.id      = hold.id, 
                            hold.morph   = hold.morph,
                            check.segments = ring.ck,
                            ReSample.Size  = ReSample.Size)
          vote.det      <- Vote.det.list$hold.det
          vote.sum.list <- Vote.det.list$hold.sum.list
          
          # plot.check
          # k = 1
          # plot(hold.keep[[k]]$theta, hold.keep[[k]]$dist, cex = 0.3)
          # points(hold.keep[[vote.pos]]$theta, hold.keep[[vote.pos]]$dist, 
          #        cex = 0.3, col = "lightblue")
          
          #'[Keys]
          vote.pos
          vote.det
          multiple_votes <- 
            vote.sum %>%
            dplyr::filter(hold.vote != 0, count == vote.hold) %>%
            nrow() > 1
          vote.pos.min <- 
            if (any(vote.det >= 0.9)) {
              hold.dist[vote.pos] < min(hold.dist[vote.det >= 0.9])
            } else {TRUE}
          
          # ------------------------------------------------------------------ #
          
          #'[Decision]
          if(all(vote.pos.min == TRUE)){
            
            # Signal
            message("[Det]: Determined by vote's distance to reference")
            
            
            
          }else if(any(vote.det >= 0.9)){
            
            #'[hold.det >= 0.9, Other independent ring structure]
            
            # Independent others (Tree Ring Structures)
            .ck.pos <- which(vote.det >= 0.9)
            
            # Signal
            message("[Det]: Other Independent Tree Ring Structure Found")
            print(.ck.pos)
            
            # More than 1 structures recognized:
            if(length(.ck.pos)>1){
              
              # Independent Positions:
              ind.pos <- append(vote.pos, .ck.pos) %>% sort()
              
              # Start Coding with which.max only return 1 value:
              .ck.pos <- which.max(vote.det) 
              
              # Update det value based on .ck.pos
              .det.list <- segment.DetHold(target.ring  = hold.keep [[.ck.pos]],
                                           target.morph = hold.morph[[.ck.pos]],
                                           hold.id      = hold.id, 
                                           hold.morph   = hold.morph,
                                           check.segments = ring.ck,
                                           ReSample.Size  = ReSample.Size)
              .det      <- .det.list$hold.det
              .sum.list <- .det.list$hold.sum.list
              
              # Re-select independent segments based on .ck.pos
              .pick.pos <- which(.det >= 0.9)
              .pick.ind <- append(.ck.pos, .pick.pos) %>% sort()
              
              # More than 1 structures recognized:
              if(identical(as.integer(ind.pos), as.integer(.pick.ind))){
                stop("[Note]: More than 1 structures recognized")}
              
            }
            
            # Final hold.dist Assure
            if(hold.dist[.ck.pos] < hold.dist[vote.pos]){
              
              # Unique clst structure of .ck.pos
              hold.id     [[.ck.pos]] <- 
                dplyr::setdiff(hold.id[[.ck.pos]], hold.id[[vote.pos]])
              hold.keep   [[.ck.pos]] <- 
                hold.keep [[.ck.pos]] %>% 
                dplyr::filter(clst %in% hold.id[[.ck.pos]])
              hold.morph  [[.ck.pos]] <- 
                hold.keep [[.ck.pos]] %>% 
                segment.smooth(ring.ref = ring.2,
                               coverage = 0.7,
                               name.write = "ref",
                               Morph.step = ReSample.Size)
              
              # if vote.pos has over 50% votes of the remaining votes
              vote.ck <- 
                hold.clst.ck[hold.clst.ck != 0] %>%
                table() %>%
                as.numeric() %>%
                .[.ck.pos]
              vote.remain <- length(which(hold.clst.ck != 0)) - vote.ck
              vote.ratio  <- vote.hold / vote.remain
              if(vote.ratio > 0.5){
                
                # Signal
                message("# Both Tree Ring Structures Acquired...")
                
                .update.2 <- TRUE
                .ring.keep  <- hold.keep [[.ck.pos]]
                .ring.morph <- hold.morph[[.ck.pos]]
                
              }else{
                
                # Signal
                message("# Inner Tree Ring Structures Acquired...")
                
                # Update vote.pos
                vote.pos <- .ck.pos
                
              }
              
            }else{
              
              #' *NORMALLY* This should not happen
              stop("[Note]: Other Independent Structure is further...")
              
            }
            
            
            
          }else if(multiple_votes == TRUE){
            
            #'[multiple_votes, Multiple options have same vote]
            
            # Signal
            message("[Det]: Overlapped Tree Ring Structures Encountered")
            print("[Note]: Overlaps are hard to judge: Potential Error")
            
            # Overlapped Tree Ring Structures + min hold.dist
            .ck.pos <- 
              vote.sum %>%
              dplyr::filter(hold.vote != 0, count == vote.hold) %>% 
              dplyr::pull(hold.vote) %>%
              {.[which.min(hold.dist[.])]}
            
            # Update vote.pos
            vote.pos <- .ck.pos
            
            
            
          }else{
            
            #'[Other Possibilities Removed, select vote.pos]
            
            # Signal
            message("[Det]: Determined by exclusion of other possibilities")
            
          }
          
          
          # Signal & Update
          .update <- TRUE
          ring.keep  <- hold.keep [[vote.pos]]
          ring.morph <- hold.morph[[vote.pos]]
          
        }
        
      }
      
      
      
      #'[Final-Result-Update for Morph & Search Object (Vote)]
      #' @description
      #' (1) Update the ring.2
      #' (2) Update the ring.ck to avoid strange twist (Remove found segments)
      if(.update == TRUE){
        
        # check j = 4
        # if(j == 3){stop()}
        
        #'[(1) Update the "ring.2"]
        # Save previous morph search base "ring.2"
        if(.update.2 == FALSE){
          
          ring.base.list [[j]] <- ring.2
          ring.keep.list [[j]] <- ring.keep
          ring.morph.list[[j]] <- ring.morph
          
          ring.2 <- ring.morph
          ring.2$clst <- "ref"
          ring.2 <- ring.2 %>% ring.spline()
          
          j = j + 1
          
        }else{ # .update.2 == TRUE
          
          # Signal-off
          .update.2 <- FALSE
          
          #'[Inner Tree Ring Structure]
          ring.base.list [[j]] <- ring.2
          ring.keep.list [[j]] <- .ring.keep
          ring.morph.list[[j]] <- .ring.morph
          
          ring.2 <- .ring.morph
          ring.2$clst <- "ref"
          ring.2 <- ring.2 %>% ring.spline()
          
          j = j + 1
          
          #'[Outer Tree Ring Structure]
          ring.base.list [[j]] <- ring.2
          ring.keep.list [[j]] <- ring.keep
          ring.morph.list[[j]] <- ring.morph
          
          ring.2 <- ring.morph
          ring.2$clst <- "ref"
          ring.2 <- ring.2 %>% ring.spline()
          
          j = j + 1
          
        }
        
        # A grouped Ring Segment ID-List
        group.ID <- hold.id[[vote.pos]]
        
        #'[(2) Update the searching list "Ring.SUM" + "ring.ck"]
        # I. "Ring.SUM" + "ring.ck"
        ring.ck$dist.ref <- NULL
        ring.ck <- 
          ring.ck %>% 
          ring.join(., ring.segment = ring.2)
        
        ring.SUM <-
          ring.ck %>% 
          dplyr::group_by(clst) %>%
          dplyr::summarise(delta = (dist - dist.ref) %>% mean()) %>% 
          dplyr::filter(., delta > 1) # Break for step product
        
        # Step stop to update ring.ck
        temp.clst <- dplyr::setdiff(ring.SUM$clst, group.ID)
        ring.ck   <- ring.ck %>% dplyr::filter(clst %in% temp.clst)
        
        ring.SUM <- # Continue processing
          ring.SUM %>% 
          dplyr::filter(clst %in% temp.clst) %>% 
          dplyr::left_join(., lookup.SUM, by = "clst") %>% 
          dplyr::filter(., segment.l > segment.size) 
        
        # .lookup.update <- ring.SUM
        
        if(purrr::is_empty(ring.SUM$clst) != TRUE){
          
          .ring.sum <-
            ring.SUM %>% 
            filter(., delta < search.range) %>% 
            arrange(delta)
          
          if(purrr::is_empty(.ring.sum$clst)){
            
            ring.SUM <- ring.SUM %>% dplyr::arrange(delta)
            ring.SUM <- ring.SUM[1,]
            
          }else{
            
            ring.SUM <- .ring.sum
            
          }
          
        }
        
        
        # II. Update ni for while-loop
        i  = 1
        ni = nrow(ring.SUM)
        
        # Reset .update
        .update       <- FALSE
        
        # Reset "hold-related" values
        l = 1
        hold.id        <- list()
        hold.keep      <- list()
        hold.morph     <- list()
        hold.value     <- vector()
        hold.dist      <- vector()
        hold.clst.ck   <- vector()
        
        Full.segment   <- FALSE
        
        # PAUSE
        # break
        
        # Check
        # if(j == 12){break}
        
      }
      
      print(paste("i = ", i-1, sep = ""))
      print(paste("ni = ", ni, sep = ""))
      print(nrow(ring.keep) / nrow(ring.2))
      print(hold.check.morph / nrow(ring.2))
      
      # Plot
      if(.plotcheck == TRUE){
        
        plot(ring.morph$theta,
             ring.morph$dist,
             col = as.factor(ring.morph$clst),
             cex = 0.3)
        
      }
      
      # i-loop control
      if(ni == 0){
        
        ring.base.list [[j]] <- ring.2
        break
        
      }
      
      # Debug
      # if(j == 10){break}
      
    }
    
    # ------------------------------------------------------------------------ #
    # ___I_END -------------------------------------------------------------- 
    # ------------------------------------------------------------------------ #
    
    # Output Complete-Ring Result
    ring.k.list <- append(ring.k.list, ring.keep.list)
    ring.m.list <- append(ring.m.list, ring.morph.list)
    ring.b.list <- append(ring.b.list, ring.base.list)
    
    # ii-loop Control
    ii = ii + 1
    if(ii > length.clst.trust){break}
    
  }
  
  # -------------------------------------------------------------------------- #
  # ___II_END -------------------------------------------------------------
  # -------------------------------------------------------------------------- #
  
  # Output
  Res.list <- list()
  Res.list$ring.k <- ring.k.list
  Res.list$ring.m <- ring.m.list
  Res.list$ring.b <- ring.b.list
  return(Res.list)
  
}
