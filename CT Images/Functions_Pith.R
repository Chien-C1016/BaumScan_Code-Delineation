# ---------------------------------------------------------------------------- #
# Functions of Pith Extraction #
# ---------------------------------------------------------------------------- #
# General ----
# Based on summaries of the clusters (which identified by dbscan pkg),
# The target is to sort clusters to see which ones are the most "inner" ones
# The "inner ones will be called "PITH".
#' (1) pith.CenterCheck: Summary Table for potential pith segments
#' (2) pith.GetPith: Get Pith positions (df / CImage)
#' (3) pith.ExtractPith: Fit Circle to Pith segment to get the center point
#' (4) dbscan.eps: Offer a self-defined eps suggested by dbscan package.

# ---------------------------------------------------------------------------- #
# Functions ----
## Tree Ring Segments Center Check
#'@description
#' Sort and order the ring segments based on their relative positions.
#' Aim at targeting the most central segment cluster.
#'@import
#  A dbscan processed data frame with cluster tag as "clst".
#'@export
#  Summary table
#'@source
pith.CenterCheck <- function(df.db){
  
  require(dplyr)
  
  #'[Preprocessing:]
  # Noise filtering
  #'@note Clust == 0 is noise
  df.db <- df.db %>% dplyr::filter(clst != 0)
  
  # Summarize
  df <- df.db %>% 
    dplyr::as_tibble() %>% 
    dplyr::group_by(clst) %>% 
    dplyr::summarise(.,
                     xmin = min(x),
                     xmax = max(x),
                     ymin = min(y),
                     ymax = max(y),
                     x.c  = mean(x),
                     y.c  = mean(y),
                     n = n())
  
  #'[Clean_breaking_ring_clusters:]
  #'@note
  #' (1) width 2 Length Ratio: df.check$ratio < 0.7
  #' (2) Convex Hull 2 circle ratio: df.check$a.ratio < 0.7
  #' (3) Segment too small: df.check$n < 100
  
  chull.area <- function(i){
    
    i <- as.matrix(i[, c("x", "y")])
    h  = chull(i)
    j  = c(h,h[1])
    ca = round((i[h,1]+i[j[-1],1])%*%diff(-i[j,2])/2)
    return(ca)
    
  }
  
  c.area <-
    df.db[, c("clst", "x", "y")] %>% 
    dplyr::group_by(clst) %>%
    dplyr::group_modify(~chull.area(.x) %>% data.frame(c.area = .))
  
  df.check <- 
    df %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(dx = xmax - xmin,
                  dy = ymax - ymin,
                  dr = mean(dx, dy)/2,
                  ratio = min(dx, dy) / max(dx, dy),
                  m.area = pi*dr^2) %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(., c.area, by = "clst") %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(a.ratio = min(m.area, c.area) / max(m.area, c.area))
  ck <- which(is.na(df.check$a.ratio))
  df.check$a.ratio[ck] <- 0
  
  # [Hot.Fix]
  # To improve the computational time,
  # We do not need any breaking clusters (not circular clusters)
  n.c <- which(df.check$ratio < 0.7 | 
                 df.check$a.ratio < 0.7 |
                 df.check$n < 100)
  df.break <- df.check[ n.c,]
  df.check <- df.check[-n.c,]
  
  # Adjust 2 data frames
  # df.check$ratio <- ratio[-n.c]
  # df.break$ratio <- ratio[ n.c]
  df.break$order <- NA
  df.break$core  <- NA
  
  # Report
  message("# Check Cluster's Relative Position #")
  
  #'[Sort:Relative Positions]
  # Start sorting for df.check
  core  <- vector()
  order <- vector()
  keep  <- vector()
  
  # [Note]
  # i controls the current comparison target cluster
  # j controls the feed in from S to O[,i]
  # k is just the total length of the comparison, length.S
  
  i = 1
  j = 1
  O <- df.check
  S <- O[-i,] # Fix.from S <- O[-(1:i),]
  k <- length(S$clst) # Fix.from length(O$clst)-1, much straight forward
  while(TRUE){
    
    if(all(c(S$xmin[j]:S$xmax[j]) %in% c(O$xmin[i]:O$xmax[i])) == T &
       all(c(S$ymin[j]:S$ymax[j]) %in% c(O$ymin[i]:O$ymax[i])) == T){
      
      j = j+1
      
      # Stop J-loop
      if(j > k){ # Fix.from cond.: max(length(S$clst))
        
        # [Note]
        # is.core is always F when loop end in here
        # That is, for all the clusters compared in S,
        # all clusters in S are in O[i]
        
        core [i] <- 0
        order[i] <- k - length(keep)
        
        # Report Status
        # message("Check cluster no. ", i, " / ", k)
        
        # i
        i = i+1
        
        # k
        k = k       # Fix.from k-1 (S is not shrinking anymore)
        
        # S
        S <- O[-i,] # Fix.from S <- O[-(1:i),]
        
        # Re-set the parameter
        j = 1
        keep <- vector()
        
      }
      
    }else{
      
      # [Note]
      # This means this j cluster either 
      # (1) has the same ring sequence, or
      # (2) is the outer ring sequence.
      
      keep <- append(keep, j)
      j = j+1
      
      # Stop J-loop
      if(j > k){ # Fix.from cond.: max(length(S$clst))
        
        # is.core??
        if(length(keep) == length(S$clst)){
          
          # [Note]
          # All clusters in S are 
          # (1) outside or 
          # (2) having same sequence number
          
          core[i]  <- 1
          
        }else{
          
          core[i]  <- 0
          
        }
        
        order[i] <- k - length(keep)
        
        # Report Status
        # message("Check cluster no. ", i, " / ", k)
        
        # i
        i = i+1
        
        # k
        k = k       # Fix.from k-1 (S is not shrinking anymore)
        
        # Stop Condition:
        if(i > (k+1)){ # Fix.from nrow(S) == 0
          # core[i]  <- 1
          # order[i] <- 0
          break}
        
        # S
        S <- O[-i,] # Fix.from S <- O[-(1:i),]
        
        # Re-set the parameter
        j = 1
        keep <- vector()
        
      }
      
    }
    
    # if(i > (k+1)){ # Fix.from nrow(S) == 0
    #   # core[i]  <- 1
    #   # order[i] <- 0
    #   break}
    
  }
  
  #'[Decision: Base on the ordered tree ring segments]
  #'@note Decision based on the Big Rings that Surround Each other
  message("# Pith Decision by Mutually Surrounded Clusters #")
  
  # Validate Position
  v.p <- which(order != 0)
  
  # db_scan clustering the center position
  df.valid <- df.check[v.p, c("x.c", "y.c")]
  # r.search <-
  #   df.check[v.p,] %>% 
  #   mutate(dx = xmax - xmin,
  #          dy = ymax - ymin,
  #          rd = min(dx, dy)) %>% 
  #   select(rd) %>% 
  #   unique()
  # eps_scale = r.search/2
  
  #'[Method: kNNdist]
  #'[Find better eps_scale]
  #'@note See: dbscan.eps
  eps_scale <- df.valid %>% dbscan.eps()
  frnn_dat  <- dbscan::frNN(df.valid, eps = eps_scale)
  db_dat    <- dbscan::dbscan(frnn_dat, minPts = 3)
  
  # Extract the clusters
  # in case double pith
  # 2 db_scan groups are available
  message("# Extract Inner-most Clusters as Pith #")
  
  n.clt <- unique(db_dat$cluster)
  n.c   <- n.clt[which(n.clt != 0)]
  n.cp  <- n.c %>% length()
  
  i = 1
  tar.p <- vector()
  while (TRUE) {
    
    ext.p <- which(db_dat$cluster == n.c[i])
    tem.p <- v.p[ext.p[which.min(order[v.p[ext.p]])]]
    tar.p <- append(tar.p, tem.p)
    
    i = i+1
    if(i > n.cp){break}
    
  }
  
  # Combine the result
  df.check$order <- order
  df.check$core  <- 0
  df.check$core[tar.p] <- 1
  
  # Report
  message("# Done #")
  
  # df.list <- list(df.check = df.check,
  #                 df.break = df.break)
  # return(df.list)
  df.res <- rbind(df.check, df.break)
  return(df.res)
  
}



## Get Pith position
#'@description
#' Sort and order the ring segments based on their relative positions.
#' Aim at targeting the most central segment cluster.
#'@import
#' (1) df: df.summary for "core" from pith.CenterCheck
#' (2) df.db: data frame with dbscan clustering
#'@export
#' (1) data frame with pith information
#' (2) CImage Object with Pith drawn (im.make = T, im required)
#'@source
pith.GetPith <- function(df, df.db, im.make = F, im = NULL){
  
  #'[Pith_Center:]
  # Get the Pith center:
  range.p = 
    df$clst[which(df$core == 1)] %>% 
    as.character() %>% 
    as.integer()
  
  # Extract positions:
  # Due to the nature of "which"
  # the input can only be one, (Not useful for more than 1 pith situation)
  # Loop the feed in process:
  i = 1
  n.pith = length(range.p)
  core.loc <- list()
  while(TRUE){
    
    loc.temp      <- df.db[which(df.db$clst == range.p[i]),]
    core.loc[[i]] <- loc.temp
    i = i+1
    
    # Stop the loop
    if(i > n.pith){break}
    
  }
  
  # Remove the useless initiating row from core.loc
  core.loc <- core.loc %>% bind_rows()
  
  if(im.make == T){
    
    if(purrr::is_empty(im) == TRUE){
      stop("The 'im.make' is Activate without a defined im Input")
    }
    
    if(is.matrix(im) != TRUE){
      stop("The Input im Requires a 'Matrix' Structure")
    }
    
    # Re-create the image:
    im.pith    <- im
    im.pith[,] <- 0
    im.pith[as.matrix(core.loc[,c("x", "y")])] <- 1
    im.pith <- im.pith %>% imager::as.cimg()
    
  }
  
  # Output Results:
  if(im.make == T){
    message("Image Made & List Results Returned")
    res <- list(df.pith = core.loc, im.pith = im.pith)
  }else{
    res <- core.loc
  }
  return(res)
  
}



## Extract Pith Position
#'@description
#' Fit circle-ellipse to the segment cluster of pith. 
#' Aiming at deriving pith center and radius.
#'@note
#' This was originally designed for getting multiple pith centers.
#' However, as the update method for pith.CenterCheck only valid for 
#' single pith, the processing pipe line will only derive one pith center.
#' The while statements are preserved just in case, in the future, 
#' there are better ways for capturing multiple pith.
#'@import
#' (1) df.getpith: input from GetPith data frame of 
#'@export
#' (1) data frame with pith center position and radius
#'@source
pith.ExtractPith <- function(df.getpith){
  
  # As the input will be a list by im.make = TRUE
  df.getpith <- df.pith
  if(is.list(df.getpith) == TRUE){
    df.getpith <- df.getpith$df.pith
  }
  
  #'[Circular_Fitting:]
  # For multiple Pith Situation
  # The output is a table with cluster ID, x, y, radius
  pith.plc <- df.getpith$clst %>% unique() %>% as.character() %>% as.integer()
  pith.num <- pith.plc %>% length()
  
  # Fit positions to circles
  #'@note 
  # Target: Get fitted circle center and its radius
  # Cite: https://stackoverflow.com/questions/27169122/
  #               how-to-find-a-best-fit-circle-ellipse-using-r
  # Author: Spacedman
  fitSS <- function(xy,
                    a0=mean(xy[,1]),
                    b0=mean(xy[,2]),
                    r0 = mean(sqrt((xy[,1]-a0)^2 + (xy[,2]-b0)^2)),
                    ...){
    SS <- function(abr){
      sum((abr[3] - sqrt((xy[,1]-abr[1])^2 + (xy[,2]-abr[2])^2))^2)
    }
    optim(c(a0,b0,r0), SS, ...)
  }
  
  # Process Inner Tree ring data from GetPith()
  pith.df <- list()
  i = 1
  while (TRUE) {
    
    temp.df  <- df.getpith[which(df.getpith$clst == pith.plc[i]), c("x", "y")]
    temp.fit <- fitSS(temp.df)
    temp.dat <- data.frame(clst = pith.plc[i],
                           x = temp.fit$par[1], 
                           y = temp.fit$par[2],
                           r = temp.fit$par[3])
    pith.df[[i]] <- temp.dat
    
    i = i+1
    if(i > pith.num){break}
    
  }
  
  pith.df <- pith.df %>% bind_rows()
  
  return(pith.df)
  
}



## Automate eps decision making for dbscan
#'@description
#' This process is suggested in the dbscan package.
#' This function is only for simplifying the usage.
#'@import
#' data frame with x & y
#'@export
#' A recommended eps_scale value for dbscan process.
#'@source
dbscan.eps <- function(data){
  
  # -------------------------------------------------------------------------- #
  #'[Structural_Check:]
  if(ncol(data) == 2){
    
    print("Input as X-Y position format")
    
  }else{
    
    print("Input as Multi-dimentional format")
    
  }
  
  print("k = dim + 1")
  k = ncol(data) + 1
  # -------------------------------------------------------------------------- #
  
  #  [Method 2. kNNdist]
  #' [Find better eps_scale]
  print("kNNdist Processing...")
  df.kNN       <-  dbscan::kNNdist(as.matrix(data), k = k) %>% sort()
  df.kNN_dist  <- df.kNN - (df.kNN %>% dplyr::lag())
  df.kNN_dist  <- df.kNN_dist[-1] # Remove NA
  
  # Create mean kNN distance vector
  # [Slow]
  # i  = 1
  # ni = length(df.kNN)
  # df.kNN_mdist <- vector()
  # while(TRUE){
  #   
  #   print(i)
  #   
  #   df.kNN_mdist[i] <- mean(df.kNN_dist[1:i]) # Slow
  #   
  #   i = i + 1
  #   
  #   if(i > ni){break}
  #   
  # }
  df.kNN_mdist <- cumsum(df.kNN_dist) / seq_along(df.kNN_dist)
  
  # Use "Scale" to check whether the distance is large enough
  print("kNNplot knee Search")
  .scale = 100 # Distance change should > 5 mean dist
  i  = 1
  ni = length(df.kNN_mdist)
  eps_select.n <- vector()
  while (purrr::is_empty(eps_select.n)) {
    
    eps_select.n <- which(df.kNN_dist > df.kNN_mdist * .scale)
    
    if(length(eps_select.n)>1){eps_select.n <- eps_select.n[1]}
    
    print(i)
    i = i + 1
    .scale <- .scale - 1
    
    if(.scale < 1){stop("Strange Situation")}
    
  }
  
  eps_scale <- 
    df.kNN[eps_select.n:(eps_select.n+1)] %>% 
    mean() %>% round(digits = 2)
  
  return(eps_scale)
  
}