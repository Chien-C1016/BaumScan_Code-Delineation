# ---------------------------------------------------------------------------- #
# Order Clusters by Angular position
# ---------------------------------------------------------------------------- #
#'@usage 
#' The Whole Code is too complicated to revise and review.
#' The Newer version will handle this issue by separate the code into
#' different sections.
#' And finally approach them with a handler function.
#' [Goal]
#' Each section should only do few and simple things.
#'@section 01
#' Angular check the existence and Angular distances
#'@concept 
#' This function is checking clusters by wished searching Angle.
#' Additionally, give the distance computations in a list format for successive
#' functions.
#'@import 
#' (1) Pith position: (xc, yc)
#' (2) dbscanned df with clst, theta, det.y (df.mt)
#' (3) theta.range controls how big the angular search should be.(in degree)
#' (4) clust.size controls how big the cluster size should be considered.
#'@export
#' (1) df.clst     : New df.mt
#' (2) summary     : A unique clusters in different Angular position
#' (3) visibility  : to show how are the clusters distribution
#' (4) clst.trust  : Trusted clusters for further process.
#' (5) theta.range : Record the input searching range (in degree)

AngularCompute <- function(xc, yc, df.mt, theta.range){
  
  require(dplyr)
  
  # Data
  df.clst <- df.mt[, c("x", "y", "clst", "theta")]
  xc
  yc
  
  # Dist computation and input ----
  edge.vec <- df.clst[,c("x", "y")] %>% as.matrix()
  pith.vec <- matrix(c(xc, yc), nrow = dim(edge.vec)[1], ncol = 2, byrow = T)
  dist     <- sqrt(rowSums((edge.vec - pith.vec)^2)) %>% as.numeric()
  df.clst$dist <- dist
  
  # -------------------------------------------------------------------------- #
  # Cluster Angular Frequency.V1 (3s) - [NOT RUN] ----
  #' clst.ct <- vector()
  #' i = 1
  #' theta.range
  #' n.theta = 360/theta.range
  #' while (TRUE) {
  #'   
  #'   # Theta limits
  #'   theta.upb <- (0+theta.range*i)/180*pi
  #'   theta.lwb <- (0+theta.range*(i-1))/180*pi
  #'   
  #'   # substract the data frame
  #'   df.temp <- df.clst[df.clst$theta >= theta.lwb,]
  #'   df.temp <- df.temp[df.temp$theta <  theta.upb,]
  #'   
  #'   clst.tp <- df.temp$clst %>% unique() # Temp
  #'   clst.ct <- append(clst.ct, clst.tp)  # theta cluster count
  #'   
  #'   # Move to next theta range
  #'   i = i+1
  #'   
  #'   # Stop condition
  #'   if(i > n.theta){break}
  #'   
  #' }
  #' 
  #' clst.summary <- clst.ct %>% table() %>% as.data.frame()
  #' names(clst.summary) <- c("clst", "Freq")
  #' 
  #' # Visibility check to determine the trusted clusters
  #' trust.freq <- n.theta
  #' if(max(clst.summary$Freq) < trust.freq){
  #'   trust.freq <- max(clst.summary$Freq)
  #' }
  #' clst.p     <- which(clst.summary$Freq >= trust.freq)
  #' clst.trust <- clst.p
  #' # clst.trust <- clst.summary$clst[clst.p] %>% as.character() %>% as.integer()
  #' #'@note 
  #' #' As the clst.summary is "tabled", 
  #' #' the sequence of the row is ordered as the cluster number.
  #' #' Therefore, trusted cluster = positions in the summary
  
  # -------------------------------------------------------------------------- #
  # Cluster Angular Frequency.V2 (<1s) ----
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
  clst.summary    <- df.clst %>% 
    group_by(clst, Ang.pos) %>%
    summarise(n.point = n(),
              mn.dist = mean(dist),
              .groups = "drop")
  
  clst.visibility <- clst.summary %>%
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
  
  # Results ----
  res.list <- list()
  res.list$df.clst     <- df.clst
  res.list$summary     <- clst.summary
  res.list$visibility  <- clst.visibility
  res.list$clst.trust  <- clst.trust
  res.list$theta.range <- theta.range
  return(res.list)
  
}