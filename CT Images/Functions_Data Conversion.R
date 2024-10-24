# ---------------------------------------------------------------------------- #
# Functions of Tree Ring Structure Convertion between Softeares #
# ---------------------------------------------------------------------------- #
# Convert from R to ImageJ system ----
## Function to reorder points within each cluster
reorder_points <- function(df) {
  
  # Create a TSP object from the x and y coordinates
  etsp <- ETSP(df[, c("x", "y")])
  
  # Solve the TSP to find an optimal path
  tour <- solve_TSP(etsp)
  
  # Reorder the dataframe based on the TSP solution
  return(df[tour, ])
}



## Convert ring structure Rdata to ImageJ data 
R2IJ <- function(ring.df,
                 im.matrix){
  
  require(dplyr)
  require(TSP)
  
  temp <-
    ring.MatrixBuild(im.ref  = im.matrix, 
                     ring.df = ring.df) %>% 
    t() 
  
  ring.id <- temp %>% as.vector() %>% unique()
  ring.id <- ring.id[ring.id != 0]
  
  i = 1
  ring.list <- list()
  while (TRUE) {
    
    .temp <- temp
    .temp[.temp[,]!=ring.id[i]] <- 0
    
    ring.list[[i]] <-
      .temp %>% 
      imager::as.cimg() %>% 
      imager::get.locations(., function(v) v>0) %>% 
      dplyr::select("x", "y") %>% 
      reorder_points() %>% 
      dplyr::mutate(clst = ring.id[i])
    
    i = i + 1
    if(i > length(ring.id)){break}
    
  }
  
  ring.df <- ring.list %>% dplyr::bind_rows()
  return(ring.df)
  
}



# Convert from ImageJ to R system ----
Ij2R <- function(manual.data, 
                 im.matrix = NULL,
                 im.rebuid = "linear"){
  
  # Rings to be converted sequentially
  Label.name <- manual.data$Label %>% unique()
  n.label = length(Label.name)
  
  # ---------------------------------------------------------------------------#
  rings.df <- list()
  i = 1
  while (TRUE) {
    
    temp.ID <- Label.name[i]
    temp.df <- manual.data[manual.data$Label == temp.ID, c("X", "Y")]
    names(temp.df) <- c("x", "y")
    rings.df[[i]]  <- temp.df %>% dplyr::distinct()
    i = i+1
    
    if(i>n.label){break}
    
  }
  
  # Add names by Labels from ImageJ
  names(rings.df) <- Label.name
  
  # ---------------------------------------------------------------------------#
  # When ImageJ is recording the ROI,
  # It normally just records the edge of the ROI.
  # To acquire the complete ring, we need to re-build them.
  
  if(purrr::is_empty(im.rebuid) != TRUE){
    
    message("# Re-build ROI from recorded edges by ", im.rebuid, " method #")
    print("Process Ring ID:")
    
    i = 1
    df.i <- list()
    while (i <= n.label) {
      
      print(i)
      temp.df  <- rings.df[[i]]
      temp.seq <- rbind(temp.df, temp.df[1,])
      delta.x  <- diff(temp.seq$x) %>% abs()
      delta.y  <- diff(temp.seq$y) %>% abs()
      
      j = 1
      df.j   <- list()
      temp.x <- vector()
      temp.y <- vector()
      # Any duplicates produced by "seq" will be cleaned-out by dplyr
      while (j <= nrow(temp.df)) {
        
        # Determine the seq.length
        delta.length <- max((delta.x[j]+1), (delta.y[j]+1))
        
        # Generate Sequence
        if(delta.x[j] > 1){
          
          temp.x <- seq(temp.seq[j,]$x,
                        temp.seq[j+1,]$x,
                        length.out = delta.length)
          
          if(delta.y[j] > 1){
            
            temp.y <- 
              seq(temp.seq[j,]$y,
                  temp.seq[j+1,]$y,
                  length.out = delta.length) %>% 
              round()
            
          }else{
            
            temp.y <- rep.int(temp.seq[j,]$y, delta.length)
            
          }
          
        }else{ # delta.x[j] <= 1
          
          if(delta.y[j] > 1){
            
            temp.y <- 
              seq(temp.seq[j,]$y,
                  temp.seq[j+1,]$y,
                  length.out = delta.length) %>% 
              round()
            temp.x <- rep.int(temp.seq[j,]$x, delta.length)
            
          }else{
            
            # delta.x & delta.y <= 1 do not need to be extended.
            temp.x <- temp.seq[j,]$x
            temp.y <- temp.seq[j,]$y
            
          }
          
        }
        
        # j-loop control
        df.j[[j]] <- data.frame(x = temp.x, y = temp.y)
        j = j + 1
        
      }
      
      # i-loop control
      df.i[[i]] <- df.j %>% dplyr::bind_rows() %>% dplyr::distinct()
      i = i + 1
      
    }
    
    rings.df <- df.i
    
    # Add names by Labels from ImageJ
    names(rings.df) <- Label.name
    
  }
  
  # Save Image Position
  df.ROI <- rings.df
  
  # ---------------------------------------------------------------------------#
  # When Image is read by R,
  # The row and Column switched as transverse matrix.
  # To cope with the coordinates from ImageJ, t() is needed
  
  rings.im <- list()
  rings.R  <- list()
  if (purrr::is_empty(im.matrix) != TRUE) {
    
    message("# Convert positions into R environment")
    
    # Empty t() matrix image preparation
    temp.empty    <- im.matrix
    temp.empty[,] <- 0
    temp.empty    <- temp.empty %>% t()
    
    i = 1
    while (TRUE) {
      
      temp.im             <- temp.empty
      temp.check          <- rings.df[[i]] %>% as.matrix()
      temp.im[temp.check] <- 1
      rings.im[[i]]       <- temp.im %>% t() %>% imager::as.cimg()
      
      # Convert into R format
      rings.R.i <- 
        rings.im[[i]] %>% 
        imager::get.locations(., function(v) v > 0) %>% 
        dplyr::mutate(clst = i %>% as.character()) %>% 
        dplyr::select(c("x", "y", "clst"))
      rings.R[[i]] <- rings.R.i
      
      i = i+1
      
      if(i>n.label){break}
      
    }
    
    # Add names by Labels from ImageJ
    names(rings.im) <- Label.name
    names(rings.R)  <- Label.name
    
  }else{
    
    message("# Computed rings.df is NOT proper for R, t() is needed")
    
  }
  
  # ---------------------------------------------------------------------------#
  
  rings <- list(rings.df = rings.df,
                rings.im = rings.im,
                rings.R  = rings.R)
  
  # Return
  return(rings)
  
}