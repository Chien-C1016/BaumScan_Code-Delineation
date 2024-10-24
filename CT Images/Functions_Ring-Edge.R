# ---------------------------------------------------------------------------- #
# Functions  #
# ---------------------------------------------------------------------------- #
# General ----
#'@description
#' 

# ---------------------------------------------------------------------------- #
# Functions ----





edge.vectorclean <- function(xc,
                             yc,
                             df.gradient,
                             cos.sigma,
                             mag.thres,
                             mag.method = c("quantile"),
                             as.mask = TRUE,
                             im.ref = NULL,
                             .plot = FALSE){
  
  xc
  yc
  if(purrr::is_empty(im.ref)){stop("Reference Image is Missing")}
  
  pith.vec <- df.gr %>% mutate(.,
                               .keep = c("none"),
                               pi.dx = xc - x,
                               pi.dy = yc - y)
  edge.vec <- df.gr[,c(4,5)]
  edge.cos <- rowSums(pith.vec*edge.vec) / 
    ((sqrt(rowSums(pith.vec^2)))*(sqrt(rowSums(edge.vec^2))))
  head(edge.cos)
  
  # ' - Angular Thresholding: ----
  # acos(0.5)*180/pi # 60 degrees
  hist(edge.cos)
  quantile(edge.cos, na.rm = T)
  sigma = cos.sigma
  df.edge.clean <- df.gr[which(edge.cos > sigma), c(1:3)]
  head(df.edge.clean)
  
  # ' - mag Thresholding: ----
  # Pre-check about the distribution of the magnitudes
  hist(df.edge.clean$mag)
  quantile(df.edge.clean$mag)
  if(mag.method == "quantile"){
    mag.thres = quantile(df.edge.clean$mag, mag.thres)
  }else{
    mag.thres = mag.thres
  }
  df.edge.clean.mt <- df.edge.clean[df.edge.clean$mag > mag.thres,]
  
  # Back conversion to image
  # [Below slow NOT RUN]
  # df.edge.mt <- merge(df.cim, df.edge.clean.mt, all.x = T)
  # df.edge.mt <- data.frame(x = df.edge.mt$x,
  #                          y = df.edge.mt$y,
  #                          value = 1 - as.numeric(is.na(df.edge.mt$mag)))
  # im.edge.mt <- df.edge.mt %>% as.cimg(., dims = dim(cim))
  
  if(imager::is.cimg(im.ref)){
    
    im.ref <- im.ref %>% as.matrix()
    
  }
  
  if(as.mask == TRUE){
    
    edge.matrix <- im.ref
    edge.matrix[,] <- 0
    edge.matrix[as.matrix(df.edge.clean.mt[,c("x", "y")])] <- 1
    im.edge.mt <- edge.matrix %>% as.cimg()
    
    # Plot check & Visualization
    if(.plot == TRUE){
      plot(cim)
      highlight((im.edge.mt > 0))
    }
    
  }else{
    
    im.edge.mt <- 
      df.edge.clean.mt %>% 
      imager::as.cimg(., v.name = "mag", dims = dim(im.ref))
    
  }
  
  return(im.edge.mt)
  
}



edge.vonMask <- function(mask, output = "df"){
  
  require(dplyr)
  
  # Mask data format
  if(imager::is.cimg(mask)){
    
    mask <- mask
    
  }else if(is.matrix(mask)){
    
    mask <- mask %>% imager::as.cimg()
    
  }else if(imager::is.pixset(mask)){
    
    mask <- mask %>% imager::as.cimg()
    
  }
  
  contour.df <-
    mask %>% 
    imager::contours(., nlevels = 1) %>% 
    as.data.frame()
  
  # Output
  if(output == "img" | output == "All"){
    
    contour.img <- 
      contour.df %>% 
      imager::as.cimg(., v.name = "level", dims = dim(mask))
    
    if(output == "All"){
      
      res.list <- list(contour.df  = contour.df,
                       contour.img = contour.img)
      return(res.list)
      
    }else{
      
      return(contour.img)
      
    }
    
  }else{
    
    if(output != "df"){
      
      message("# No such output method defined, return data frame instead #")
    }
    return(contour.df)
    
  }
  
}


edge.FeatureCompute <- function(xc, yc, df.edge){
  
  # Computation:
  vec.edge.mt <- df.edge %>% mutate(.,
                                    .keep = c('none'),
                                    dx = .$x - xc,
                                    dy = .$y - yc)
  vec.pith.x  <- matrix(c(1,0), 
                        nrow = dim(vec.edge.mt)[1], 
                        ncol = 2, 
                        byrow = T)
  
  # Determinants "y"
  det.y <- (vec.edge.mt[,2] < 0) %>% as.numeric()
  
  # Angle calculation
  cos.theta <- rowSums(vec.edge.mt*vec.pith.x)/(sqrt(rowSums(vec.edge.mt^2))*1)
  theta.2pi <- abs(det.y*2*pi - acos(cos.theta))
  
  # data frame illustrates:
  df.mt <- df.edge %>% dplyr::mutate(theta = theta.2pi)
  
  # Dist computation and input
  edge.vec <- df.mt[,c("x", "y")] %>% as.matrix()
  pith.vec <- matrix(c(xc, yc), nrow = dim(edge.vec)[1], ncol = 2, byrow = T)
  dist     <- sqrt(rowSums((edge.vec - pith.vec)^2)) %>% as.numeric()
  df.mt$dist <- dist
  
  return(df.mt)
  
}


edge.DistResample <- function(df.FeatureCompute,
                              resample.digits = 2){
  
  df.DistResample <-
    df.FeatureCompute %>% 
    dplyr::mutate(theta = 
                    theta %>% 
                    round(., digits = resample.digits) %>% 
                    as.character() %>% 
                    as.numeric()) %>% 
    dplyr::group_by(clst, theta) %>% 
    dplyr::summarise(.pos = which.max(dist),
                     x = x[.pos],
                     y = y[.pos],
                     dist = dist[.pos]) %>% 
    dplyr::select(c("clst", "x", "y", "theta", "dist"))
  
  return(df.DistResample)
  
}


edge.MagnitudeClean <- function(df.edge,
                                ReSample.Size = 0.0001){
  
  # Secure the data structure
  ReSample.Size = log10(ReSample.Size) %>% abs
  df.edge <- df.edge[, c("x", "y", "value", "clst", "theta")]
  
  df.edge <- 
    df.edge %>% 
    dplyr::mutate(theta = 
                    theta %>% 
                    round(., digits = ReSample.Size) %>% 
                    as.character() %>% 
                    as.numeric()) %>% 
    dplyr::group_by(clst, theta) %>% 
    dplyr::summarise(.pos = which.max(value),
                     x = x[.pos],
                     y = y[.pos],
                     value = max(value)) %>% 
    dplyr::select(c("clst", "theta", "x", "y", "value"))
  
  return(df.edge)
  
}


edge.DistGrouping <- function(df.edge,
                              group.dist = 2,
                              edge.size.thres = 30){
  
  require(dplyr)
  require(purrr)
  
  if(purrr::is_empty(edge.size.thres)){
    
    message("ALL Edge-Segments are Considered")
    message("[Note] TIME Consuming ...")
    
    # Control
    clst.ID <- df.edge$clst %>% unique() %>% as.integer()
    
  }else{
    
    # Signal
    size.sep = TRUE
    
    df.edge.sum <- 
      df.edge$clst %>% 
      table() %>% 
      dplyr::as_tibble()
    
    # Cluster Separation
    clst.select <- 
      df.edge.sum %>% 
      dplyr::filter(n > edge.size.thres) %>% 
      dplyr::select(".") %>% 
      unlist() %>% 
      unname() %>% 
      as.integer()
    
    # Control
    clst.ID <- clst.select
    
  }
  
  # Save df.edge
  .df.edge <- df.edge
  
  # Clust.ID
  #' @description 
  #' This is the conotrol of which clusters to be checked.
  clst.ID # <- df.edge$clst %>% unique() %>% as.integer()
  
  i = 1
  clst.i <- list()
  clst.group <- df.edge$clst
  while (TRUE) {
    
    # message("# i = ", i)
    print(paste("i =", i))
    
    ta.clst <- df.edge %>% dplyr::filter(clst == clst.ID[i])
    df.clst <- df.edge %>% dplyr::filter(clst != clst.ID[i])
    
    # V1
    # j = 1
    # group.list.j = list()
    # temp.df.clst <- df.clst
    # temp.df      <- df.clst[, c("x", "y")] %>% as.matrix()
    # temp.clst    <- df.clst$clst
    # while (TRUE) {
    # 
    #   temp.ta <- matrix(data = ta.clst[j, c("x", "y")] %>% as.matrix(),
    #                     nrow = nrow(temp.df),
    #                     ncol = 2,
    #                     byrow = TRUE)
    # 
    #   temp.dist     <- sqrt(rowSums((temp.df - temp.ta)^2))
    #   temp.dist.sum <-
    #     dplyr::tibble(clst = temp.clst,
    #                   dist = temp.dist) %>%
    #     dplyr::group_by(clst) %>%
    #     dplyr::summarise(dist.min = min(dist),
    #                      npts     = n()) %>%
    #     dplyr::filter(dist.min < group.dist)
    # 
    #   # j-loop updates
    #   print(j)
    #   group.list.j[[j]] <- temp.dist.sum
    #   temp.clst.ID <- dplyr::setdiff(df.clst$clst, temp.dist.sum$clst)
    #   temp.df.clst <- temp.df.clst %>% dplyr::filter(clst %in% temp.clst.ID)
    #   temp.df      <- temp.df.clst[, c("x", "y")] %>% as.matrix()
    #   temp.clst    <- temp.df.clst$clst
    #   j = j + 1
    # 
    #   if(j > nrow(ta.clst)){break}
    # 
    # }
    
    # V2
    # j = 1
    # clst.j  <- vector()
    # dist.j  <- vector()
    # temp.df <- df.clst[, c("x", "y")] %>% as.matrix()
    # while (TRUE) {
    # 
    #   temp.ta <- matrix(data = ta.clst[j, c("x", "y")] %>% as.matrix(),
    #                     nrow = nrow(temp.df),
    #                     ncol = 2,
    #                     byrow = TRUE)
    # 
    #   temp.dist     <- sqrt(rowSums((temp.df - temp.ta)^2))
    #   temp.dist.sum <-
    #     dplyr::tibble(clst = df.clst$clst,
    #                   dist = temp.dist) %>%
    #     dplyr::group_by(clst) %>%
    #     dplyr::summarise(dist.min = min(dist),
    #                      npts     = n()) %>%
    #     dplyr::filter(dist.min < group.dist)
    # 
    #   # j-loop updates
    #   print(j)
    #   clst.j <- append(clst.j, temp.dist.sum$clst)
    #   dist.j <- append(dist.j, temp.dist.sum$dist.min)
    #   j = j + 1
    # 
    #   if(j > nrow(ta.clst)){break}
    # 
    # }
    
    # V3
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
    
    clst.i[[i]] <- 
      data.frame(clst = as.numeric(levels(df.i$.))[df.i$.],
                 dist = sum.i / df.i$Freq)
    
    # New Group Cluster ID
    .pos <- which(clst.group %in% clst.j)
    clst.group[.pos] <- clst.ID[i]
    
    # PAULSE
    # if(length(clst.j %>% unique()) > 1){break}
    
    # i-loop Control
    i = i + 1
    if(i > length(clst.ID)){break}
    
  }
  
  clst.i %>% bind_rows(., .id = "clst.ID") %>% as_tibble()
  
  # Plot Check
  plot(df.edge$x, 
       df.edge$y,
       cex = 0.3,
       col = as.factor(clst.group))
  
  df.group      <- df.edge
  df.group$clst <- clst.group
  
  return(df.group)
  
}



edge.filter <- function(df.clst, min.npts = 1){
  
  #'[Summary of df.clst]
  df.sum <-
    df.clst %>% 
    dplyr::group_by(clst) %>% 
    dplyr::summarise(n = n())
  
  #'[Filtering df.clst by min.npts]
  clst.filtered <- df.sum$clst[which(df.sum$n != min.npts)]
  df.clst <- df.clst %>% dplyr::filter(clst %in% clst.filtered)
  
  # Output
  return(df.clst)
  
}

