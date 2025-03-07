
# ---------------------------------------------------------------------------- #
# Environment ----
#'[Dir_n_Funs:]
library(dplyr)
library(purrr)
library(ggplot2)

Dir <- "~/Data"
setwd ("~/Data")

source("~/Data/Scripts/CT Images/Functions_Pith.R")
source("~/Data/Scripts/CT Images/Functions_Ring-Clust.R")
source("~/Data/Scripts/CT Images/Functions_Ring-Edge.R")
source("~/Data/Scripts/CT Images/Ring Distance Computation.R")
source("~/Data/Scripts/CT Images/Functions_Data Conversion.R")


res <- 0.15 # mm

# ---------------------------------------------------------------------------- #
#'[Read_Data:]
#'*Sequence Read_Data: Fi-Ki-La*
#'[IMG_Moisture-Preserved Stem Discs]
# Specify the directory path
data_path <- paste0(Dir, "/IMG")
# Sapwood Data
list_files <- list.files(path = data_path, 
                         pattern = "^[A-Z][a-z]-0[1-5]\\.tif$",
                         full.names = TRUE)
print(list_files)

IMG_data <- lapply(list_files, OpenImageR::readImage)

#'[Pith_Algorithm]
# Specify the directory path
data_path <- paste0(Dir, "/Algorithm/Pith")
# Pith Data
list_files <- list.files(path = data_path, 
                         pattern = "^[A-Z][a-z]_0[1-5]_Pith\\.RData$",
                         full.names = TRUE)
print(list_files)

PA_data <- lapply(list_files, readRDS)

#'[Pith_Manual]
# Specify the directory path
data_path <- paste0(Dir, "/Manual/Pith")
# Pith Data
list_files <- list.files(path = data_path, 
                        pattern = "^Pith_[A-Z][a-z]_0[1-5]\\.tif$",
                        full.names = TRUE)
print(list_files)

PM_IMG <- lapply(list_files, OpenImageR::readImage)

PM_data <-
  PM_IMG %>%
  purrr::map(~ .x %>% 
               imager::as.cimg() %>%
               mask.as.Ring() %>% 
               pith.ExtractPith
  )

#'[Tree_Ring_Algorithm]
# Specify the directory path
data_path <- paste0(Dir, "/Algorithm/Tree_Ring")
# Tree_Ring Data
list_files <- list.files(path = data_path, 
                        pattern = "^[A-Z][a-z]_0[1-5]_Results_Auto\\.csv$",
                        full.names = TRUE)
print(list_files)

TRA_data <- purrr::map2(list_files, IMG_data, ~.x %>% 
                          read.csv() %>% 
                          Ij2R(., im.matrix = .y) %>% 
                          purrr::pluck(., "rings.R"))

#'[Tree_Ring_Manual]
# Specify the directory path
data_path <- paste0(Dir, "/Manual/Tree_Ring")
# Tree_Ring Data
list_files <- list.files(path = data_path, 
                         pattern = "^[A-Z][a-z]_0[1-5]_Results\\.csv$",
                         full.names = TRUE)
print(list_files)

TRM_data <- purrr::map2(list_files, IMG_data, ~.x %>% 
                          read.csv() %>% 
                          Ij2R(., im.matrix = .y) %>% 
                          purrr::pluck(., "rings.R"))

#'[Sapwood_Algorithm]
# Specify the directory path
data_path <- paste0(Dir, "/Algorithm/Sapwood")
# Sapwood Data
list_files <- list.files(path = data_path, 
                         pattern = "^[A-Z][a-z]_0[1-5]_Results\\.csv$",
                         full.names = TRUE)
print(list_files)

SWA_data <- purrr::map2(list_files, IMG_data, ~.x %>% 
                          read.csv() %>% 
                          Ij2R(., im.matrix = .y) %>% 
                          purrr::pluck(., "rings.R"))

#'[Sapwood_Manual]
# Specify the directory path
data_path <- paste0(Dir, "/Manual/Sapwood")
# Sapwood Data
list_files <- list.files(path = data_path, 
                         pattern = "^[A-Z][a-z]_0[1-5]_Results\\.csv$",
                         full.names = TRUE)
print(list_files)

SWM_data <- purrr::map2(list_files, IMG_data, ~.x %>% 
                          read.csv() %>% 
                          Ij2R(., im.matrix = .y) %>% 
                          purrr::pluck(., "rings.R"))

# ---------------------------------------------------------------------------- #
#'[Sequence_DiscID:]
seq_disc <- 
  rep("Spruce", 5) %>% 
  append(., rep("Pine",  5)) %>% 
  append(., rep("Larch", 5))

num_disc <- 
  c(1:5) %>% 
  append(., c(1:5)) %>% 
  append(., c(1:5))

# ---------------------------------------------------------------------------- #
#'[Species_Code:]
Species.ID <- c("All", "Spruce", "Pine", "Larch")
Color.code <- c("dodgerblue4", "tomato3", "aquamarine4", "goldenrod1")
Species.col <- setNames(Color.code, Species.ID)

# ---------------------------------------------------------------------------- #
# Performance_Pith ----
#'[Performance_Pith:]
performance_Pith <-
  purrr::map2(PA_data, PM_data, 
              ~ data.frame(
                dx = .x$x - .y$x,
                dy = .x$y - .y$y) %>% 
                dplyr::mutate(dist = sqrt(dx^2 + dy^2))
  ) %>% 
  dplyr::bind_rows() %>% 
  cbind(., seq_disc) %>% 
  cbind(., num_disc) %>% 
  dplyr::mutate(dx = dx * res,
                dy = dy * res,
                dist = dist * res)

saveRDS(performance_Pith, file = "performance_Pith.RData")

# Descriptive stats
mean_x <- performance_Pith$dx %>% mean()
mean_y <- performance_Pith$dy %>% mean()
mean_d <- performance_Pith$dist %>% mean()
sd_d   <- performance_Pith$dist %>% sd()
sum_species <-
  performance_Pith %>% 
  dplyr::group_by(seq_disc) %>% 
  summarise(mean_x = mean(dx),
            mean_y = mean(dy),
            mean_d = mean(dist),
            sd_x = sd(dx),
            sd_y = sd(dy),
            sd_d = sd(dist))

# Figure
performance_Pith %>% 
  dplyr::mutate(seq_disc = factor(seq_disc, levels = Species.ID)) %>% 
  ggplot(., aes(x = dx, y = dy, color = seq_disc)) +
  # Semi-transparent points for visibility
  geom_point(size = 3, alpha = 0.7) + 
  # Self-define color-code
  scale_color_manual(values = Species.col) +
  theme_minimal() +  # Minimalistic theme
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = mean_y, linetype = "dotdash", color = "tomato3") +
  geom_vline(xintercept = mean_x, linetype = "dotdash", color = "tomato3") +
  labs(x     = "Difference in X (mm)", 
       y     = "Difference in Y (mm)", 
       color = "Species") +
  coord_fixed(ratio = 1) +  # Ensures equal scaling of x and y axes+
  theme(legend.position = "right",
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.text  = element_text(size = 12),  
        axis.title = element_text(size = 12, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 

# ---------------------------------------------------------------------------- #
# Performance_Tree_Rings ----
#'[Performance_Tree_Rings:]
# Calculating the Accuracy of each Tree Rings
Disc.ID = 1
n.disc  = length(TRA_data)
performance_TR <- list()
while (TRUE) {
  
  data.list <- list(TRA_data[[Disc.ID]], 
                    TRM_data[[Disc.ID]], 
                    seq_along(TRA_data[[Disc.ID]]))
  
  performance_TR[[Disc.ID]] <-
    purrr::pmap(data.list,
                ~ RingDistance(Target.ring = ..1, 
                               Next.ring   = ..2, 
                               center = c(PM_data[[Disc.ID]]$x, 
                                          PM_data[[Disc.ID]]$y)) %>% 
                  dplyr::mutate(Ring.ID = ..3)) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(Species = seq_disc[Disc.ID],
                  Disc.ID = num_disc[Disc.ID],
                  Ring.dist = Ring.dist * res)
  
  # loop-control
  Disc.ID = Disc.ID + 1
  if(Disc.ID > n.disc){break}
  
} 

saveRDS(performance_TR, file = "performance_TR.RData")

## Figure ----
df_PTR <-
  performance_TR %>% 
  purrr::map(., ~ .x %>% 
               dplyr::mutate(clst = Ring.ID) %>% 
               ring.RelativeDist(ring.TRW = .,
                                 Target.ID = 1,
                                 .keep = TRUE) %>% 
               dplyr::select(-clst)) %>% 
  dplyr::bind_rows()

.df <- df_PTR %>% dplyr::mutate(Species = "All")
.df_PTR <- 
  rbind(.df, df_PTR) %>% 
  dplyr::mutate(Species = factor(Species, levels = Species.ID)) 

# Overall
Error_TRW <- 
  .df_PTR %>% 
  ggplot(.) +
  
  # Density Plot
  ggdensity::geom_hdr(aes(x = relative.position, 
                          y = Ring.dist,
                          fill = Species),
                      probs = c(0.85, 0.70, 0.5, 0.20, 0.10),
                      show.legend = c(alpha = TRUE)) +
  
  # Add horizontal line at the mean of y for each facet
  geom_hline(data = .df_PTR %>%
               dplyr::group_by(Species) %>%
               dplyr::summarise(y_mean = mean(Ring.dist)), 
             aes(yintercept = y_mean), 
             color = "tomato3", linetype = "dashed") +
  
  # Self-define color-code
  scale_fill_manual(values = Species.col) + 
  
  # Fix Y-Axis from 0 to 3
  scale_y_continuous(limits = c(0, 3)) +  
  
  theme_minimal() + 
  labs(x = "Relative Position to Pith", 
       y = "Absolute Error (mm)") +
  theme(legend.position = "right",
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        #panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.text  = element_text(size = 12),  
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_text(margin = margin(r = 5)),
        
        strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(2, "lines")) +
  
  facet_wrap(~Species, ncol = 2)

Error_TRW

# ---------------------------------------------------------------------------- #
# Performance_Sapwood ----
#'[Performance_Sapwood:]
# Calculating the Accuracy of each Sapwood Boundaries
Disc.ID = 1
n.disc  = length(SWA_data)
performance_SW <- list()
while (TRUE) {
  
  data.list <- list(SWA_data[[Disc.ID]], 
                    SWM_data[[Disc.ID]], 
                    seq_along(SWA_data[[Disc.ID]]))
  
  performance_SW[[Disc.ID]] <-
    purrr::pmap(data.list,
                ~ RingDistance(Target.ring = ..1, 
                               Next.ring   = ..2, 
                               center = c(PM_data[[Disc.ID]]$x, 
                                          PM_data[[Disc.ID]]$y)) %>% 
                  dplyr::mutate(Ring.ID = ..3)) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(Species = seq_disc[Disc.ID],
                  Disc.ID = num_disc[Disc.ID],
                  Ring.dist = Ring.dist * res)
  
  # loop-control
  Disc.ID = Disc.ID + 1
  if(Disc.ID > n.disc){break}
  
} 

saveRDS(performance_SW, file = "performance_SW.RData")

## Figure ----
df_PSW <-
  performance_SW %>% 
  purrr::map(., ~ .x %>% 
               dplyr::mutate(clst = Ring.ID) %>% 
               ring.RelativeDist(ring.TRW = .,
                                 Target.ID = 2,
                                 .keep = TRUE) %>% 
               dplyr::select(-clst)) %>% 
  dplyr::bind_rows()

## ' Case-01: Pool "All" data ----
#'[Size:751*908]
.df <- df_PSW %>% dplyr::mutate(Species = "All")
.df_PSW <- 
  rbind(.df, df_PSW) %>% 
  dplyr::mutate(Species = factor(Species, levels = Species.ID)) 

# Overall
Error_SW <- 
  .df_PSW %>% 
  ggplot(.) +
  
  # Density Plot
  ggdensity::geom_hdr(aes(x = relative.position, 
                          y = Ring.dist,
                          fill = Species),
                      probs = c(0.85, 0.70, 0.5, 0.20, 0.10),
                      show.legend = c(alpha = TRUE)) +
  
  # Add horizontal line at the mean of y for each facet
  geom_hline(data = .df_PSW %>%
               dplyr::group_by(Species) %>%
               dplyr::summarise(y_mean = mean(Ring.dist)), 
             aes(yintercept = y_mean), 
             color = "tomato3", linetype = "dashed") +
  
  # Self-define color-code
  scale_fill_manual(values = Species.col) + 
  
  # Fix Y-Axis from 0 to 3
  # scale_y_continuous(limits = c(0, 3)) +  
  
  theme_minimal() + 
  labs(x = "Relative Position to Pith", 
       y = "Absolute Error (mm)") +
  theme(legend.position = "right",
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        #panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.text  = element_text(size = 12),  
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_text(margin = margin(r = 5)),
        
        strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(2, "lines")) +
  
  facet_wrap(~Species, ncol = 2)

Error_SW

## ' Case-02: Only display Species data ----
#'[Size:1075*465]
.df_PSW <- 
  df_PSW %>% 
  dplyr::mutate(Species = factor(Species, levels = Species.ID))

# Overall
Error_SW <- 
  .df_PSW %>% 
  ggplot(.) +
  
  # Density Plot
  ggdensity::geom_hdr(aes(x = relative.position, 
                          y = Ring.dist,
                          fill = Species),
                      probs = c(0.85, 0.70, 0.5, 0.20, 0.10),
                      show.legend = c(alpha = TRUE)) +
  
  # Add horizontal line at the mean of y for each facet
  geom_hline(data = .df_PSW %>%
               dplyr::group_by(Species) %>%
               dplyr::summarise(y_mean = mean(Ring.dist)), 
             aes(yintercept = y_mean), 
             color = "tomato3", linetype = "dashed") +
  
  # Self-define color-code
  scale_fill_manual(values = Species.col) + 
  
  # Fix Y-Axis from 0 to 3
  # scale_y_continuous(limits = c(0, 3)) +  
  
  theme_minimal() + 
  labs(x = "Relative Position to Pith", 
       y = "Absolute Error (mm)") +
  theme(legend.position = "right",
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        #panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.text  = element_text(size = 12),  
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_text(margin = margin(r = 5)),
        
        strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(2, "lines")) +
  
  facet_wrap(~Species, ncol = 3) +
  
  # Manually reorder legend: Species first, Probs second
  guides(fill  = guide_legend(order = 1),  # Species first
         alpha = guide_legend(order = 2))  # Probs second

Error_SW






















