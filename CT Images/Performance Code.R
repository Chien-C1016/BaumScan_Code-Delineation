
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

res <- 0.15 # mm

# ---------------------------------------------------------------------------- #
#'[Read_Data:]
#'*Sequence Read_Data: Fi-Ki-La*
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

TRA_data <- lapply(list_files, read.csv)

#'[Tree_Ring_Manual]
# Specify the directory path
data_path <- paste0(Dir, "/Manual/Tree_Ring")
# Tree_Ring Data
list_files <- list.files(path = data_path, 
                         pattern = "^[A-Z][a-z]_0[1-5]_Results\\.csv$",
                         full.names = TRUE)
print(list_files)

TRM_data <- lapply(list_files, read.csv)

#'[Sapwood_Algorithm]
# Specify the directory path
data_path <- paste0(Dir, "/Algorithm/Sapwood")
# Sapwood Data
list_files <- list.files(path = data_path, 
                         pattern = "^[A-Z][a-z]_0[1-5]_Results\\.csv$",
                         full.names = TRUE)
print(list_files)

SWA_data <- lapply(list_files, read.csv)

#'[Sapwood_Manual]
# Specify the directory path
data_path <- paste0(Dir, "/Manual/Sapwood")
# Sapwood Data
list_files <- list.files(path = data_path, 
                         pattern = "^[A-Z][a-z]_0[1-5]_Results\\.csv$",
                         full.names = TRUE)
print(list_files)

SWM_data <- lapply(list_files, read.csv)

#'[IMG_Moisture-Preserved Stem Discs]
# Specify the directory path
data_path <- paste0(Dir, "/IMG")
# Sapwood Data
list_files <- list.files(path = data_path, 
                         pattern = "^[A-Z][a-z]-0[1-5]\\.tif$",
                         full.names = TRUE)
print(list_files)

IMG_data <- lapply(list_files, OpenImageR::readImage)

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
  ggplot(., aes(x = dx, y = dy, color = seq_disc)) +
  geom_point(size = 3, alpha = 0.7) +  # Semi-transparent points for visibility
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
































