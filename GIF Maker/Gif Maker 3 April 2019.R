#### Packages #####
library(fs)
library(magick)
library(tidyverse)
###################

#### Directories ####
output_folder <- "/Users/aita2652/Cell Profiler Material/LiveCell_2May2019/Python/GIFs"
#^ Directory where you would like to store your GIFs

parent_wd <- "/Users/aita2652/Cell Profiler Material/LiveCell_2May2019/Python/U2OS-20X-BpAHDDB-60uM30uM15uM-LIVE"
#^ Directory where you Python code's image output are stored. Note: The level is important!. Should be right 
#above the folders which correspond to the wells.

setwd(parent_wd) #Sets parent directory as the working directory.

######################

#### User Defined Functions ####
gif_maker <- function(y)
{
  list_png <- list.files(recursive = FALSE, pattern = "*.tif")
  #to make a list of filenames only! Not their paths. It is to name the images in gif and let map() walk.
  
  a <- getwd()
  gif_smartname <- paste0(basename(dirname(a)),basename(a))
  
  list_png %>%
    
    map(function(fnme){
      image_read(fnme)%>%
        image_annotate(fnme, color = "white", gravity = "south", size = 100)
    })%>%
    
    image_join() %>%  # joins image
    image_annotate("Time Point:", size = 70, color = "white", gravity = "southwest") %>%
    image_annotate(gif_smartname, size = 50, color = "white", gravity = "north") %>%
    image_scale("800") %>% 
    image_animate(fps = 1) %>%
    image_write(path = file.path(y,(paste0(gif_smartname,".gif"))))
}

#################################

#### Main Body ####

newdir <- list.dirs(list.dirs(getwd(), recursive = F), recursive = F) 
#^ Stores the paths of all the directories and sub-directories

for (i in 1:length(newdir))
{
  setwd(newdir[i])
  gif_maker(output_folder)
}

setwd(parent_wd)

print("All GIFs made!")
#################################

#####Credits:
###https://community.rstudio.com/t/how-to-import-a-lot-of-images-from-a-folder/19240/8
###https://stackoverflow.com/questions/22069095/r-get-list-of-files-but-not-of-directories
###https://stackoverflow.com/questions/46923970/r-dplyr-read-list-of-files-and-use-filename-as-a-variable

