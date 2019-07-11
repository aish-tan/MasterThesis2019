#Aishvarya Thesis Main R Pipeline


##Installing necessary packages 

# install.packages("devtools")
# install.packages("fs")
# install.packages("tidyverse")
# devtools::install_github("kassambara/factoextra")
# install.packages("corrplot")
# install.packages('gridExtra')
# install.packages('directlabels')

##Reading packages, run everytime at initialization
library(tidyverse)
library(readr)
#library(ggbiplot) #devtools::install_github("vqv/ggbiplot")
library(stats)
library(tidyr)
library(dplyr)
library(factoextra) #devtools::install_github("kassambara/factoextra")
library(fs)
library(ggplot2)
library(corrplot) 
library(gridExtra)
library(directlabels) 
library(lattice) 


##########################
##Loading Directories

#Input directory, where the CellProfiler outputs are stored (USER PROVIDED)
data_dir <- "/Users/aita2652/Cell Profiler Material/15MayLiveCell/New CP Pipeline Output"

#Output directory, where the plots and spreadsheets are to be stored (USER PROVIDED)
out_dir <- "/Users/aita2652/Cell Profiler Material/15MayLiveCell/ROut11Jul"

#Stores the paths of all the directories and sub-directories
newdir <- list.dirs(data_dir, recursive = F)

##########################



##Basic Code: Does PCA and also lays the foundation for the rest of the code
#Might take upto a couple of minutes
for (i in 1:length(newdir))
{
  
  f1 <- paste(newdir[i], "Nucleus.csv", sep = "/")
  f2 <- paste(newdir[i], "Cytoplasm.csv", sep = "/")
  
  raw_nuc_data <- read.csv(file = f1, header = T, stringsAsFactors = F)
  raw_cyt_data <- read.csv(file = f2, header = T, stringsAsFactors = F)
  
  ##Renaming columns to avoid name clashes, at nucleus data
  raw_nuc_data <- rename_at(raw_nuc_data, vars(contains("AreaShape")), funs(paste0("Nuc_",.)))
  
  ##Renaming columns to avoid name clashes, here at cytoplasm data
  raw_cyt_data <- rename_at(raw_cyt_data, vars(contains("AreaShape")), funs(paste0("Mito_",.)))
  
  
  ##Joining nuclear readings and cytoplasm readings together
  raw_nuc_cyt <- merge(raw_nuc_data, raw_cyt_data, by = c("ImageNumber","ObjectNumber"))
  
  raw_nuc_cyt$Column <-  parse_number(raw_nuc_cyt$Metadata_Well.x)
  raw_nuc_cyt$Row <- substr(raw_nuc_cyt$Metadata_Well.x, 0,1)
  
  
  #Making a column of treatments
  
  raw_nuc_cyt <- raw_nuc_cyt %>%
    mutate(Treatment = case_when(Column == 3 ~ "BPA",
                                 Column == 4 ~ "BPA",
                                 Column == 5 ~ "Control",
                                 Column == 6 ~ "HD",
                                 Column == 7 ~ "HD",
                                 Column == 8 ~ "Control",
                                 Column == 9 ~ "DB",
                                 Column == 10 ~ "DB"))
  
  #Making a column of concentration
  
  raw_nuc_cyt <- raw_nuc_cyt %>%
    mutate(TreatConc = case_when(Row == "C" ~ "100 ug",
                                 Row == "D" ~ "050 ug",
                                 Row == "E" ~ "025 ug",
                                 Row == "F" ~ "012.5 ug"))
  
  raw_nuc_cyt$TreatConc[raw_nuc_cyt$Treatment == "Control"] <- "Control"
  
  ############
  ##Preparing the big, big dataframe!
  raw_nuc_cyt$TimePoint <- as.numeric(i)
  
  if(i == 1)
  {
    full_nuc_cyt <- raw_nuc_cyt
  }else{
    full_nuc_cyt <- rbind(full_nuc_cyt, raw_nuc_cyt)
  }
  
  raw_nuc_cyt$TimePoint <- NULL
  
  ############
  
  
  
  raw_nuc_cyt <- unite(raw_nuc_cyt, col = "TreatWithConc", c("Treatment", "TreatConc"), sep = '-')
  raw_nuc_cyt$TreatWithConc <- as.factor(raw_nuc_cyt$TreatWithConc)
  
  #write.csv(raw_nuc_cyt, paste(out_dir, "raw_nuc_cyt.csv", sep = "/"), quote = T, row.names = F)
  #to make a readable file for verification or reading of a single timepoint^
  
  #############
  #Cell counter
  if(i ==1) {
    TreatLevel <- levels(as.factor(raw_nuc_cyt$TreatWithConc)) 
    CellCounter <- count(raw_nuc_cyt, TreatWithConc)
    names(CellCounter)[2] <- paste0("TP",sprintf("%02d", i))
  }else
  {
    tv1 <- count(raw_nuc_cyt, TreatWithConc)
    CellCounter <- left_join(CellCounter, tv1)
    names(CellCounter)[i+1] <- paste0("TP",sprintf("%02d", i))
  }
  ##############
  
  ##############
  ##Making workable data
  
  #Works only if the custom CellProfiler pipeline output is used to generate the data files.
  #If the CellProfiler pipeline's output is modified, modify this portion of the code as well!
  #This is important so that the PCA doesn't take constants or strings into consideration
  
  nuc_cyt <- raw_nuc_cyt
  
  nuc_cyt[,3:17] <- NULL
  nuc_cyt[,52:66] <- NULL
  nuc_cyt$Row <- NULL
  nuc_cyt$Column <- NULL
  nuc_cyt$TreatWithConc <- NULL
  ##removes descriptors having the particular text in their name
  nuc_cyt <- nuc_cyt[,!grepl("Center", colnames(nuc_cyt))] 
  nuc_cyt <- nuc_cyt[,!grepl("Location", colnames(nuc_cyt))]
  nuc_cyt[,c('Number_Object_Number')] <- NULL
  nuc_cyt <- nuc_cyt[,!grepl("EulerNumber", colnames(nuc_cyt))]
  
  #############
  
  
  #############
  ##PCA!!!!
  
  ##Scaling relevant descriptors, here first 2 descriptors are ignored since they are address descriptors
  nuc_cyt.pca <- scale(nuc_cyt[,c(3:ncol(nuc_cyt))])
  nuc_cyt.pca[is.na(nuc_cyt.pca)] <- 0
  
  nuc_cyt.pca2 <- prcomp(nuc_cyt.pca, center = TRUE)
  
  #############
  
  #############
  ##Contributors plot
  
  ## Top 15 contributors to PC1 and PC2 combined
  cont1 <- fviz_contrib(nuc_cyt.pca2, choice = "var", axes = 1:2, top = 15) 
  
  #png(filename = paste(out_dir, paste0("TP", i ,"tryTop 15 contributors to PC1 and PC2.pdf")
  #                    , sep = "/"), width = 1920, height = 1080)
  #plot(cont1, main = 'Top 15 contributors to PC1 and PC2 at TimePoint 1')
  #dev.off()
  #Plots can be saved or ignored
  
  ##The importance of PCs visualized
  cont2 <- fviz_screeplot(nuc_cyt.pca2, ncp = 15, addlabels = T, hjust = -0.3) + 
    labs(title = paste0("TP",i,"- Importance of top 15 PCs"))
  # png(filename = paste(out_dir, paste0("TP",i,".Importance of top 15 PCs.png"),
  #  sep = "/"), width = 1920, height = 1080)
  #plot(cont2, main = 'Importance of top 15 PCs at TimePoint 1')
  #dev.off()
  
  png(filename = paste0(out_dir,"/PCAplot_TP",i,".png"), width = 1920, height = 1080, res = 100)
  PCnetplot <- grid.arrange(cont1, cont2, nrow = 1,
                            top = paste0("TP ",i))
  dev.off()
  # ggsave(filename = paste0(out_dir,"/PCAplot_TP",i,".pdf"),
  #      device = 'pdf', plot = PCnetplot)
  
  #############
  
  
  #############
  ##Finding top contributors for every timepoint
  #Storing contribution data for first two PCs' top 15 contributors
  cont3 <- fviz_contrib(nuc_cyt.pca2, choice = "var", axes = 1:2, top = 15)
  #Storing just the numbers from the array in a dataframe
  x <- cont3$data
  
  
  ##Every feature ranked according to its contribution
  x <- mutate(x, !!paste0("Rank.TP",sprintf("%02d", i)) := rank(-x$contrib, ties.method = "first"))
  if(i == 1)
  {
    y <- select(x, name, Rank.TP01)
  }else{
    foo <- select(x, starts_with("Rank.TP"))
    y <- cbind(y, foo)
  }
  
  #############
  
  
  #############
  #Nulling values to keep things clear and clean
  nuc_cyt.pca2 <- NULL
  nuc_cyt.pca <- NULL
  nuc_cyt <- NULL
  
  raw_nuc_data <- NULL
  raw_cyt_data <- NULL
  raw_nuc_cyt <- NULL
  
  
  f1 <- NULL
  f2 <- NULL
  
  cont1 <- NULL
  cont2 <- NULL
  x <- NULL
  cont3 <- NULL
  foo <- NULL
  
}

##Descriptor ranking table preparation, 

#To generate the ranking table of all the descriptors across all the timepoints, uncomment the below code
#write.csv(y, paste(out_dir, "DescriptorRanks.csv", sep = "/"), quote = T, row.names = F)

y <- column_to_rownames(y, var = "name")
y <- cbind(y, NetRank = rowSums(y))
y <- rownames_to_column(y, var = "Descriptors")

#########################

#Here, one can either use the ranking generated- NetRank which is sum-based, from the above dataframe (y),
#or use the the one below which is sum-based as well, although provides clearer insight on when and which 
#descriptor performed better since it classifies descriptors based on their ranks.


##Ranking algorithm, ranks the final descriptors based on their performance across all the timepoints  
#Storing rank classes


y.class <- y
y.class <- column_to_rownames(y.class, var = "Descriptors")

y.class[y.class <= 10] <- 1
y.class[y.class > 10 & y.class <= 20] <- 2
y.class[y.class > 20] <- 3
y.class$NetRank <- NULL
y.class <- rownames_to_column(y.class, var = "Descriptors")

y.class$Class1 <- rowSums(y.class[,2:31] == 1)
y.class$Class2 <- rowSums(y.class[,2:31] == 2)
y.class$Class3 <- rowSums(y.class[,2:31] == 3)

y.class <- arrange(y.class, -Class1, -Class2, -Class3)

#Select 'n' as the number of top descriptors necessary (USER PROVIDED)
n <- 10

TopDes <- c(y.class[1:n,1])


###Plotting cell counts

#Control count is divided by 4 in all cases to have the same scale of y-axis for every descriptor

angled.boxes <-
  list("far.from.others.borders","calc.boxes","enlarge.box","draw.rects")
#Angled boxes are the boxes with the name of the curve, allowing graphs to be more descriptive

#plot1 BPA

CellCounterBPA <- CellCounter[1:4,]
CellCounterBPA <- rbind(CellCounterBPA, CellCounter[5,])
TidyBPA <- gather(data = CellCounterBPA, key = TP, value = Count, -TreatWithConc)
TidyBPA$Count <- ifelse(TidyBPA$TreatWithConc == 'Control-Control', 
                        TidyBPA$Count/4, TidyBPA$Count)
x3 <- ggplot(TidyBPA, aes(x= TP, y = Count, group = TreatWithConc)) +
  geom_line(aes(color=TreatWithConc))+
  geom_point(aes(color=TreatWithConc)) +
  ggtitle(label = "Count of cells treated with Bisphenol A (BPA), over the period of time") +
  scale_x_discrete(breaks = c('TP05', 'TP10', 'TP15','TP20', 'TP25', 'TP30'))
x31 <- direct.label(x3, 'angled.boxes')
x31

ggsave(filename = paste0(out_dir,"/CellCountBPA.pdf"),
       device = 'pdf', dpi = 320)

#Plot2 HD
CellCounterHD <- CellCounter[10:13,]
CellCounterHD <- rbind(CellCounterHD, CellCounter[5,])
TidyHD <- gather(data = CellCounterHD, key = TP, value = Count, -TreatWithConc)
TidyHD$Count <- ifelse(TidyHD$TreatWithConc == 'Control-Control', 
                       TidyHD$Count/4, TidyHD$Count)
x4 <- ggplot(TidyHD, aes(x= TP, y = Count, group = TreatWithConc)) +
  geom_line(aes(color=TreatWithConc))+
  geom_point(aes(color=TreatWithConc)) +
  ggtitle(label = "Count of cells treated with Cetrimonium bromide (HD), over the period of time") +
  scale_x_discrete(breaks = c('TP05', 'TP10', 'TP15','TP20', 'TP25', 'TP30'))
x41 <- direct.label(x4, 'angled.boxes')
x41

ggsave(filename = paste0(out_dir,"/CellCountHD.pdf"),
       device = 'pdf', dpi = 320)

#Plot3 DB
CellCounterDB <- CellCounter[6:9,]
CellCounterDB <- rbind(CellCounterDB, CellCounter[5,])
TidyDB <- gather(data = CellCounterDB, key = TP, value = Count, -TreatWithConc)
TidyDB$Count <- ifelse(TidyDB$TreatWithConc == 'Control-Control', 
                       TidyDB$Count/4, TidyDB$Count)
x5 <- ggplot(TidyDB, aes(x= TP, y = Count, group = TreatWithConc)) +
  geom_line(aes(color=TreatWithConc))+
  geom_point(aes(color=TreatWithConc)) +
  ggtitle(label = "Count of cells treated with Dibutyltin dilaurate (DB), over the period of time") +
  scale_x_discrete(breaks = c('TP05', 'TP10', 'TP15','TP20', 'TP25', 'TP30'))
x51 <- direct.label(x5, 'angled.boxes')
x51

ggsave(filename = paste0(out_dir,"/CellCountDB.pdf"),
       device = 'pdf', dpi = 320)

png(filename = paste0(out_dir,"/NetplotCellCount.png"), width = 1920, height = 1080, res = 100)
netplot <- grid.arrange(x31, x41, x51, nrow = 2) 
dev.off()
# ggsave(filename = paste0(out_dir,"/CellCountNetPlot.pdf"),
#        device = 'pdf', dpi = 320, plot = netplot)


#################################

####Time series per treatment
spreadHD <- separate(TidyHD, col = TreatWithConc, sep = '-', into = c('Treatment', 'TreatConc'))
spreadHD$TP <- parse_number(spreadHD$TP)
colnames(spreadHD)[3] <- 'TimePoint'

#Here the descriptor values are multiplied with the cell counts to have clear visualization
#This multiplication can be ignored, but the graphs may NOT display the trend


#To save the individual graphs, uncomment the save graph option

for(i in 1:length(TopDes))
{
  HDnuccyt <- subset(full_nuc_cyt, Treatment == 'HD' | Treatment == 'Control', 
                     select = c(Treatment, TimePoint, TreatConc))
  HDnuccyt[,4] <- subset(full_nuc_cyt, Treatment == 'HD' | Treatment == 'Control', 
                         select = TopDes[i])
  
  spreadHD <- separate(TidyHD, col = TreatWithConc, sep = '-', into = c('Treatment', 'TreatConc'))
  spreadHD$TP <- parse_number(spreadHD$TP)
  colnames(spreadHD)[3] <- 'TimePoint'
  
  HDnuccyt <- left_join(HDnuccyt, spreadHD, by = c('Treatment', 'TreatConc', 'TimePoint'))
  
  HDnuccyt <- mutate(HDnuccyt, !!paste0("Cell Count * ",TopDes[i]) :=  HDnuccyt[,4]* Count)  
  
  t1 <- ggplot(HDnuccyt, aes(x = TimePoint, y = HDnuccyt[,6], color = TreatConc)) + 
    ylab(paste0("Cell Count x ",TopDes[i])) +
    geom_smooth( se = T, fullrange = T, aes(fill = TreatConc)) +
    labs(title = "Treatment HD") +
    theme(axis.title=element_text(size=7))
  
  #ggsave(filename = paste0(out_dir,"/NetTreatmentPlot_",TopDes[i],"HD.pdf"),
  #             device = 'pdf', plot = t1)
  
 
  
  #DB
  DBnuccyt <- subset(full_nuc_cyt, Treatment == 'DB' | Treatment == 'Control', 
                     select = c(Treatment, TimePoint, TreatConc))
  DBnuccyt[,4] <- subset(full_nuc_cyt, Treatment == 'DB' | Treatment == 'Control', 
                         select = TopDes[i])
  
  spreadDB <- separate(TidyDB, col = TreatWithConc, sep = '-', into = c('Treatment', 'TreatConc'))
  spreadDB$TP <- parse_number(spreadDB$TP)
  colnames(spreadDB)[3] <- 'TimePoint'
  
  DBnuccyt <- left_join(DBnuccyt, spreadDB, by = c('Treatment', 'TreatConc', 'TimePoint'))
  
  DBnuccyt <- mutate(DBnuccyt, !!paste0("Cell Count * ",TopDes[i]) :=  DBnuccyt[,4]* Count)  
  
  t2 <- ggplot(DBnuccyt, aes(x = TimePoint, y = DBnuccyt[,6], color = TreatConc)) + 
    ylab(paste0("Cell Count x ",TopDes[i])) +
    geom_smooth( se = T, fullrange = T, aes(fill = TreatConc)) +
    labs(title = "Treatment DB") +
    theme(axis.title=element_text(size=7))
  
  # ggsave(filename = paste0(out_dir,"/NetTreatmentPlot_",TopDes[i],"DB.pdf"),
  #        device = 'pdf', plot = t2)
  
  
  
  #BPA
  BPAnuccyt <- subset(full_nuc_cyt, Treatment == 'BPA' | Treatment == 'Control', 
                      select = c(Treatment, TimePoint, TreatConc))
  BPAnuccyt[,4] <- subset(full_nuc_cyt, Treatment == 'BPA' | Treatment == 'Control', 
                          select = TopDes[i])
  
  spreadBPA <- separate(TidyBPA, col = TreatWithConc, sep = '-', into = c('Treatment', 'TreatConc'))
  spreadBPA$TP <- parse_number(spreadBPA$TP)
  colnames(spreadBPA)[3] <- 'TimePoint'
  
  BPAnuccyt <- left_join(BPAnuccyt, spreadBPA, by = c('Treatment', 'TreatConc', 'TimePoint'))
  
  BPAnuccyt <- mutate(BPAnuccyt, !!paste0("Cell Count * ",TopDes[i]) :=  BPAnuccyt[,4]* Count)  
  
  t3 <- ggplot(BPAnuccyt, aes(x = TimePoint, y = BPAnuccyt[,6], color = TreatConc)) + 
    ylab(paste0("Cell Count x ",TopDes[i])) +
    geom_smooth( se = T, fullrange = T, aes(fill = TreatConc)) +
    labs(title = "Treatment BPA") +
    theme(axis.title=element_text(size=7))
  
 
  tnetplot <- grid.arrange(t1, t2, t3, nrow = 2,
                           top = paste0("Descriptor '" ,TopDes[i],"'  over the period of time"))
   
  ggsave(filename = paste0(out_dir,"/NetTreatmentPlot_",TopDes[i],".pdf"),
         device = 'pdf', plot = tnetplot)
}


####################################Credits: (unordered)##########################################
#https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
#https://stackoverflow.com/questions/5677885/ignore-outliers-in-ggplot2-boxplot
#https://digibio.blogspot.com/2016/09/box-plots-and-connect-by-median.html
#http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization
#https://stackoverflow.com/questions/5824173/replace-a-value-in-a-data-frame-based-on-a-conditional-if-statement
#https://stackoverflow.com/questions/37072844/select-subset-of-columns-based-on-vector-r
#https://stackoverflow.com/questions/5812493/how-to-add-leading-zeros
#https://stackoverflow.com/questions/26319567/use-grepl-to-search-either-of-multiple-substrings-in-a-text
#https://stackoverflow.com/questions/41815039/remove-columns-that-contain-a-specific-word
##https://stackoverflow.com/questions/43487773/how-to-rename-selected-columns-using-dplyr-with-new-column-names-as-strings
#https://stackoverflow.com/questions/34617883/how-to-remove-multiple-columns-in-r-dataframe)
#https://tolstoy.newcastle.edu.au/R/devel/06/05/5377.html
#https://rstudio-pubs-static.s3.amazonaws.com/3364_d1a578f521174152b46b19d0c83cbe7e.html
#http://directlabels.r-forge.r-project.org/docs/lineplot/posfuns/angled.boxes.html
#https://rpubs.com/williamsurles/293454


#Lecture slides by Maria 'Mia' Kjellson and Benjamin Guiastrennec
#Special thanks to Alboukadel Kassambara for his online tutorials on different websites and multiple R packages. 


