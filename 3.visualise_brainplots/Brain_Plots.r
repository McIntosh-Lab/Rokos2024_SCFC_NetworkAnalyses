options(repos = c(ggseg='https://ggseg.r-universe.dev',
                  CRAN = 'https://cloud.r-project.org'))
install.packages('ggplot2')
install.packages('wesanderson')
install.packages('readxl')
install.packages('tidyverse')
install.packages('dplyr')
install.packages('ggseg')
install.packages('ggsegSchaefer')

library(readxl)
library(ggseg)
library(ggsegSchaefer)
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(dplyr)

### Step 1: Set file paths and names
data_dir <- "~/path/to/directory/Rokos2024_SCFC_NetworkAnalyses" #MODIFY
setwd(data_dir)
plot_data_file <- file.path(data_dir, "outputs", "FIGURE4A_bsrs.csv") # MODIFY
output_file <- file.path(data_dir, "outputs", "figures", "FIGURE4A_brain.png") # MODIFY
labels_file <- file.path(data_dir, "3.visualise_brainplots", "ST193_NetLabels.xlsx")

## Step 2: Load in the Data
labels <- read_excel(labels_file)
results = data.frame(labels)
colnames(results)[2]  <- "region"
colnames(results)[3]  <- "hemi"

PlotData <- read_csv(plot_data_file,col_names=c('bsrs')) # Data with a value in each row.
plot.df<-cbind(results, PlotData[,])
plot.df$bsrs[plot.df$bsrs == 0] <- NA
plot.df[ plot.df == "NaN" ] <- NA #setting the NA excel values to actual NAs

pal <- wes_palette("Zissou1", 50, type = "continuous")

newdata <- subset(plot.df, bsrs!= "NA")
someData <- tibble(
  region = c(newdata$region), 
  p = c(as.double(newdata$bsrs)),
  groups = c(newdata$hemi)
)

## Step 3: Plot
sp<-someData%>%
  ggseg(atlas = schaefer17_200,
        mapping=aes(fill=as.double(p)),
        show.legend = TRUE,
        position="stacked", colour="black")+
  scale_color_manual(values = pal)+
  theme(legend.title=element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 24))
sp+scale_fill_gradientn(colours = pal, limits=c(-5,5)) 
ggsave(output_file, sp+scale_fill_gradientn(colours = pal, limits=c(-5,5)), width = 8, height=6, dpi=600)


