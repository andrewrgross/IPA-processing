######################################################################################
###### Proteomic analysis for networks, 2016-02-17

########################################################################
### Header
library(ggplot2)
library(Hmisc)

######################################################################################
### Functions

convert.char.column.to.num <- function(dataframe,column) {
  new.vector <- c()
  for (value in dataframe[,column]) {
    new.value <- as.numeric(strsplit(as.character(value),"/")[[1]][1])
    new.vector <- c(new.vector,new.value)
  }
  dataframe[column] <- new.vector
  return(dataframe)
}
######################################################################################
### Upload data

setwd(dir = "Z:/Data/zika/Andrew/IPA results/")
diseases <- read.csv('Disease_5dpi_all.csv')

######################################################################################
### Format

### Convert p-value to log10 scale
diseases <- diseases[c(2,3,5,4,7)]
diseases$p.Value <- log10(diseases$p.Value)

### Remove unwanted rows
diseases <- diseases[-c(5,6),]

### Rename rows
row.names(diseases) <- seq(1,16)

### Capitalize functions
diseases$Diseases.or.Functions.Annotation <- capitalize(as.character(diseases$Diseases.or.Functions.Annotation))

### Label unknown as mixed
diseases$Predicted.Activation.State <- as.character(diseases$Predicted.Activation.State)
diseases$Predicted.Activation.State[12:15] <- 'Mixed'

### Specify text positions
diseases$vjust <- 2
diseases$hjust <- 0.5
diseases$vjust[c(4,6,7,8)] <- 0.2
diseases$hjust[c(4,6,7,8)] <- -.05
diseases$vjust[c(4)] <- 0

diseases$vjust[c(10,11)] <- 0.3
diseases$hjust[c(10,11)] <- 1.1
diseases$vjust[c(5,9)] <- 2
diseases$hjust[c(5,9)] <- -0.03

### Specify text size
diseases$size <- 120
diseases$size[c(4:9)] <- 40

######################################################################################
### Plot simple circles

g <- ggplot(data = diseases, aes(x = abs(Activation.z.score), y = p.Value)) +
  geom_point(aes(size = X..Molecules, color = Predicted.Activation.State)) + 
  scale_color_manual(values = c('Mixed' = 'grey', 'Decreased' = 'forestgreen', 'Increased' = 'red'),
                     name = 'Activation State') +
  xlim(-.25,4) +
  ylim(-2, -20) +
  theme(panel.grid = element_blank(), panel.background = element_rect(size = 2, color = 'black', fill = 'white'),
        plot.title = element_text(size = 16, face = 'bold', margin = margin(10,0,15,0), hjust = 0.5),
        axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12),
        legend.key = element_rect(fill = 'white'))

h <- g + geom_text(data = diseases, aes(x = abs(Activation.z.score), y = p.Value, label = Diseases.or.Functions.Annotation,
                       hjust = hjust, vjust = vjust, size = size)) + scale_radius(range = c(2,6), name = 'Genes Detected') +
  labs(title = 'Diseases and Functions',
       x = 'Activation Z-Score',
       y = 'p-Value')
h



########################################################################
### Export plot
########################################################################

setwd("z:/Data/zika/Andrew/")

filename <- "disease_function-14"

png(filename=paste0(filename,".png"), 
    type="cairo",
    units="in", 
    width=8, 
    height=8, 
    pointsize=12, 
    res=300)
h
dev.off()


