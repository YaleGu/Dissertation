install.packages('knitr')
install.packages('RCurl')
install.packages("forecast")
install.packages("ggfortify")
install.packages('nnet')
install.packages('spatialEco')
install.packages("corrplot")
install.packages("Hmisc")
library(Hmisc)
library(corrplot)
library(spatialEco)
library(nnet)
library(forecast)
library(ggfortify)
library(sp)
library(sf)
library(ggplot2)
library(raster)
library(dplyr)
library(tmaptools)
library(here)
library(maptools)
library(spdep)
library(rgdal)
library(tmap)
library(gridExtra)
library(gstat)
library(OpenStreetMap)
library(spacetime)
library(knitr)
library(RCurl)
library(reshape)
library(lattice)

map <- st_read(here::here("D:/Spatial Data Science/Dissertation/SHP/GroceryPoints_single.shp"))

map2 <- sf:::as_Spatial(map)

# Nearest Neighbour Index
nni(map2, win = "hull")

nni_result <- list()
for(i in unique(map2$NAME)){
  nni_result[[i]] <- nni(map2[map2$NAME == i,], win = "hull")
}
nni_summary <- as.data.frame(matrix(unlist(nni_result), nrow = length(unique(map2$NAME)), byrow = T))
names(nni_summary) <- c('NNI','z.score','p','expected.mean.distance','observed.mean.distance')

result <- as.data.frame(do.call(rbind, nni_result))
result <- as.matrix(result)

#write.csv(result,"D:/Spatial Data Science/Dissertation/SHP/NNIresult.csv")

rawdata <- read.csv("D:/Spatial Data Science/i2p/dissertation/correff1.csv")
df = subset(rawdata, select = -c(NAME) )



#import the map and merge attributes
map3 <- st_read(here::here("D:/Spatial Data Science/Dissertation/SHP/Borough_with_MobilityChange.shp"))
map3$time_Treat <- as.numeric(map3$time_Treat)
map3$Stage2 <- as.numeric(map3$Stage2)
map3$Stage3 <- as.numeric(map3$Stage3)

#plot maps(mobility change)
tm_shape(map3)+ 
  tm_fill("Stage3", style="order", palette="Spectral",title='Variation of mobility(Stage 3)')+
  tm_borders("grey")+
  tm_compass(position=c("left","top"))+
  tm_layout(attr.position=c("left", "bottom"),legend.title.size=1.5,
            legend.text.size = 1.1,legend.position=c('right','bottom'),frame=FALSE,legend.outside = TRUE,legend.outside.position='left')+tm_scale_bar(text.size = 0.7)

ggplot(map3) + 
  geom_sf(aes(fill=time_Treat))+
  scale_fill_brewer(palette = "RdYlBu") +
  coord_sf()+
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

#spatial autocorrelation
W <- nb2listw(poly2nb(map3))
W2 <- listw2mat(W)


map4 <- as_Spatial(map3)
data_matrix2<-data.matrix(map4@data)
data_matrix3 <- map4@data
time_treat <- as.numeric(data_matrix3$time_Treat)
moran.test(x=time_treat, listw=W)
moran.mc(x=time_treat, listw=W, nsim=9999)
lm <- localmoran(x=time_treat, listw=W)

stage2 <- as.numeric(data_matrix3$Stage2)
moran.test(x=stage2, listw=W)
moran.mc(x=stage2, listw=W, nsim=9999)
lm2 <- localmoran(x=stage2, listw=W)

stage3 <- as.numeric(data_matrix3$Stage3)
moran.test(x=stage3, listw=W)
moran.mc(x=stage3, listw=W, nsim=9999)
lm3 <- localmoran

lmMatrx <- read.csv("D:/Spatial Data Science/i2p/dissertation/localMoran.csv")
map5 <- merge(map3,lmMatrx,by.x='NAME',by.y='Borough')

# plot maps(Moran's I)

plot1 <- tm_shape(map5)+ 
  tm_fill("Stage1", style="order", palette="seq",title="Local Moran's I(Stage 1)")+
  tm_borders("grey")+
  tm_compass(position=c("left","top"))+
  tm_layout(attr.position=c("left", "bottom"),legend.title.size=1.5,
            legend.text.size = 1.1,legend.position=c('right','bottom'),frame=FALSE,legend.outside = TRUE,legend.outside.position='left',aes.palette = list(seq = "-RdBu"))+tm_scale_bar(text.size = 0.7)
plot2 <- tm_shape(map5)+ 
  tm_fill("Stage2.y", style="order", palette="seq",title="Local Moran's I(Stage 2)")+
  tm_borders("grey")+
  tm_compass(position=c("left","top"))+
  tm_layout(attr.position=c("left", "bottom"),legend.title.size=1.5,
            legend.text.size = 1.1,legend.position=c('right','bottom'),frame=FALSE,legend.outside = TRUE,legend.outside.position='left',aes.palette = list(seq = "-RdBu"))+tm_scale_bar(text.size = 0.7)
plot3 <- tm_shape(map5)+ 
  tm_fill("Stage3.y", style="order", palette="seq",title="Local Moran's I(Stage 3)")+
  tm_borders("grey")+
  tm_compass(position=c("left","top"))+
  tm_layout(attr.position=c("left", "bottom"),legend.title.size=1.5,
            legend.text.size = 1.1,legend.position=c('right','bottom'),frame=FALSE,legend.outside = TRUE,legend.outside.position='left',aes.palette = list(seq = "-RdBu"))+
  tm_scale_bar(text.size = 0.7) 
tmap_arrange(plot1,plot2,plot3,nrow=2)

# Correlation Matrix
CM1 <- read.csv("D:/Spatial Data Science/i2p/dissertation/correffStage3.csv")
cor_matrix <- cor(CM1,method = "pearson")
cor_matrix1<- as.data.frame(cor_matrix)

res2 <- rcorr(as.matrix(CM1))
# Generate the matrix with coefficient and graph
corrplot(cor_matrix,type="upper",tl.pos="ld",tl.col = "black",cl.pos = 'b')
corrplot(res2$r, type="lower", p.mat = res2$P, sig.level = 0.05, insig = "blank",tl.pos="ld",tl.col = "black",cl.pos = 'b',tl.srt=45)
#corrplot(cor_matrix,add=TRUE, type="lower", method="number", col="black",diag=FALSE,tl.pos="n", cl.pos="n")


# Plot maps(mobility change LSOA level)
LSOAmap <- st_read(here::here("D:/Spatial Data Science/Dissertation/SHP/LSOApredata.shp"))
LSOAmap$Stage1 <- as.numeric(LSOAmap$Stage1)
LSOAmap$Stage2 <- as.numeric(LSOAmap$Stage2)
LSOAmap$Stage3 <- as.numeric(LSOAmap$Stage3)
plota <- tm_shape(LSOAmap)+ 
  tm_fill("Stage1", style="cont", palette="Spectral",title="Mobility to groceries(Stage 1)",midpoint =0)+
  tm_borders("grey",alpha=0.3)+
  tm_compass(position=c("left","top"))+
  tm_layout(attr.position=c("left", "bottom"),legend.title.size=1.5,
            legend.text.size = 1.1,legend.position=c('right','bottom'),frame=FALSE,legend.outside = TRUE,legend.outside.position='left')+tm_scale_bar(text.size = 0.7)
plota
plotb <- tm_shape(LSOAmap)+ 
  tm_fill("Stage2", style="cont", palette="Spectral",title="Mobility to groceries(Stage 2)",midpoint =0)+
  tm_borders("grey",alpha=0.3)+
  tm_compass(position=c("left","top"))+
  tm_layout(attr.position=c("left", "bottom"),legend.title.size=1.5,
            legend.text.size = 1.1,legend.position=c('right','bottom'),frame=FALSE,legend.outside = TRUE,legend.outside.position='left')+tm_scale_bar(text.size = 0.7)
plotb
plotc <- tm_shape(LSOAmap)+ 
  tm_fill("Stage3", style="cont", palette="Spectral",title="Mobility to groceries(Stage 3)",midpoint =0)+
  tm_borders("grey",alpha=0.3)+
  tm_compass(position=c("left","top"))+
  tm_layout(attr.position=c("left", "bottom"),legend.title.size=1.5,
            legend.text.size = 1.1,legend.position=c('right','bottom'),frame=FALSE,legend.outside = TRUE,legend.outside.position='left')+tm_scale_bar(text.size = 0.7)
plotc
tmap_arrange(plota,plotb,plotc,nrow=3)

# plot maps(K-means Clustering)
kmeansR <- read.csv("D:/Spatial Data Science/i2p/dissertation/KmeansResult.csv")
LSOAmap <- merge(LSOAmap,kmeansR,by.x='LSOA11CD',by.y='LSOA')
tm_shape(LSOAmap)+ 
  tm_fill("Kmeans", style="cat", palette="Set3",title="Kmeans Clustering")+
  tm_borders("grey",alpha=0.3)+
  tm_compass(position=c("left","top"))+
  tm_layout(attr.position=c("left", "bottom"),legend.title.size=1.5,
            legend.text.size = 1.1,legend.position=c('right','bottom'),frame=FALSE,legend.outside = TRUE,legend.outside.position='left')+tm_scale_bar(text.size = 0.7)
tm_shape(LSOAmap)+ 
  tm_fill("LAD11CD", palette="Set3")+
  tm_borders("black",alpha=0.4)+
  tm_compass(position=c("left","top"))+
  tm_layout(attr.position=c("left", "bottom"),legend.title.size=0.5,
            legend.text.size = 0.5,legend.position=c('right','bottom'),frame=FALSE,legend.outside = TRUE,legend.outside.position='left')+tm_scale_bar(text.size = 0.5)
