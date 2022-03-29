library(readxl)
library(CRTspillover)
df1 <- read_excel('C:/Users/smith/Downloads/clusters/Bago(12).xls')
df2 <- read_excel('C:/Users/smith/Downloads/clusters/Chasimba(5).xls')
df3 <- read_excel('C:/Users/smith/Downloads/clusters/Kiwangwa(15).xls')
df4 <- read_excel('C:/Users/smith/Downloads/clusters/Kiwangwa(16).xls')
df5 <- read_excel('C:/Users/smith/Downloads/clusters/Kiwangwa(17).xls')
df6 <- read_excel('C:/Users/smith/Downloads/clusters/Kiwangwa(18).xls')
df7 <- read_excel('C:/Users/smith/Downloads/clusters/Kiwangwa(19).xls')
df8 <- read_excel('C:/Users/smith/Downloads/clusters/Kiwangwa(20).xls')
df9 <- read_excel('C:/Users/smith/Downloads/clusters/Kongo(3).xls')
df10 <- read_excel('C:/Users/smith/Downloads/clusters/Kongo(4).xls')
df11 <- read_excel('C:/Users/smith/Downloads/clusters/Masuguru(14).xls')
df12 <- read_excel('C:/Users/smith/Downloads/clusters/Matimbwa(1).xls')
df13 <- read_excel('C:/Users/smith/Downloads/clusters/Matimbwa(2).xls')
df14 <- read_excel('C:/Users/smith/Downloads/clusters/Msata (21).xls')
df15 <- read_excel('C:/Users/smith/Downloads/clusters/Msata (22).xls')
df16 <- read_excel('C:/Users/smith/Downloads/clusters/Msata (23).xls')
df17 <- read_excel('C:/Users/smith/Downloads/clusters/Msata (24).xls')
df18 <- read_excel('C:/Users/smith/Downloads/clusters/Msinune(7).xls')
df19 <- read_excel('C:/Users/smith/Downloads/clusters/Msinune(8).xls')
df20 <- read_excel('C:/Users/smith/Downloads/clusters/Mwavi(10).xls')
df21 <- read_excel('C:/Users/smith/Downloads/clusters/Mwavi(11).xls')
df22 <- read_excel('C:/Users/smith/Downloads/clusters/Mwavi(9).xls')
df23 <- read_excel('C:/Users/smith/Downloads/clusters/Yombo(6).xls')

addcluster=function(dfx,df=df){
  #standardise the variable names
  names(dfx)[names(dfx) == "v05southg"] <- "lat"
  names(dfx)[names(dfx) == "v05eastgp"] <- "long"
  names(dfx)[names(dfx) == "southgp"] <- "lat"
  names(dfx)[names(dfx) == "eastgp"] <- "long"
  #select only the required variables
  dfx=dfx[,c('FID','lat','long')]
  cno=max(df$cluster)+1
  #add cluster number
	dfx$cluster=cno
	dfx=as.data.frame(dfx)
	#append to data frame
	df=rbind(df,dfx)
return(df)}

cno=1
names(df1)[names(df1) == "v05southg"] <- "lat"
names(df1)[names(df1) == "v05eastgp"] <- "long"
df = df1[,c('FID','lat','long')]
  df$cluster=1
  df=addcluster(dfx=df2,df=df)
  df=addcluster(dfx=df3,df=df)
  df=addcluster(dfx=df4,df=df)
  df=addcluster(dfx=df5,df=df)
  df=addcluster(dfx=df6,df=df)
  df=addcluster(dfx=df7,df=df)
  df=addcluster(dfx=df8,df=df)
  df=addcluster(dfx=df9,df=df)
  df=addcluster(dfx=df10,df=df)
  df=addcluster(dfx=df11,df=df)
  df=addcluster(dfx=df12,df=df)
  df=addcluster(dfx=df13,df=df)
  df=addcluster(dfx=df14,df=df)
  df=addcluster(dfx=df15,df=df)
  df=addcluster(dfx=df16,df=df)
  df=addcluster(dfx=df17,df=df)
  df=addcluster(dfx=df18,df=df)
  df=addcluster(dfx=df19,df=df)
  df=addcluster(dfx=df20,df=df)
  df=addcluster(dfx=df21,df=df)
  df=addcluster(dfx=df22,df=df)
  df=addcluster(dfx=df23,df=df)
  df$lat=unlist(df$lat)
  plot(df$long,df$lat)

 # Correct outlying points that are obvious typos in entering Lat Long
 df$lat[df$lat > -6] = -6.30998
 df$long[df$long < 36] = 38.59039
 df$long[df$long < 37] = 38.61693

 # Create a variable for flagging outliers for further investigation
 df$flag = FALSE

plotcluster = function(df,cno){
  plot(df$long[df$cluster==cno & !df$flag],df$lat[df$cluster==cno & !df$flag])
}
plotcluster(df=df,cno = 23)

df$flag[df$cluster==2 & df$lat > -6.5]= TRUE
df$flag[df$cluster==5 & df$long < 38.4]= TRUE
df$flag[df$cluster==7 & df$long < 38.4]= TRUE
df$flag[df$cluster==9 & df$long < 38.6]= TRUE
df$flag[df$cluster==9 & df$lat > -6.52]= TRUE
df$flag[df$cluster==10 & df$long < 38.7]= TRUE
df$flag[df$cluster==10 & df$lat < -6.6]= TRUE
df$flag[df$cluster==10 & df$long > 38.87]= TRUE # may be outlying houses, check distance
df$flag[df$cluster==11 & df$long < 38.3]= TRUE # may be outlying houses, check distance
df$flag[df$cluster==12 & df$lat > -6.4]= TRUE
df$flag[df$cluster==12 & df$long < 38.8]= TRUE
df$flag[df$cluster==13 & df$long < 38.7]= TRUE
df$flag[df$cluster==19 & df$long < 38.45]= TRUE
df$flag[df$cluster==19 & df$lat > -6.3]= TRUE
df$flag[df$cluster==19 & df$lat < -6.48]= TRUE
df$flag[df$cluster==21 & df$long > 38.8]= TRUE
df$flag[df$cluster==23 & df$long < 38.6]= TRUE

df = Convert_LatLong(df)
plot(df$x[!df$flag],df$y[!df$flag])
plot(df$x[df$flag],df$y[df$flag])

centroids = df %>%
  group_by(cluster) %>%
  summarise(
    mean_x = mean(x, na.rm=T),
    mean_y = mean(y, na.rm=T))
merged = merge(df[, c("FID","x","y","flag","cluster")], centroids, by="cluster")
remoteness = with(merged,sqrt((x-mean_x)^2 + (y-mean_y)^2))
angle = asin((merged$y-merged$mean_y)/remoteness)

# Analysis of the direction and distance does not suggest any clear explanation for
# outliers
# For the moment remove the points that are flagged
df=df[!df$flag,]

# use Plot_CRTmap;

Plot_CRTmap(trial=df,showLocations=TRUE,showClusterLabels = TRUE,colourClusters=FALSE)
Plot_CRTmap(trial=df,showLocations=FALSE,colourClusters=TRUE)

