library(sensitivity)
library(rgl)
library(Rglpk)
library(pse)
library(FME)
require(ggplot2)
require(reshape)


## Real data
Data_peru<-read.csv("tb_peru.csv")
Data1<-data.frame(years=Data_peru[52:75,2],abandonos=Data_peru[52:75,29])
Data_abandono_tb_mdr<-read.csv("abandono_MDR.csv")


quartz(width=10, height=6, pointsize=10)
par(bg = 'cadetblue1')
plot(Data1$years, Data1$abandonos, ylab="% Abandonos Recuperados",col="red", pch=15, 
     xlab="Años", main="MINSA y otras instituciones", ylim=c(0,13),cex=1.5)

reduction_change<-Data1$abandonos[-1] - Data1$abandonos[-length(Data1$abandonos)]
years<-c("1991-1992", "1992-1993", "1993-1994","1994-1995", "1995-1996","1996-1997",
         "1997-1998","1998-1999","1999-2000", "2000-2001","2001-2002", "2002-2003",
         "2003-2004","2004-2005","2005-2006","2006-2007","2007-2008","2008-2009",
         "2009-2010","2010-2011","2011-2012","2012-2013","2013-2014")
change_abandono<-data.frame(years,reduction_change)

quartz(width=8, height=6, pointsize=9)
#par(las=2) # make label text perpendicular to axis
#par(mar=c(5,8,4,2)) # increase y-axis margin.
par(bg = 'cadetblue1')
barplot(change_abandono$reduction_change, width = 0.5, ylab="Cambio de la tasa de abandono recuperado",
        xlab="Años",ylim=c(-3,2), xlim=c(1,23),
        col="blue",cex.names=0.7,space=1)
axis(side = 2)
#ticks = c(seq(from=1,to=23,by=1))
axis(side = 1, at =change_abandono$years, labels=change_abandono$years, cex=0.7)



