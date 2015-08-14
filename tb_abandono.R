library(sensitivity)
library(rgl)
library(Rglpk)
library(pse)
library(FME)
require(ggplot2)
require(reshape)
require(shape)
require(epicalc)
#install.packages("sqldf")
library(sqldf)

## Abandono recuperado
Data_peru<-read.csv("tb_peru.csv")
Data1<-data.frame(years=Data_peru[52:75,2],abandonos=Data_peru[52:75,29])


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



### Abandono TB MDR
Data_abandono_tb_mdr<-read.csv("abandono_MDR.csv")
Data_abandono_tb_mdr$region<-factor(Data_abandono_tb_mdr$region)
list(Data_abandono_tb_mdr$region)



quartz (height=6,width=10, pointsize=9)
par(bg = 'cadetblue1')
xyplot(Data_abandono_tb_mdr$p_abandonos~ Data_abandono_tb_mdr$year| Data_abandono_tb_mdr$region , 
       data=Data_abandono_tb_mdr, ylab="Tasa de abandono TB MDR", 
       xlab="Year", main="",
       par.settings=list(fontsize=list(text=8, points=9)),xlim=c(2010:2014), add=T) 


s <- split(Data_abandono_tb_mdr, Data_abandono_tb_mdr$region)
lapply(s, function(x) colMeans(x[, c("p_abandonos","Abandonos")],na.rm = TRUE))


Data_abandono_tb_mdr$region <- reorder(Data_abandono_tb_mdr$region, Data_abandono_tb_mdr$p_abandonos, median, na.rm = TRUE)
quartz(width=10, height=6, pointsize=8)
par(bg = 'cadetblue1')
par(mar=c(5,4,4,4)+3)
boxplot(Data_abandono_tb_mdr$p_abandonos ~ Data_abandono_tb_mdr$region,axes=FALSE, 
        ylab="Tasa de abandonos TB MDR", xlab="",
        main="",col="blue", ylim=c(0,13))
orderVtr <- levels(Data_abandono_tb_mdr$region["scores"])
axis(1, at=seq(1, 34, by=1), labels=FALSE)
axis(2, at=seq(0, 13, by=1))
text(x = seq(1, 34, by=1), par("usr")[3] - 1, srt = 45,
     labels = as.vector(orderVtr), pos = 2, xpd = TRUE)
abline(h=6.766667, col="red")


newdata <- Data_abandono_tb_mdr[c(-34,-68,-102),]

newdata$region <- reorder(newdata$region, newdata$Abandonos, median, na.rm = TRUE)
quartz(width=10, height=6, pointsize=8)
par(bg = 'cadetblue1')
par(mar=c(5,4,4,4)+3)
boxplot(newdata$Abandonos ~ newdata$region,axes=FALSE, 
        ylab="Número de abandonos TB MDR", xlab="",
        main="",col="blue", ylim=c(0,250))
orderVtr <- levels(newdata$region["scores"])
axis(1, at=seq(1, 34, by=1), labels=FALSE)
axis(2, at=seq(0, 250, by=10))
text(x = seq(1, 34, by=1), par("usr")[3] - 5, srt = 45,
     labels = as.vector(orderVtr), pos = 2, xpd = TRUE)

### abandono TB MDR por establecimientos de salud

Data_abandono_eess<-read.csv("Cohorte MDR_PERU_ 2009-2011_EESS.csv")
str(Data_abandono_eess)
tab1(Data_abandono_eess$DIRECCION_SALUD)

Data_abandono_eess_Lima <- Data_abandono_eess[ which(Data_abandono_eess$DIRECCION_SALUD=='CALLAO' 
                         | Data_abandono_eess$DIRECCION_SALUD=='LIMA CIUDAD' 
                         | Data_abandono_eess$DIRECCION_SALUD=='LIMA ESTE'
                         | Data_abandono_eess$DIRECCION_SALUD=='LIMA SUR' 
                         | Data_abandono_eess$DIRECCION_SALUD=='REGION LIMA'), ]


head(Data_abandono_eess_Lima)
tail(Data_abandono_eess_Lima)


join_string <- "select
                Data_abandono_eess_Lima.DIRECCION_SALUD
              , Data_abandono_eess_Lima.RED_SALUD
              , Data_abandono_eess_Lima.MICRO_RED
              , sum(Data_abandono_eess_Lima.Abandono) as Sum_abandono
              from Data_abandono_eess_Lima 
              group by Data_abandono_eess_Lima.DIRECCION_SALUD , Data_abandono_eess_Lima.RED_SALUD, Data_abandono_eess_Lima.MICRO_RED"

sqldf(join_string)


attach(Data_abandono_eess_Lima)
aggdata <-aggregate(Data_abandono_eess_Lima$Abandono, by=list(RED_SALUD), 
                    FUN=sum, na.rm=TRUE)
print(aggdata)

aggdata <-aggregate(Data_abandono_eess_Lima$Abandono, by=list(RED_SALUD, MICRO_RED), 
                    FUN=sum, na.rm=TRUE)
print(aggdata)


aggdata <-aggregate(Data_abandono_eess_Lima$Abandono, by=list(RED_SALUD, MICRO_RED,DESC_ESTAB), 
                    FUN=sum, na.rm=TRUE)
print(aggdata)
detach(Data_abandono_eess_Lima)



tab1(Data_abandono_eess$DIRECCION_SALUD)





