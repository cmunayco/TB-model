library(sensitivity)
library(rgl)
library(Rglpk)
library(pse)
library(FME)
require(ggplot2)
require(reshape)
require(shape)
require(epicalc)



##databases

tb_mdr_cox<-read.csv("rep_200_Abandono TB Resistente 2009-2014 260615.csv")
tb_mdr_cox$periodo_estudio[tb_mdr_cox$anoini<2012]<-"Cohortes 2009-2011"
tb_mdr_cox$periodo_estudio[tb_mdr_cox$anoini>=2012]<-"Cohortes 2012-2014"

tab1(tb_mdr_cox$periodo_estudio)
### Descriptive analysis

tab1(tb_mdr_cox$edad_inicio, graph=F)
tab1(tb_mdr_cox$sexo, graph=F)
tab1(tb_mdr_cox$condicion_egreso)
tab1(tb_mdr_cox$condicion_egreso_ultimo_tto)



quartz(width=8, height=6, pointsize=9)
par(bg = 'cadetblue1')
with(subset(tb_mdr_cox,tb_mdr_cox$sexo=="m"), hist(edad_actual,col="blue", main="", xlab="Edad (años)", ylab="Frecuencia"))
with(subset(tb_mdr_cox,tb_mdr_cox$sexo=="f"), hist(edad_actual,col="red", add = TRUE))
legend(65, 500, c("Hombres", "Mujeres"),
       fill = c("blue", "red"))

quartz(width=8, height=6, pointsize=9)
par(bg = 'cadetblue1')
par(mfrow=c(2,1))
with(subset(tb_mdr_cox,tb_mdr_cox$sexo=="m"), plot(density(edad_actual,na.rm =T),col="blue", main="Hombres", xlab="Edad (años)", ylab="Frecuencia"))
with(subset(tb_mdr_cox,tb_mdr_cox$sexo=="f"), plot(density(edad_actual,na.rm =T),col="red", main="Mujeres", xlab="Edad (años)", ylab="Frecuencia"))


quartz(width=8, height=6, pointsize=9)
par(bg = 'cadetblue1')
with(subset(tb_mdr_cox,tb_mdr_cox$anoini<2012), hist(edad_actual,col="blue", main="", xlab="Edad (años)", ylab="Frecuencia"))
with(subset(tb_mdr_cox,tb_mdr_cox$anoini>=2012), hist(edad_actual,col="red", add = TRUE))
legend(70, 400, c("Cohortes 2009-2011", "Cohortes 2012-2014"),
       fill = c("blue", "red"))

quartz(width=8, height=6, pointsize=9)
par(bg = 'cadetblue1')
par(mfrow=c(2,1))
with(subset(tb_mdr_cox,tb_mdr_cox$anoini<2012), plot(density(edad_actual,na.rm =T),col="blue", main="Cohortes 2009-2011", xlab="Edad (años)", ylab="Frecuencia"))
with(subset(tb_mdr_cox,tb_mdr_cox$anoini>=2012), plot(density(edad_actual,na.rm =T),col="red", main="Cohortes 2012-2014", xlab="Edad (años)", ylab="Frecuencia"))


tab1(tb_mdr_cox$lab_vih)
tab1(tb_mdr_cox$p_resistencia)
tab1(tb_mdr_cox$finicio_rafa)
tab1(tb_mdr_cox$fecapr_rafa)
tab1(tb_mdr_cox$initial_treatment_type)
tab1(tb_mdr_cox$num_tto_previo)
tab1(tb_mdr_cox$dosis_egreso)
tab1(tb_mdr_cox$ultimo_tto)
tab1(tb_mdr_cox$condicion_egreso)
tab1(tb_mdr_cox$condicion_egreso_ultimo_tto)
tab1(tb_mdr_cox$traant)
tab1(tb_mdr_cox$tto_actual_special_status.1)
tab1(tb_mdr_cox$dosis_transferencia)
tab1(tb_mdr_cox$estatus)
tab1(tb_mdr_cox$REGION.DE.SALUD)
tab1(tb_mdr_cox$DESC_DPTO)


table(tb_mdr_cox$lab_vih,tb_mdr_cox$periodo_estudio)
quartz(width=8, height=6, pointsize=9)
barplot(prop.table(table(tb_mdr_cox$lab_vih,tb_mdr_cox$periodo_estudio),2)*100, ylim=c(0,100), 
        col=c("green", "blue", "red"),xlab="Cohortes",space=0.5)
legend("top", c("No data", "VIH -", "VIH +"),
       fill = c("green", "blue", "red"), box.lty = 0, cex=1.5)
