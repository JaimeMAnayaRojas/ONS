# Clean environment
rm(list=ls(all=TRUE))
require(MuMIn)
require(partR2)
require(stringr)
require(plyr)
require(lme4)
require(lmerTest)
require(Hmisc)



#Upload Data
data = read.csv("data/data.csv")


##################################################################################################################################
##################################################################################################################################
#Clean up stream names and process Nitrogen data
##################################################################################################################################
##################################################################################################################################

data$count <- 1

for (i in 1:dim(data)[1]){
  data[i,'stream'] <- str_split(data$Full_name[i],'-')[[1]][1]
}
unique(data$stream)
data[which(data$stream=='El_Cedro'),'drainage'] <- 'El_Cedro'

data$diff.d15 <- (data$d15N - data$d15Nb)

data$TP <- data$diff.d15 


##################################################################################################################################
##################################################################################################################################
#Process delta 13C
##################################################################################################################################
##################################################################################################################################

D <- 7.018
I <- 0.048

data$L <- 93 / ( 1 + (0.246*data$CNratio - 0.775)^(-1)) 
data$d13Cc <- data$d13C + D * (I + (3.9/(1+287/data$L))) 

data$diff.d13Cc <- data$d13Cc - data$d13Cb


#Center on 10 mm
center.size <- 10

data$Length.cen <- (data$Length - center.size)

##################################################################################################################################
##################################################################################################################################
#Question 1 and 2. KG-N sites
##################################################################################################################################
##################################################################################################################################

#Do guppies and killifish have ONS's and is the killifish one stronger than guppies in KG nat?
#Nitrogen
#Full model
modN.full <- lmer(TP ~ 0 + FG + FG:Length.cen + (1 + Length.cen|drainage/FG) , data=data[which(data$Habitat=='KG.nat'),] )
summary(modN.full)
#Failed to converge

#Removed random effect of length
modN <- lmer(TP ~ 0 + FG + FG:Length.cen + (1 |drainage/FG) , data=data[which(data$Habitat=='KG.nat'),] )
summary(modN)
#Singular random effects


# #Removed nested random effects
# modN <- lmer(TP ~ 0 + FG + FG:Length.cen + (1 |drainage:FG) , data=data[which(data$Habitat=='KG.nat'),] )
# summary(modN)
# #Identical to singular random effects above.



randN <- ranef(modN)$`FG:drainage`
randN$FG <- sapply(strsplit(dimnames(randN)[[1]],split=':'), "[", 1)
randN$drainage <- sapply(strsplit(dimnames(randN)[[1]],split=':'), "[", 2)
names(randN) <- c('rand','FG','drainage')

coef <- summary(modN)$coef[,'Estimate']


#natural site slope comparison
L <- array(0,c(length(coef)))
L[which(names(coef)=='FGGuppy:Length.cen')] <- 1
L[which(names(coef)=='FGKillifish:Length.cen')] <- -1


#This is the test using lmerTest
contest(modN,t(L))
#This converts the reported F to a t-value
sqrt(contest(modN,t(L))$`F value`)


#natural site int comparison
L <- array(0,c(length(coef)))
L[which(names(coef)=='FGGuppy')] <- 1
L[which(names(coef)=='FGKillifish')] <- -1

#This is the test
contest(modN,t(L))
#This converts the reported F to a t-value
sqrt(contest(modN,t(L))$`F value`)



#Carbon
#Full Model with length random effect
modC.full <- lmer(diff.d13Cc ~ 0 + FG + FG:Length.cen + (1 + Length.cen|drainage/FG) , data=data[which(data$Habitat=='KG.nat'),] )
summary(modC.full)
#does not converge

#Without length random effect
modC <- lmer(diff.d13Cc ~ 0 + FG + FG:Length.cen + (1 |drainage/FG) , data=data[which(data$Habitat=='KG.nat'),] )
summary(modC)


r.squaredGLMM(modC)

randC <- ranef(modC)$`FG:drainage`

randC$FG <- sapply(strsplit(dimnames(randC)[[1]],split=':'), "[", 1)
randC$drainage <- sapply(strsplit(dimnames(randC)[[1]],split=':'), "[", 2)
names(randC) <- c('randC','FG','drainage')



coef <- summary(modC)$coef[,'Estimate']



#natural site slope comparison
L <- array(0,c(length(coef)))
L[which(names(coef)=='FGGuppy:Length.cen')] <- 1
L[which(names(coef)=='FGKillifish:Length.cen')] <- -1

#This is the test using lmerTest
contest(modC,t(L))
sqrt(contest(modC,t(L))$`F value`)


#natural site int comparison
L <- array(0,c(length(coef)))
L[which(names(coef)=='FGGuppy')] <- 1
L[which(names(coef)=='FGKillifish')] <- -1


#This is the test using lmerTest
contest(modC,t(L))
sqrt(contest(modC,t(L))$`F value`)


## Plot and save figure 1
jpeg(file = "Figures/Figure 1.jpeg",width = 9, height = 5, units='in', res=600)
par(mfrow=c(1,2), mar = c(4, 5 ,2, 1), oma = c(0.5, 1, 1, 0.5))

toplot <- data[which(data$Habitat=='KG.nat'),]

toplot <- merge(x=toplot,y=randN,by=c('drainage','FG'),all.x=TRUE)
toplot$TPc <- toplot$TP - toplot$rand

min <- -0.1#min(toplot$TPc) - 0.1*min(toplot$TPc)
max <- 7#max(toplot$TPc) + 0.1*max(toplot$TPc)

plot(toplot[which(toplot$FG=='Killifish'),'Length'],toplot[which(toplot$FG=='Killifish'),'TPc'],
     ylim=c(min,max),xlim=c(min(toplot$Length,na.rm=TRUE),max(toplot$Length,na.rm=TRUE)),pch=16,col='orange',bty='L',
     ylab=expression(paste(Delta^{15}, "N (\u2030)")),xlab='Fish Standard Length (mm)')
points(toplot[which(toplot$FG=='Guppy'),'Length'],toplot[which(toplot$FG=='Guppy'),'TPc'],pch=16,col='grey')

gsizes <- seq(10,max(toplot[which(toplot$FG=='Guppy'),'Length'],na.rm=TRUE),by=1)
ksizes <- seq(10,max(toplot[which(toplot$FG=='Killifish'),'Length'],na.rm=TRUE),by=1)
lines(ksizes,
      (summary(modN)$coef['FGKillifish','Estimate'] + summary(modN)$coef['FGKillifish:Length.cen','Estimate']*(ksizes-center.size)),
      lwd=3,col='orange')

lines(gsizes,
      (summary(modN)$coef['FGGuppy','Estimate'] + summary(modN)$coef['FGGuppy:Length.cen','Estimate']*(gsizes-center.size)),
      lwd=3,col='black',lty=3)
title(main='A) Trophic position',adj = 0.1, line = 0)
legend('topleft',c('Guppy','Killifish'),lty=c(1,1),pch=c(16,16),col=c('black','orange'),bty='n')



min <- -2.5
max <- 7

toplot <- merge(x=toplot,y=randC,by=c('drainage','FG'),all.x=TRUE)
toplot$d13Ccc <- toplot$diff.d13Cc - toplot$randC
plot(toplot[which(toplot$FG=='Killifish'),'Length'],toplot[which(toplot$FG=='Killifish'),'d13Ccc'],
     ylim=c(min,max),xlim=c(min(toplot$Length,na.rm=TRUE),max(toplot$Length,na.rm=TRUE)),pch=16,col='orange',bty='L',
     ylab=expression(paste(Delta^{13}, "C (\u2030)")),xlab='Fish Standard Length (mm)')
points(toplot[which(toplot$FG=='Guppy'),'Length'],toplot[which(toplot$FG=='Guppy'),'d13Ccc'],pch=16,col='grey')

gsizes <- seq(10,max(toplot[which(toplot$FG=='Guppy'),'Length'],na.rm=TRUE),by=1)
ksizes <- seq(10,max(toplot[which(toplot$FG=='Killifish'),'Length'],na.rm=TRUE),by=1)
lines(ksizes,
      (summary(modC)$coef['FGKillifish','Estimate'] + summary(modC)$coef['FGKillifish:Length.cen','Estimate']*(ksizes-center.size)),
      lwd=3,col='orange')

lines(gsizes,
      (summary(modC)$coef['FGGuppy','Estimate'] + summary(modC)$coef['FGGuppy:Length.cen','Estimate']*(gsizes-center.size)),
      lwd=3,col='black',lty=3)
title(main='B) Carbon source',adj = 0.1, line = 0)

dev.off()




##################################################################################################################################
##################################################################################################################################
#Question 3. Compare Communities
##################################################################################################################################
##################################################################################################################################

#Guppies
#Nitrogen
#Full model with length random effect
modN.full <- lmer(TP ~ 0 + FG:Habitat + FG:Habitat:Length.cen +  (1 + Length.cen|drainage/Habitat/FG), data=data[which(data$stream!='Endler' & data$stream!= 'Taylor' & data$stream!= 'Caigual' & data$stream!='El_Cedro' & (data$Habitat=='KG.nat' | data$Habitat=='KGP' | data$Habitat=='KO')),] )
summary(modN.full)


#Reduced model without length random effect
modN <- lmer(TP ~ 0 + FG:Habitat + FG:Habitat:Length.cen +  (1 |drainage/Habitat/FG), data=data[which(data$stream!='Endler' & data$stream!= 'Taylor' & data$stream!= 'Caigual' & data$stream!='El_Cedro' & (data$Habitat=='KG.nat' | data$Habitat=='KGP' | data$Habitat=='KO')),] )
summary(modN)
#Singular, but okay

anova(modN,modN.full,test='LRT')


r.squaredGLMM(modN)

coefN <- summary(modN)$coef[,'Estimate']

#natural site int comparison
L <- array(0,c(length(coefN)))
L[which(names(coefN)=='FGGuppy:HabitatKG.nat')] <- 1
L[which(names(coefN)=='FGGuppy:HabitatKGP')] <- -1


#This is the test using lmerTest
contest(modN,t(L))


#natural site slope comparison
L <- array(0,c(length(coefN)))
L[which(names(coefN)=='FGGuppy:HabitatKG.nat:Length.cen')] <- 1
L[which(names(coefN)=='FGGuppy:HabitatKGP:Length.cen')] <- -1


#This is the test using lmerTest
contest(modN,t(L))
sqrt(contest(modN,t(L))$`F value`)




#Killifish
#natural site int comparison
L <- array(0,c(length(coefN),2))
L[which(names(coefN)=='FGKillifish:HabitatKG.nat'),] <- 1
L[which(names(coefN)=='FGKillifish:HabitatKGP'),1] <- -1
L[which(names(coefN)=='FGKillifish:HabitatKO'),2] <- -1

#This is the test using lmerTest
contest(modN,t(L))


#natural site slope comparison
L <- array(0,c(length(coefN),2))
L[which(names(coefN)=='FGKillifish:HabitatKG.nat:Length.cen'),] <- 1
L[which(names(coefN)=='FGKillifish:HabitatKGP:Length.cen'),1] <- -1
L[which(names(coefN)=='FGKillifish:HabitatKO:Length.cen'),2] <- -1


#This is the test using lmerTest
contest(modN,t(L))








#Carbon
#Full model with random slopes
modC.full <- lmer(diff.d13Cc ~ 0 + FG:Habitat + FG:Habitat:Length.cen +  (1 + Length.cen|drainage/Habitat/FG),data=data[which(data$stream!='Endler' & data$stream!= 'Taylor' & data$stream!= 'Caigual' & data$stream!='El_Cedro' & (data$Habitat=='KG.nat' | data$Habitat=='KGP' | data$Habitat=='KO')),] )
summary(modC.full)
#Singular, but converged


modC <- lmer(diff.d13Cc ~ 0 + FG:Habitat + FG:Habitat:Length.cen +  (1 |drainage/Habitat/FG),data=data[which(data$stream!='Endler' & data$stream!= 'Taylor' & data$stream!= 'Caigual' & data$stream!='El_Cedro' & (data$Habitat=='KG.nat' | data$Habitat=='KGP' | data$Habitat=='KO')),] )
summary(modC)

#LRT test
anova(modC,modC.full,test='LRT')
#LRT significant, so keep full model


r.squaredGLMM(modC.full)

coefC <- summary(modC.full)$coef[,'Estimate']

#natural site int comparison
L <- array(0,c(length(coefC)))
L[which(names(coefC)=='FGGuppy:HabitatKG.nat')] <- 1
L[which(names(coefC)=='FGGuppy:HabitatKGP')] <- -1


#This is the test using lmerTest
contest(modC.full,t(L))
sqrt(contest(modC.full,t(L))$`F value`)

#natural site slope comparison
L <- array(0,c(length(coefC)))
L[which(names(coefC)=='FGGuppy:HabitatKG.nat:Length.cen')] <- 1
L[which(names(coefC)=='FGGuppy:HabitatKGP:Length.cen')] <- -1


#This is the test using lmerTest
contest(modC.full,t(L))
sqrt(contest(modC.full,t(L))$`F value`)



#Killifish
#natural site int comparison
L <- array(0,c(length(coefC),2))
L[which(names(coefC)=='FGKillifish:HabitatKG.nat'),] <- 1
L[which(names(coefC)=='FGKillifish:HabitatKGP'),1] <- -1
L[which(names(coefC)=='FGKillifish:HabitatKO'),2] <- -1

#This is the test using lmerTest
contest(modC.full,t(L))


#natural site slope comparison
L <- array(0,c(length(coefC),2))
L[which(names(coefC)=='FGKillifish:HabitatKG.nat:Length.cen'),] <- 1
L[which(names(coefC)=='FGKillifish:HabitatKGP:Length.cen'),1] <- -1
L[which(names(coefC)=='FGKillifish:HabitatKO:Length.cen'),2] <- -1


#This is the test using lmerTest
contest(modC.full,t(L))












## Plot and Save figure 2
jpeg(file = "Figures/Figure 2.jpeg",width = 9, height = 8, units='in', res=600)

par(mfrow=c(2,2), mar = c(4, 5 ,2, 1), oma = c(0.5, 1, 1, 0.5))

y <- summary(modN)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Estimate']
yplus <- y +  summary(modN)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Std. Error']
yminus <- y -  summary(modN)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Std. Error']

errbar(c(1,3,5)-0.25,y=c(y),yplus=c(yplus),yminus=c(yminus),bty='L',ylab=expression(Intercept~of~paste(Delta^{15}, "N (\u2030)")),
       xlim=c(0.5,5.5),xaxt = "n", xlab="",cex=2,col='orange')
title(main='A)',adj = 0.1, line = 0)
axis(1, at=c(1,3,5), labels=c('','',''),cex.axis=1)
text(x = c(1,3,5),
     y = par("usr")[3] - 0.2,
     labels = c('KGP','KG-Nat','KO'),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1)

y <- summary(modN)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.nat'),'Estimate']
yplus <- y +  summary(modN)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.nat'),'Std. Error']
yminus <- y -  summary(modN)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.nat'),'Std. Error']

errbar(c(1,3)+0.25,y=y,yplus=yplus,yminus=yminus,add=TRUE,pch=17,cex=2)
legend(3,3,c('Killifish','Guppy'),pch=c(16,17),cex=1.5,bty='n',col=c('orange','black'))




y <- summary(modN)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Estimate']
yplus <- y +  summary(modN)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Std. Error']
yminus <- y -  summary(modN)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Std. Error']

errbar(c(1,3,5)-0.25,y=c(y),yplus=c(yplus),yminus=c(yminus),bty='L',ylab=expression(Slope~of~paste(Delta^{15}, "N (\u2030)")),
       xlim=c(0.5,5.5),xaxt = "n", xlab="",cex=2,ylim=c(-0.05,0.125),col='orange')
title(main='B)',adj = 0.1, line = 0)
axis(1, at=c(1,3,5), labels=c('','',''),cex.axis=1)
text(x = c(1,3,5),
     y = par("usr")[3] - 0.02,
     labels = c('KGP','KG-Nat','KO'),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1)

y <- summary(modN)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Estimate']
yplus <- y +  summary(modN)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Std. Error']
yminus <- y -  summary(modN)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Std. Error']

errbar(c(1,3)+0.25,y=y,yplus=yplus,yminus=yminus,add=TRUE,pch=17,cex=2)





y <- summary(modC.full)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Estimate']
yplus <- y +  summary(modC.full)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Std. Error']
yminus <- y -  summary(modC.full)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Std. Error']

errbar(c(1,3,5)-0.25,y=c(y),yplus=c(yplus),yminus=c(yminus),bty='L',ylab=expression(Intercept~of~paste(Delta^{13}, "C (\u2030)")),
       xlim=c(0.5,5.5),xaxt = "n", xlab="Community Type",cex=2,ylim=c(-0.25,5.7),col='orange')
title(main='C)',adj = 0.1, line = 0)
axis(1, at=c(1,3,5), labels=c('','',''),cex.axis=1)
text(x = c(1,3,5),
     y = par("usr")[3] - 0.6,
     labels = c('KGP','KG-Nat','KO'),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1)

y <- summary(modC.full)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.nat'),'Estimate']
yplus <- y +  summary(modC.full)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.nat'),'Std. Error']
yminus <- y -  summary(modC.full)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.nat'),'Std. Error']

errbar(c(1,3)+0.25,y=y,yplus=yplus,yminus=yminus,add=TRUE,pch=17,cex=2)




y <- summary(modC.full)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Estimate']
yplus <- y +  summary(modC.full)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Std. Error']
yminus <- y -  summary(modC.full)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Std. Error']

errbar(c(1,3,5)-0.25,y=c(y),yplus=c(yplus),yminus=c(yminus),bty='L',ylab=expression(Slope~of~paste(Delta^{13}, "C (\u2030)")),
       xlim=c(0.5,5.5),xaxt = "n", xlab="Community Type",cex=2,ylim=c(-0.05,0.075),col='orange')
title(main='D)',adj = 0.1, line = 0)
axis(1, at=c(1,3,5), labels=c('','',''),cex.axis=1)
text(x = c(1,3,5),
     y = par("usr")[3] - 0.012,
     labels = c('KGP','KG-Nat','KO'),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1)

y <- summary(modC.full)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Estimate']
yplus <- y +  summary(modC.full)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Std. Error']
yminus <- y -  summary(modC.full)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Std. Error']

errbar(c(1,3)+0.25,y=y,yplus=yplus,yminus=yminus,add=TRUE,pch=17,cex=2)


dev.off()



#========================================================================================================
#========================================================================================================
###Experimental Communities
#========================================================================================================
#========================================================================================================

#This is the full model with random slopes
modN.full <- lmer(TP ~ 0 + FG:Habitat + FG:Habitat:Length.cen + (1 + Length.cen|drainage/Habitat/FG), data=data )
summary(modN.full)


#Reduced model without random slopes
modN <- lmer(TP ~ 0 + FG:Habitat + FG:Habitat:Length.cen + (1|drainage/Habitat/FG), data=data )
summary(modN)
#Converged

anova(modN,modN.full,test='LRT')
#Reduced model better


coefN <- summary(modN)$coef[,'Estimate']

write.csv(summary(modN)$coef,"./Results/TP Parameters.csv")
write.csv(as.matrix(summary(modN)$vcov),"./Results/d15N vcov.csv")


#Slopes All
L <- array(0,c(length(coefN),3))
L[which(names(coefN)=='FGGuppy:HabitatKGP:Length.cen'),] <- 1
L[which(names(coefN)=='FGGuppy:HabitatKG.nat:Length.cen'),1] <- -1
L[which(names(coefN)=='FGGuppy:HabitatKG.exp:Length.cen'),2] <- -1
L[which(names(coefN)=='FGGuppy:HabitatKG.old:Length.cen'),3] <- -1

contest(modN,t(L))


#ints All
L <- array(0,c(length(coefN),3))
L[which(names(coefN)=='FGGuppy:HabitatKGP'),] <- 1
L[which(names(coefN)=='FGGuppy:HabitatKG.nat'),1] <- -1
L[which(names(coefN)=='FGGuppy:HabitatKG.exp'),2] <- -1
L[which(names(coefN)=='FGGuppy:HabitatKG.old'),3] <- -1

contest(modN,t(L))





#Slopes All
L <- array(0,c(length(coefN),4))
L[which(names(coefN)=='FGKillifish:HabitatKGP:Length.cen'),] <- 1
L[which(names(coefN)=='FGKillifish:HabitatKG.nat:Length.cen'),1] <- -1
L[which(names(coefN)=='FGKillifish:HabitatKO:Length.cen'),2] <- -1
L[which(names(coefN)=='FGKillifish:HabitatKG.exp:Length.cen'),3] <- -1
L[which(names(coefN)=='FGKillifish:HabitatKG.old:Length.cen'),4] <- -1


contest(modN,t(L))


#ints All
L <- array(0,c(length(coefN),4))
L[which(names(coefN)=='FGKillifish:HabitatKGP'),] <- 1
L[which(names(coefN)=='FGKillifish:HabitatKG.nat'),1] <- -1
L[which(names(coefN)=='FGKillifish:HabitatKO'),2] <- -1
L[which(names(coefN)=='FGKillifish:HabitatKG.exp'),3] <- -1
L[which(names(coefN)=='FGKillifish:HabitatKG.old'),4] <- -1


contest(modN,t(L))




#Full model with random slopes
modC.full <- lmer(diff.d13Cc ~ 0 + FG:Habitat + FG:Habitat:Length.cen +  (1+ Length.cen|drainage/Habitat/FG), data=data )
summary(modC.full)
#Failed to converge

#Reduced model without random slopes
modC <- lmer(diff.d13Cc ~ 0 + FG:Habitat + FG:Habitat:Length.cen + (1 |drainage/Habitat/FG), data=data )
summary(modC)


coefC <- summary(modC)$coef[,'Estimate']

write.csv(summary(modC)$coef,"./Results/d13C Parameters.csv")
write.csv(as.matrix(summary(modC)$vcov),"./Results/d13C vcov.csv")

#Slopes All
L <- array(0,c(length(coefC),3))
L[which(names(coefC)=='FGGuppy:HabitatKGP:Length.cen'),] <- 1
L[which(names(coefC)=='FGGuppy:HabitatKG.nat:Length.cen'),1] <- -1
L[which(names(coefC)=='FGGuppy:HabitatKG.exp:Length.cen'),2] <- -1
L[which(names(coefC)=='FGGuppy:HabitatKG.old:Length.cen'),3] <- -1

contest(modC,t(L))


#ints All
L <- array(0,c(length(coefC),3))
L[which(names(coefC)=='FGGuppy:HabitatKGP'),] <- 1
L[which(names(coefC)=='FGGuppy:HabitatKG.nat'),1] <- -1
L[which(names(coefC)=='FGGuppy:HabitatKG.exp'),2] <- -1
L[which(names(coefC)=='FGGuppy:HabitatKG.old'),3] <- -1

contest(modC,t(L))


#Slopes All
L <- array(0,c(length(coefC),4))
L[which(names(coefC)=='FGKillifish:HabitatKGP:Length.cen'),] <- 1
L[which(names(coefC)=='FGKillifish:HabitatKG.nat:Length.cen'),1] <- -1
L[which(names(coefC)=='FGKillifish:HabitatKO:Length.cen'),2] <- -1
L[which(names(coefC)=='FGKillifish:HabitatKG.exp:Length.cen'),3] <- -1
L[which(names(coefC)=='FGKillifish:HabitatKG.old:Length.cen'),4] <- -1


contest(modC,t(L))



#ints All
L <- array(0,c(length(coefC),4))
L[which(names(coefC)=='FGKillifish:HabitatKGP'),] <- 1
L[which(names(coefC)=='FGKillifish:HabitatKG.nat'),1] <- -1
L[which(names(coefC)=='FGKillifish:HabitatKO'),2] <- -1
L[which(names(coefC)=='FGKillifish:HabitatKG.exp'),3] <- -1
L[which(names(coefC)=='FGKillifish:HabitatKG.old'),4] <- -1


contest(modC,t(L))
















jpeg(file = "Figures/Figure 3.jpeg",width = 9, height = 8, units='in', res=600)

par(mfrow=c(2,2), mar = c(4, 5 ,2, 1), oma = c(0.5, 1, 1, 0.5))

y <- summary(modN)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.exp','FGKillifish:HabitatKG.old','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Estimate']
yplus <- y +  summary(modN)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.exp','FGKillifish:HabitatKG.old','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Std. Error']
yminus <- y -  summary(modN)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.exp','FGKillifish:HabitatKG.old','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Std. Error']

errbar(c(1,3,5,7,9)-0.25,y=c(y),yplus=c(yplus),yminus=c(yminus),bty='L',ylab=expression(Intercept~of~paste(Delta^{15}, "N (\u2030)")),
       xlim=c(0.5,9.5),xaxt = "n", xlab="",cex=2,col='orange',ylim=c(2,5))
title(main='A)',adj = 0.1, line = 0)
axis(1, at=c(1,3,5,7,9), labels=c('','','','',''),cex.axis=1)
text(x = c(1,3,5,7,9),
     y = par("usr")[3] - 0.35,
     labels = c('KGP','KG-New','KG-Old','KG-Nat','KO'),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1)

y <- summary(modN)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.exp','FGGuppy:HabitatKG.old','FGGuppy:HabitatKG.nat'),'Estimate']
yplus <- y +  summary(modN)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.exp','FGGuppy:HabitatKG.old','FGGuppy:HabitatKG.nat'),'Std. Error']
yminus <- y -  summary(modN)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.exp','FGGuppy:HabitatKG.old','FGGuppy:HabitatKG.nat'),'Std. Error']

errbar(c(1,3,5,7)+0.25,y=y,yplus=yplus,yminus=yminus,add=TRUE,pch=17,cex=2)
legend(6.75,5,c('Killifish','Guppy'),pch=c(16,17),cex=1.5,bty='n',col=c('orange','black'))




y <- summary(modN)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.exp:Length.cen','FGKillifish:HabitatKG.old:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Estimate']
yplus <- y +  summary(modN)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.exp:Length.cen','FGKillifish:HabitatKG.old:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Std. Error']
yminus <- y -  summary(modN)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.exp:Length.cen','FGKillifish:HabitatKG.old:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Std. Error']

errbar(c(1,3,5,7,9)-0.25,y=c(y),yplus=c(yplus),yminus=c(yminus),bty='L',ylab=expression(Slope~of~paste(Delta^{15}, "N (\u2030)")),
       xlim=c(0.5,9.5),xaxt = "n", xlab="",cex=2,ylim=c(-0.05,0.125),col='orange')
title(main='B)',adj = 0.1, line = 0)
axis(1, at=c(1,3,5,7,9), labels=c('','','','',''),cex.axis=1)
text(x = c(1,3,5,7,9),
     y = par("usr")[3] - 0.02,
     labels = c('KGP','KG-New','KG-Old','KG-Nat','KO'),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1)

y <- summary(modN)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.exp:Length.cen','FGGuppy:HabitatKG.old:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Estimate']
yplus <- y +  summary(modN)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.exp:Length.cen','FGGuppy:HabitatKG.old:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Std. Error']
yminus <- y -  summary(modN)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.exp:Length.cen','FGGuppy:HabitatKG.old:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Std. Error']

errbar(c(1,3,5,7)+0.25,y=y,yplus=yplus,yminus=yminus,add=TRUE,pch=17,cex=2)





y <- summary(modC)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.exp','FGKillifish:HabitatKG.old','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Estimate']
yplus <- y +  summary(modC)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.exp','FGKillifish:HabitatKG.old','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Std. Error']
yminus <- y -  summary(modC)$coef[c('FGKillifish:HabitatKGP','FGKillifish:HabitatKG.exp','FGKillifish:HabitatKG.old','FGKillifish:HabitatKG.nat','FGKillifish:HabitatKO'),'Std. Error']

errbar(c(1,3,5,7,9)-0.25,y=c(y),yplus=c(yplus),yminus=c(yminus),bty='L',ylab=expression(Intercept~of~paste(Delta^{13}, "C (\u2030)")),
       xlim=c(0.5,9.5),xaxt = "n", xlab="",cex=2,col='orange',ylim=c(0,5))
title(main='C)',adj = 0.1, line = 0)
axis(1, at=c(1,3,5,7,9), labels=c('','','','',''),cex.axis=1)
text(x = c(1,3,5,7,9),
     y = par("usr")[3] - 0.7,
     labels = c('KGP','KG-New','KG-Old','KG-Nat','KO'),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1)

y <- summary(modC)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.exp','FGGuppy:HabitatKG.old','FGGuppy:HabitatKG.nat'),'Estimate']
yplus <- y +  summary(modC)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.exp','FGGuppy:HabitatKG.old','FGGuppy:HabitatKG.nat'),'Std. Error']
yminus <- y -  summary(modC)$coef[c('FGGuppy:HabitatKGP','FGGuppy:HabitatKG.exp','FGGuppy:HabitatKG.old','FGGuppy:HabitatKG.nat'),'Std. Error']

errbar(c(1,3,5,7)+0.25,y=y,yplus=yplus,yminus=yminus,add=TRUE,pch=17,cex=2)
# legend(7,3,c('Killifish','Guppy'),pch=c(16,17),cex=1.5,bty='n',col=c('orange','black'))




y <- summary(modC)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.exp:Length.cen','FGKillifish:HabitatKG.old:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Estimate']
yplus <- y +  summary(modC)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.exp:Length.cen','FGKillifish:HabitatKG.old:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Std. Error']
yminus <- y -  summary(modC)$coef[c('FGKillifish:HabitatKGP:Length.cen','FGKillifish:HabitatKG.exp:Length.cen','FGKillifish:HabitatKG.old:Length.cen','FGKillifish:HabitatKG.nat:Length.cen','FGKillifish:HabitatKO:Length.cen'),'Std. Error']

errbar(c(1,3,5,7,9)-0.25,y=c(y),yplus=c(yplus),yminus=c(yminus),bty='L',ylab=expression(Slope~of~paste(Delta^{13}, "C (\u2030)")),
       xlim=c(0.5,9.5),xaxt = "n", xlab="",cex=2,ylim=c(-0.12,0.07),col='orange')
title(main='D)',adj = 0.1, line = 0)
axis(1, at=c(1,3,5,7,9), labels=c('','','','',''),cex.axis=1)
text(x = c(1,3,5,7,9),
     y = par("usr")[3] - 0.02,
     labels = c('KGP','KG-New','KG-Old','KG-Nat','KO'),
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1)

y <- summary(modC)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.exp:Length.cen','FGGuppy:HabitatKG.old:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Estimate']
yplus <- y +  summary(modC)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.exp:Length.cen','FGGuppy:HabitatKG.old:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Std. Error']
yminus <- y -  summary(modC)$coef[c('FGGuppy:HabitatKGP:Length.cen','FGGuppy:HabitatKG.exp:Length.cen','FGGuppy:HabitatKG.old:Length.cen','FGGuppy:HabitatKG.nat:Length.cen'),'Std. Error']

errbar(c(1,3,5,7)+0.25,y=y,yplus=yplus,yminus=yminus,add=TRUE,pch=17,cex=2)

dev.off()


#Supplemental Figures
source('R/2 d15N Figures New.R')
source('R/2 d13C Figures New.R')



