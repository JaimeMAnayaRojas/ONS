jpeg(file = "./Figures/Figure-3.jpeg",width = 7.5, height = 6, units='in', res=600)
par(mfrow=c(2,3), mar = c(4, 4 ,2, 1), oma = c(0.5, 1, 1, 0.5))

toplot <- data[which(data$Habitat=='KGP'),]

min <- -2.5
max <- 10


toplot <- merge(x=toplot,y=rand,by='drainage',all.x=TRUE)
toplot$d13Ccc <- toplot$diff.d13Cc - toplot$randC
plot(toplot[which(toplot$FG=='Killifish'),'Length'],toplot[which(toplot$FG=='Killifish'),'d13Ccc'],
     ylim=c(min,max),xlim=c(min(toplot$Length,na.rm=TRUE),max(toplot$Length,na.rm=TRUE)),pch=16,col='orange',bty='L',
     ylab=expression(paste(Delta^{13}, "C (\u2030)")),xlab='Fish Standard Length (mm)')
points(toplot[which(toplot$FG=='Guppy'),'Length'],toplot[which(toplot$FG=='Guppy'),'d13Ccc'],pch=16,col='grey')

gsizes <- seq(10,max(toplot[which(toplot$FG=='Guppy'),'Length'],na.rm=TRUE),by=1)
ksizes <- seq(10,max(toplot[which(toplot$FG=='Killifish'),'Length'],na.rm=TRUE),by=1)
lines(ksizes,
      (summary(mod)$coef['FGKillifish:HabitatKGP','Estimate'] + summary(mod)$coef['FGKillifish:HabitatKGP:Length.cen','Estimate']*(ksizes-center.size)),
      lwd=3,col='orange')

lines(gsizes,
      (summary(mod)$coef['FGGuppy:HabitatKGP','Estimate'] + summary(mod)$coef['FGGuppy:HabitatKGP:Length.cen','Estimate']*(gsizes-center.size)),
      lwd=3,col='black')
title(main='A) KGP-Nat',adj = 0.1, line = 0)
legend('topleft',c('Guppy','Killifish'),lty=c(1,1),pch=c(16,16),col=c('grey','orange'),bty='n')





toplot <- data[which(data$Habitat=='KG.exp'),]
# 
toplot <- merge(x=toplot,y=rand,by='drainage',all.x=TRUE)
toplot$d13Ccc <- toplot$diff.d13Cc - toplot$randC
plot(toplot[which(toplot$FG=='Killifish'),'Length'],toplot[which(toplot$FG=='Killifish'),'d13Ccc'],
     ylim=c(min,max),xlim=c(min(toplot$Length,na.rm=TRUE),max(toplot$Length,na.rm=TRUE)),pch=16,col='orange',bty='L',
     ylab=expression(paste(Delta^{13}, "C (\u2030)")),xlab='Fish Standard Length (mm)')
points(toplot[which(toplot$FG=='Guppy'),'Length'],toplot[which(toplot$FG=='Guppy'),'d13Ccc'],pch=16,col='grey')

gsizes <- seq(10,max(toplot[which(toplot$FG=='Guppy'),'Length'],na.rm=TRUE),by=1)
ksizes <- seq(10,max(toplot[which(toplot$FG=='Killifish'),'Length'],na.rm=TRUE),by=1)
lines(ksizes,
      (summary(mod)$coef['FGKillifish:HabitatKG.exp','Estimate'] + summary(mod)$coef['FGKillifish:HabitatKG.exp:Length.cen','Estimate']*(ksizes-center.size)),
      lwd=3,col='orange')

lines(gsizes,
      (summary(mod)$coef['FGGuppy:HabitatKG.exp','Estimate'] + summary(mod)$coef['FGGuppy:HabitatKG.exp:Length.cen','Estimate']*(gsizes-center.size)),
      lwd=3,col='black')
title(main='B) KG-New exp',adj = 0.1, line = 0)


toplot <- data[which(data$Habitat=='KG.old'),]
# 
toplot <- merge(x=toplot,y=rand,by='drainage',all.x=TRUE)
toplot$d13Ccc <- toplot$diff.d13Cc - toplot$randC
plot(toplot[which(toplot$FG=='Killifish'),'Length'],toplot[which(toplot$FG=='Killifish'),'d13Ccc'],
     ylim=c(min,max),xlim=c(min(toplot$Length,na.rm=TRUE),max(toplot$Length,na.rm=TRUE)),pch=16,col='orange',bty='L',
     ylab=expression(paste(Delta^{13}, "C (\u2030)")),xlab='Fish Standard Length (mm)')
points(toplot[which(toplot$FG=='Guppy'),'Length'],toplot[which(toplot$FG=='Guppy'),'d13Ccc'],pch=16,col='grey')

gsizes <- seq(10,max(toplot[which(toplot$FG=='Guppy'),'Length'],na.rm=TRUE),by=1)
ksizes <- seq(10,max(toplot[which(toplot$FG=='Killifish'),'Length'],na.rm=TRUE),by=1)
lines(ksizes,
      (summary(mod)$coef['FGKillifish:HabitatKG.old','Estimate'] + summary(mod)$coef['FGKillifish:HabitatKG.old:Length.cen','Estimate']*(ksizes-center.size)),
      lwd=3,col='orange')

lines(gsizes,
      (summary(mod)$coef['FGGuppy:HabitatKG.old','Estimate'] + summary(mod)$coef['FGGuppy:HabitatKG.old:Length.cen','Estimate']*(gsizes-center.size)),
      lwd=3,col='black')
title(main='C) KG-Old exp',adj = 0.1, line = 0)



toplot <- data[which(data$Habitat=='KG.nat'),]
# 
toplot <- merge(x=toplot,y=rand,by='drainage',all.x=TRUE)
toplot$d13Ccc <- toplot$diff.d13Cc - toplot$randC
plot(toplot[which(toplot$FG=='Killifish'),'Length'],toplot[which(toplot$FG=='Killifish'),'d13Ccc'],
     ylim=c(min,max),xlim=c(min(toplot$Length,na.rm=TRUE),max(toplot$Length,na.rm=TRUE)),pch=16,col='orange',bty='L',
     ylab=expression(paste(Delta^{13}, "C (\u2030)")),xlab='Fish Standard Length (mm)')
points(toplot[which(toplot$FG=='Guppy'),'Length'],toplot[which(toplot$FG=='Guppy'),'d13Ccc'],pch=16,col='grey')

gsizes <- seq(10,max(toplot[which(toplot$FG=='Guppy'),'Length'],na.rm=TRUE),by=1)
ksizes <- seq(10,max(toplot[which(toplot$FG=='Killifish'),'Length'],na.rm=TRUE),by=1)
lines(ksizes,
      (summary(mod)$coef['FGKillifish:HabitatKG.nat','Estimate'] + summary(mod)$coef['FGKillifish:HabitatKG.nat:Length.cen','Estimate']*(ksizes-center.size)),
      lwd=3,col='orange')

lines(gsizes,
      (summary(mod)$coef['FGGuppy:HabitatKG.nat','Estimate'] + summary(mod)$coef['FGGuppy:HabitatKG.nat:Length.cen','Estimate']*(gsizes-center.size)),
      lwd=3,col='black')
title(main='D) KG-Nat',adj = 0.1, line = 0)




toplot <- data[which(data$Habitat=='KO'),]
# 
toplot <- merge(x=toplot,y=rand,by='drainage',all.x=TRUE)
toplot$d13Ccc <- toplot$diff.d13Cc - toplot$randC
plot(toplot[which(toplot$FG=='Killifish'),'Length'],toplot[which(toplot$FG=='Killifish'),'d13Ccc'],
     ylim=c(min,max),xlim=c(min(toplot$Length,na.rm=TRUE),max(toplot$Length,na.rm=TRUE)),pch=16,col='orange',bty='L',
     ylab=expression(paste(Delta^{13}, "C (\u2030)")),xlab='Fish Standard Length (mm)')

ksizes <- seq(10,max(toplot[which(toplot$FG=='Killifish'),'Length'],na.rm=TRUE),by=1)
lines(ksizes,
      (summary(mod)$coef['FGKillifish:HabitatKO','Estimate'] + summary(mod)$coef['FGKillifish:HabitatKO:Length.cen','Estimate']*(ksizes-center.size)),
      lwd=3,col='orange')
title(main='E) KO-Nat',adj = 0.1, line = 0)




graphics.off()



