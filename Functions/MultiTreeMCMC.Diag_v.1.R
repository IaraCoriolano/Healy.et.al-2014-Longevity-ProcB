##########################
#Summary MultiTreeMCMCglmm output
##########################
#Extract summary from MCMCglmm models
#v1.1
#----
#guillert(at)tcd.ie & healyke(at)tcd.ie - 3/03/2012
##########################
#After running SumMultiPGLS script

#Creat object for the model output
mods <- SumMultiPGLS("Model1_tree_")

#Extract names for both fixed and random terms
Sol.names <-names(mods$Sol[[1]][1,])
VCV.names <- c("Residuals","Phylo.term")

#creat a single names object
names.mod <- c(Sol.names,VCV.names)

#matrix to store outputs
output.est <- matrix(0,nrow = (length(mods$Sol[]))*(length(mods$Sol[[1]][,1]))
, ncol = length(names.mod), dimnames = list(c(),c(names.mod)))

#loops through all fixed and random terms creating a single matrix of the terms for the postirior outputs form all the tree
 
 for(z in 1:(mods$chain.length)){
 	for(t in 1:(length(Sol.names[]))){
		#the fixed factors
		output.est[((z*(length(mods$Sol[[1]][,1])))-(length(mods$Sol[[1]][,1]))+1):((length(mods$Sol[[1]][,1]))*z),t] <- c(as.mcmc(mods$Sol[[z]][,t]))
		
		
		#the residual term
		output.est[((z*(length(mods$Sol[[1]][,1])))-(length(mods$Sol[[1]][,1]))+1):((length(mods$Sol[[1]][,1]))*z),length(Sol.names) + 1] <- c(as.mcmc(mods$Residuals[[z]]))
		
		#the phylo terms (animal)
		output.est[((z*(length(mods$Sol[[1]][,1])))-(length(mods$Sol[[1]][,1]))+1):((length(mods$Sol[[1]][,1]))*z),length(Sol.names) + 2] <- c(as.mcmc(mods$Phylo.term[[z]]))

		
		}
}




##########################################################################################
#######################------------out put duaignostics-----------########################
##########################################################################################

#make mcmc for ploting functions
output.est <- as.mcmc(output.est)

trace(output.est[,1])

########------------calculate mode and CI---------------################
#create a little table of the summary values
sum.val <- matrix(0,ncol = 3
, nrow = length(names.mod), dimnames = list(c(names.mod),c("mode","C.I.lower","C.I.upper")))


for(i in 1:(length(names.mod))){
#mode
sum.val[i,1]    <-  as.numeric(mlv(output.est[,i],method="mfv")[1])
#C.I
sum.val[i,2:3] <- as.vector(quantile(output.est[,i],c(0.025,0.975)))

}

# P>0
#sum(mod.full$mc[,1]>=0)/length(m1$Sol[,1])


#to plot the correlations between each of the parameters in the models in a panel 
#this plot is in the r scripts and also the stability folders
#source("panel_functions.r") 

#pairs(output.est[,1:13],
#	diag.panel=panel.hist,
#	lower.panel=panel.cor,
#	upper.panel=panel.linear)
#t0 plot posterior distribution use c(as.mcmc(x1),as.mcmc(x2),etc)

#to cheack for convergane, the seperate chains are one vector now
#gelman.plot(mcmc.list(as.mcmc(output.est[1:5000]),as.mcmc(output.est[5001:10000])))

#to cheack for temporlal autocorrelation
#acf(as.numeric(output.est[,1]),lag.max=100)





#########to extract DIC###########

DIC.val <- matrix(0,ncol = 1
, nrow = length(names.mod), dimnames = list(c(names.mod),c("DIC")))


for(i in 1:(length(names.mod))){
#mode
DIC.val[i] <- mods$DIC[[i]]
DIC.m<- mlv(DIC.val[,i],method="mfv")
DIC.v <- sum(DIC.val)/(length(DIC.val))
}

##########if comparing both aves and mammals##############

mods.a <- mods #whatever the first loaded models are
mods.b <- SumMultiPGLS("Model1_tree_") #either the aves or mam models

DIC.com <- matrix(0,ncol = 1
, nrow = mods$chain.length, dimnames = list(seq(1,mods$chain.length),c("DIC.mode")))


for(i in 1:mods$chain.length){
#mode
DIC.com[i] <- mods.a$DIC[[i]] + mods.b$DIC[[i]]
DIC.c<- mlv(DIC.com[],method="mfv")
DIC.ci <- quantile(DIC.com[],c(0.025,0.975))
}


#if loading them all up and comparing
mods.v <- mods #whatever the first loaded models are
mods.nv <- SumMultiPGLS("Model1_tree_") #either the aves or mam models

mods.aves <-SumMultiPGLS("Model1_tree_")
mods.mammals <-  SumMultiPGLS("Model1_tree_")


DIC.main <- matrix(0,ncol = 1
, nrow = mods$chain.length, dimnames = list(seq(1,mods$chain.length),c("DIC.mode")))
DIC.av.ma <- matrix(0,ncol = 1
, nrow = mods$chain.length, dimnames = list(seq(1,mods$chain.length),c("DIC.mode")))
DIC.diff <- matrix(0,ncol = 1
, nrow = mods$chain.length, dimnames = list(seq(1,mods$chain.length),c("DIC.mode")))


for(i in 1:mods$chain.length){
#mode
DIC.main[i] <- mods.v$DIC[[i]] + mods.nv$DIC[[i]]
DIC.av.ma[i] <- mods.aves$DIC[[i]] + mods.mammals$DIC[[i]]
DIC.diff[i] <- DIC.main[i] - DIC.av.ma[i]
DIC.d<- mlv(DIC.diff[],method="mfv")
DIC.di <- quantile(DIC.diff[],c(0.025,0.975))
}



##########################################################################################
################---to plot the output as a density barplpot-----------####################
##########################################################################################

#source(densityplot.r)
densityplot(output.est)
abline(0,0,col="grey",lty= 2)


