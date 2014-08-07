##########################
#Run MCMCglmm on a 'mulTree' object
##########################
#Running a MCMCglmm model on a list of phylogenies and the data stored in a 'mulTree' object. The results can be written out of R environment as individual models.
#v0.1
##########################
#SYNTAX :
#<mulTree> a 'mulTree' object obtained from as.mulTree function
#<formula> an object of class 'formula'
#<priors> a series of priors to use for the mcmc
#<parameters> a series of parameters to use for the mcmc
#<verbose> whether to be verbose or not (default=FALSE)
#<output> any optional string of characters that will be used as chain name for the models output (default=FALSE)
##########################
#----
#guillert(at)tcd.ie & healyke(at)tcd.ie - 07/08/2014
##########################
#Requirements:
#-R 3
#-R package "ape"
#-R package "caper"
#-R package "mvtnorm"
#-R package "MCMCglmm"
#-R package "modeest"
#-R package "coda"
##########################


M3mcmcglmm<-function(mulTree, formula, priors, parameters, verbose=FALSE, output=FALSE)
{stop("IN DEVELOPEMENT")
#HEADER
    #libraries
    require(ape)
    require(caper)
    require(mvtnorm)
    require(coda)

#DATA
    #mulTree
    if(class(mulTree) != 'mulTree') {
        stop(as.character(substitute(mulTree))," is not a \"mulTree\" object. \nUse as.mulTree() function.", call.=FALSE)
    } else {
        if(length(mulTree) != 3) {
            stop(as.character(substitute(mulTree))," is not a \"mulTree\" object. \nUse as.mulTree() function.", call.=FALSE)
        } else {
            if(class(mulTree[[1]]) != 'multiPhylo') {
                stop(as.character(substitute(mulTree))," is not a \"mulTree\" object. \nUse as.mulTree() function.", call.=FALSE)
            } else {
                if(class(mulTree[[2]]) != 'data.frame') {
                    stop(as.character(substitute(mulTree))," is not a \"mulTree\" object. \nUse as.mulTree() function.", call.=FALSE)
                }
            }
        }
    }

#FUNCTION


#RUNING THE MODELS


#OUTPUT


#Step 3: RUNNING THE MODELS
        
    #running the models (overwriting)

        #Bayesian framework
            #set up uniformative prior. see Jarrod Hadfield's notes on prior parametrers used. 
        else { prior<-list(R = list(V = 1/2, nu=0.002), G = list(G1=list(V = 1/2, nu=0.002)))
            for (i in (1:length(trees))){

                #Model running using MCMCglmm function (MCMCglmm) on each tree [i] on two independent chains

                    #Chain 1
                    model<-MCMCglmm(formula, random=~animal,pedigree=trees[[i]],prior=prior,data=comparative.data(data=data.frame, phy=trees[[i]], names.col="sp.col", vcv=FALSE)$data,verbose=FALSE,family=c("gaussian"),nitt=ngen,burnin=as.integer(ngen/6),thin=thinn)
                
                    #Chain 2
                    model.1<-MCMCglmm(formula, random=~animal,pedigree=trees[[i]],prior=prior,data=comparative.data(data=data.frame, phy=trees[[i]], names.col="sp.col", vcv=FALSE)$data,verbose=FALSE,family=c("gaussian"),nitt=ngen/50,burnin=as.integer(ngen/300),thin=thinn/50)
                    
                    #Convergence check using Gelman and Rubins diagnoses set to return true or false based on level of scale reduction set (default == 1.1)
                    convergence<-gelman.diag(mcmc.list(as.mcmc(model.1$Sol[1:(length(model.1$Sol[,1])),]),as.mcmc(model$Sol[1:(length(model.1$Sol[,1])),])))

                #Saving the model ran on each tree [i]
                save(model, file=as.character(file.names[[i]]))

                #Printing the time and the ESS + convergence diagnosis
                cat(format(Sys.time(), "%H:%M:%S"), "-", output, "on tree", as.character(i), "done", "\n")
                cat(format(Sys.time(), "%H:%M:%S"), "-", "Effective sample size is >1000:",all(effectiveSize(model$Sol[])>1000), "\n")
                cat(format(Sys.time(), "%H:%M:%S"), "-", "All levels converged:",all(convergence$psrf[,1]<converge), "\n")
            }
        }

#Step 4: OUTPUT                     
    
    #Output
    cat(format(Sys.time(), "%H:%M:%S"), "-", output, "run and saved on all trees","\n")
    cat("Use SumMultiPGLS function to summarise the output", "\n")

#End
}