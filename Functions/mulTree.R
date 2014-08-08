##########################
#Run MCMCglmm on a 'mulTree' object
##########################
#Running a MCMCglmm model on a list of phylogenies and the data stored in a 'mulTree' object. The results can be written out of R environment as individual models.
#v0.1
##########################
#SYNTAX :
#<mulTree> a 'mulTree' object obtained from as.mulTree function
#<formula> an object of class 'formula'
#<chains> the number of independent chains for the mcmc
#<parameters> a vector of three elements to use as parameters for the mcmc. Should be respectively Number of generations, sampling and burnin.
#<priors> a series of priors to use for the mcmc (default=NULL is using the default parameters from MCMCglmm function)
#<...> any additional argument to be passed to MCMCglmm() function
#<convergence> a numerical value for assessing chain convergence (default=1.1)
#<ESS> a numerical value for assessing the effective sample size (default=1000)
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


mcmcglmm.mulTree<-function(mulTree, formula, parameters, chains=2, priors=NULL, ...,convergence=1.1, ESS=1000, verbose=FALSE, output=FALSE)
{stop("IN DEVELOPEMENT")
#HEADER
    #libraries
    require(ape)
    require(caper)
    require(mvtnorm)
    require(coda)
    require(modeest)

#DATA
    #mulTree
    if(class(mulTree) != 'mulTree') {
        stop(as.character(substitute(mulTree))," is not a \"mulTree\" object.\nUse as.mulTree() function.", call.=FALSE)
    } else {
        if(length(mulTree) != 3) {
            stop(as.character(substitute(mulTree))," is not a \"mulTree\" object.\nUse as.mulTree() function.", call.=FALSE)
        } else {
            if(class(mulTree[[1]]) != 'multiPhylo') {
                stop(as.character(substitute(mulTree))," is not a \"mulTree\" object.\nUse as.mulTree() function.", call.=FALSE)
            } else {
                if(class(mulTree[[2]]) != 'data.frame') {
                    stop(as.character(substitute(mulTree))," is not a \"mulTree\" object.\nUse as.mulTree() function.", call.=FALSE)
                }
            }
        }
    }

    #formula
    if(class(formula) != 'formula') {
        stop(as.character(substitute(formula))," is not a \"formula\" object.", call.=FALSE)
    }

    #chains
    if(class(chains) != 'numeric') {
        stop("\"chains\" must be numeric.", call.=FALSE) 
    } else {
        if(length(chains) != 1) {
            stop("\"chains\" must be numeric.", call.=FALSE) 
        }
    }

    #parameters
    if(class(parameters) != 'numeric') {
        stop(as.character(substitute(parameters))," is not a \"vector\" object.", call.=FALSE)
    } else {
        if(length(parameters) != 3) {
            stop("Wrong format for ",as.character(substitute(parameters)),", must be a vector of three elements:\nthe number of generations ; the sampling and the burnin.", call.=FALSE) 
        }
    }

    #priors
    if(class(priors) == NULL) {
        prior.default=TRUE
    } else {
        prior.default=FALSE
        if(class(priors) != 'list') {
            stop("Wrong format for ",as.character(substitute(priors)),", must be a list of three elements:\nsee ?MCMCglmm manual.", call.=FALSE) 
        }
    }

    #convergence
    if(class(convergence) != 'numeric') {
        stop("\"convergence\" must be numeric.", call.=FALSE) 
    } else {
        if(length(convergence) != 1) {
            stop("\"convergence\" must be numeric.", call.=FALSE) 
        }
    }

    #ESS
    if(class(ESS) != 'numeric') {
        stop("\"ESS\" must be numeric.", call.=FALSE) 
    } else {
        if(length(ESS) != 1) {
            stop("\"ESS\" must be numeric.", call.=FALSE) 
        }
    }

    #verbose
    if(class(verbose) != 'logical') {
        stop("\"verbose\" must be logical.", call.=FALSE) 
    }

    #output
    if(class(output) == 'logical') {
        if(output==FALSE){
            do.output=FALSE
            cat("No output option selected, this might highly decrease the performances of this function.")
        } else {
            do.output=TRUE
            output="mulTree.out"
            cat("Analysis output is set to \"mulTree.out\" by default.\nTo modify it, specify the output chain name using:\nmulTree(..., output=<OUTPUT_NAME>, ...)") 
        }
    }
    if(class(output) != 'logical') {
        if(class(output) != 'character') {
            stop(as.character(substitute(output)),", must be a chain of characters", call.=FALSE) 
        } else {
            if(length(output) != 1) {
                stop(as.character(substitute(output)),", must be a chain of characters", call.=FALSE) 
            } else {
                do.output=TRUE
            }
        }
    }
#FUNCTION


    FUN.MCMCglmm<-function(ntree, mulTree, formula, priors, parameters, ...){
        #Model running using MCMCglmm function (MCMCglmm) on each tree [i] on two independent chains
        model<-MCMCglmm(formula, random=~animal, pedigree=mulTree$phy[[ntree]], prior=priors, data=mulTree$data, verbose=FALSE, nitt=parameters[1], thin=parameters[2], burnin=parameters[3], ...)
        return(model)
    }

    FUN.convergence.test<-function(chains){
        #Creating the mcmc.list
        list.mcmc<-list()
        for (nchains in 1:chains){
            list.mcmc[[nchains]]<-as.mcmc(get(paste("model_chain",nchains,sep=""))$Sol[1:(length(get(paste("model_chain",nchains,sep=""))$Sol[,1])),])
        }

        #Convergence check using Gelman and Rubins diagnoses set to return true or false based on level of scale reduction set (default = 1.1)
        convergence<-gelman.diag(mcmc.list(list.mcmc))
        return(convergence)
    }


#RUNNING THE MODELS

    #Creating the list of file names per tree
    file.names<-vector("list", length(mulTree$phy))
    for (ntree in (1:length(mulTree$phy))){
        file.names[[ntree]]<-paste(output, as.character("_tree"), as.character(ntree), sep="")
    }
    #Adding the chain name
    for (nchains in 1:chains) {
        file.names.chain<-paste(file.names, as.character("_chain"), as.character(nchains), as.character(".R"), sep="")
        assign(paste("file.names.chain", nchains, sep=""), file.names.chain)
    }

    #Running the models n times for every trees
    for (ntree in 1:length(mulTree$phy)) {
        
        #Running the model for one tree
        for (nchains in 1:chains) {
            model_chain<-FUN.MCMCglmm(ntree, mulTree, formula, priors, parameters, ...)
            assign(paste("model_chain", nchains, sep=""), model_chain)
        }

        #Testing the convergence for one tree
        converge.test<-FUN.convergence.test(chains)

        #Saving all the chains for one tree
        #output?
        
        #verbose?

    }

    #Output
    cat(format(Sys.time(), "%H:%M:%S"), "-", output, "run and saved on all trees","\n")
    cat("Use SumMultiPGLS function to summarise the output", "\n")


for(i in 2:5) {
  thismodel <- get(paste("model", i, sep=""))
  rsq <- c(rsq, summary(thismodel)$r.squared)
} 



x=1

while (x<100)

 {
   vectorx<- rnorm(100)
    assign(paste("vector",x,sep=""),vectorx)
x=x+1

} 


for(i in 1:6) { #-- Create objects  'r.1', 'r.2', ... 'r.6' --
    nam <- paste("r", i, sep = ".")
    assign(nam, 1:i)
}



    for (ntree in 1:length(mulTree$phy)) {
        
        for (nchain in 1:chains) {
            model_chain${nchain}<-FUN.MCMCglmm(ntree, mulTree, formula, priors, parameters, ...)
        }

    }

save all chains

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