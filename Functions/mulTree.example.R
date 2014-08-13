#Example for using the functions in mulTree package
#Kind of how we got the data for the paper

#install the packages
#install.packages(c("ape","hdrcde","caper","coda","MCMCglmm"))
library(ape)
library(hdrcde)
library(caper)
library(coda)
library(MCMCglmm)

#load the functions
source("rTreeBind.R")
source("as.mulTree.R")
source("mulTree.R")
source("read.mulTree.R")
source("hdr.mulTree.R")
source("sum.mulTree.R")

#rTreeBind
#randomly binding the trees together

#loading the data
mammal_trees<-read.tree("example/10mammaltrees_329taxa.tre")
bird_trees<-read.tree("example/10avetrees_96taxa.tre")

#running the function
combined_trees<-rTreeBind(x=mammal_trees, y=bird_trees, sample=10, root.age=250)
#with
#the two first arguments being the two trees to binding
#the third being how many trees to we want in the end
#the forth being where to bind them (what age?)

#as.mulTree
#creating the 'mulTree' object prior to the analysis to be input in the main function

#loading the data
longevity_data<-read.csv("example/425sp_LongevityData.csv",header=TRUE)

#running the function
mulTree_data<-as.mulTree(data=longevity_data, trees=combined_trees, species="species.m")
#with:
#the first argument being the data to analyse
#the second argument being a multiPhylo object
#the third argument being the name of the column containing the species

#mulTree
#formula
long.formula<-long.m~BMR + mass.m + volant.m + mass.m:volant.m
#mcmc parameters (number of generations, thin/sampling, burnin)
mcmc.parameters<-c(10000, 500, 2500)
#priors (R-structure, G-structure)
mcmc.priors<-list(R = list(V = 1/2, nu=0.002), G = list(G1=list(V = 1/2, nu=0.002)))

#performing MCMCglmm on multiple trees
mulTree(mulTree_data, formula=long.formula, parameters=mcmc.parameters, priors=mcmc.priors, verbose=TRUE, output="longevity.example", warn=FALSE)
#with
#warn=FALSE disabling the warnings from the MCMCglmm function

#Obviously way more generations are needed to converge and to obtain a good ESS.
#But this example already takes 2 minutes
#The function saved all the models and the convergence diagnosis per tree using the "longevity.example" chain name in the current directory

#read.mulTree
#Checking the convergence diagnosis for the trees
read.mulTree(mulTree.mcmc="longevity.example", convergence=TRUE)

#Reading all the chains
mulTree_mcmc<-read.mulTree(mulTree.mcmc="longevity.example")

#sum.mulTree
#summarizing the results of all the chains and all the trees
sum.mulTree(mulTree_mcmc)

#