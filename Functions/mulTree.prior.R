##########################
#Generates the list of priors to use for the MCMCglmmm in mulTree()
##########################
#Generates the proper format list of priors to be used in the MCMCglmm function called in mulTree()
#v0.1
##########################
#SYNTAX :
#<R> R-structure : a list with the expected (co)variances (V) and degree of belief parameter (nu) for the inverse-Wishart, and also the mean vector (alpha.mu) and covariance matrix (alpha.V) for the redundant working parameters.
#<G> G-structure : a list with the expected (co)variances (V) and degree of belief parameter (nu) for the inverse-Wishart, and also the mean vector (alpha.mu) and covariance matrix (alpha.V) for the redundant working parameters.
#<B> fixed effects: a list with the expected value (mu) and a (co)variance matrix (V) representing the strength of belief
##########################
#----
#guillert(at)tcd.ie & healyke(at)tcd.ie - 08/08/2014
##########################
#Requirements:
#-R 3
##########################


mulTree.priors<-function(R=NULL, G=NULL, B=NULL) {
#DATA
    #R
    if(class(R) != 'list') {
        stop(as.character(substitute(R))," must be a list with the expected (co)variance (V) and the degree of belief (nu).", call.=FALSE)
    }

    #G
    if(class(G) != 'list') {
        stop(as.character(substitute(G))," must be a list with the expected (co)variance (V) and the degree of belief (nu).", call.=FALSE)
    }

    #B
    if(class(B) != 'list') {
        stop(as.character(substitute(B))," must be a list with the expected (co)variance (V) and the degree of belief (nu).", call.=FALSE)
    }

#CREATING THE LIST

#OUTPUT
}


    priors=R$V=error structure
          =R$nu=
          =
          R (random term, fixed term)
          G ()
          B ()
          covariance matrix, variance




          optional list of prior specifications having 3 possible elements:
          R (R-structure) G (G-structure) and B (fixed effects).
          B is a list containing the expected value (mu) and a (co)variance matrix (V) representing the strength of belief:
          the defaults are B$mu=0 and B$V=I*1e+10, where where I is an identity matrix of appropriate dimension.
          The priors for the variance structures (R and G) are lists with the expected (co)variances (V) and degree of belief parameter (nu) for the inverse-Wishart,
          and also the mean vector (alpha.mu) and covariance matrix (alpha.V) for the redundant working parameters. The deafults are nu=0, V=1, alpha.mu=0, and alpha.V=0.
          When alpha.V is non-zero, parameter expanded algorithms are used.