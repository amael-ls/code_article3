
#### Aim of prog: Generate n random numbers such as their sum is alpha
## Comment
# This program generates n random densities of cohort such that the basal area is
# alpha m^2/ha = alpha x 10e-4 (dimensionless number)
# For instance, to generate a basal area 25 m^2/ha = 2.5e-3 (dimensionless number)
# with 150 individuals in a plot of 400 m^2 we would do:
#	- Compute the basal area (BA) for 400 m^2: BA = 400 * 2.5e-3 = 1
#	- Generate 150 cohorts, each of them is a density that will correspond to a dbh in Ω
#	- Their sum will be one, which corresponds to a BA of 25 in a 400 m^2 plot
#
# The dirichlet generator function is from MCMCpack
# See https://reference.wolfram.com/language/ref/DirichletDistribution.html for more explanation on this distribtion

#### Probability density function of Dirichlet
ddirichlet = function(x, alpha)
{
	dirichlet1 = function(x, alpha)
	{
		logD = sum(lgamma(alpha)) - lgamma(sum(alpha))
		s = sum((alpha-1)*log(x))
		exp(sum(s)-logD)
	}

	# make sure x is a matrix
	if(!is.matrix(x))
		if(is.data.frame(x))
			x = as.matrix(x)
		else
			x = t(x)
	
	if(!is.matrix(alpha))
		alpha = matrix( alpha, ncol=length(alpha), nrow=nrow(x), byrow=TRUE)

	if( any(dim(x) != dim(alpha)) )
		stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")

	pd = vector(length=nrow(x))
	for(i in 1:nrow(x))
		pd[i] = dirichlet1(x[i,],alpha[i,])

	# Enforce 0 <= x[i,j] <= 1, sum(x[i,]) = 1
	pd[ apply( x, 1, function(z) any( z <0 | z > 1)) ] = 0
	pd[ apply( x, 1, function(z) all.equal(sum( z ),1) !=TRUE) ] = 0
	return(pd)
}

#### Random number generator for Dirichlet
rdirichlet = function(n, alpha)
{
	l = length(alpha)
	x = matrix(rgamma(l*n,alpha), ncol = l, byrow = TRUE)
	sm = x %*% rep(1,l)
	return(x/as.vector(sm))
}

#### Generate numbers
## Parameters
nbDensities = 150
nbDraws = 10

BA = 25
plotArea = 400
conversion = 1e-4

sumTo = rep(BA * plotArea * conversion, nbDensities)

## Generate tree densities
treeDensities = rdirichlet(nbDraws, sumTo)

## Check cohorts basal area
rowSums(treeDensities)
