# Uploads the raster file created by BoundaryConditions.R
# and simululates a reaction-diffusion equation in two spatial
# dimensions.

# Script created by Kevin Bell, knb226@mun.ca and modified by
# Amy Hurford ahurford@mun.ca
#******************************************************

# To solve the PDE requires the deSolve package
require(deSolve)

# Remove variables
rm(list = ls())

# Set working directory and load the raster file 
setwd("~/Dropbox/Boundary_Conditions/")
# This map is of the Burin peninsula
load("Burin.RData")

# This function is the reaction-diffusion equation to be solved
mod2D <- function (time, state, pars, Nx,Ny, dx,dy,rastermap) {
  NN <- Nx*Ny
  # State is reshaped into a matrix, since we are working in 2D space
  Popn <- matrix(nrow = Nx, ncol = Ny,state[1:NN])
  # The raster map contains information about where the shore is
  rasterMx <-matrix(nrow = Nx, ncol = Ny,rastermap@data@values)
  # These will be used for the zero-flux boundary conditions
  zerox <- rep(0, Nx)
  zeroy <- rep(0, Ny)
  
  with (as.list(pars), {
    # logistic population growth
    dPopn <- r * Popn * (1- Popn/K) 
    
    #1. Fluxes in x-direction; zero fluxes near boundaries
    FluxPopn <- -Dx * rbind(zeroy,(Popn[2:Nx,] - Popn[1:(Nx-1),]), zeroy)/dx
    
    # Terms are reaction (logistic growth), diffussion (finite differences), and advection
    # in the x-direction.
    dPopn    <- dPopn - (FluxPopn[2:(Nx+1),] - FluxPopn[1:Nx,])/dx - ax*FluxPopn[1:Nx,]

    ## 2. Fluxes in y-direction; zero fluxes near boundaries
    FluxPopn <- -Dy * cbind(zerox,(Popn[,2:Ny] - Popn[,1:(Ny-1)]), zerox)/dy
    
    ## Terms are dPopn (logistic growth and movement in x-direction), diffusion in
    # y-direction (via finite differences), and advection (in the y-direction)
    dPopn    <- dPopn - (FluxPopn[,2:(Ny+1)] - FluxPopn[,1:Ny])/dy - ay*FluxPopn[,1:Ny]
    
    # This is a tolerance threshold: it is very slow to sum across all space,
    # but much quicker if small values are set to zero.
    dPopn[dPopn<10e-8]<-0
    sum1=sum(dPopn)
    
    # The boundary is dealt with very heuristically here: any density that
    # corresponds to the shoreline is set to zero.
    dPopn[rastermap@data@values==1]<-0
    
    # The population density is renormalized to deal with having lost density
    # by setting any density that fell on the boundary to zero.
    dPopn=dPopn*sum1/sum(dPopn)
    return(list(as.vector(dPopn)))
  })
}

# Parameters
# r: net reproductive rate (logistic growth)
# K: carrying capacity
# Dx: diffussion in the x-direction
# Dy: diffussion in the y-direction
# ax: advection in the x-direction
# ay: advection in the y-direction
pars    <- c(r= 0.2,K= 1, Dx=1e-4, Dy = 1e-4,ay = 50, ax = 0)   

Rx  <- (rastermap@extent@xmax-rastermap@extent@xmin) # length in the x direction
Ry <- (rastermap@extent@ymax-rastermap@extent@ymin) # length in the y direction
Ny  <- rastermap@nrows # number of discritizations in the x-direction
Nx <- rastermap@ncols 
dx <- Rx/Nx # discretization width                    
dy <- Ry/Ny 
NN <- Nx*Ny # total number of grid cells               

## initial conditions
yini    <- rep(0,NN)
# Cell index where initial density is 1.
cc      <- 24000
yini[cc] <- 1
# Length of simulation
T = 50
times   <- seq(0, T, by = 1)

# This is solving the PDE
out <- ode.2D(y = yini, times = times, func = mod2D, parms = pars,
              dimens = c(Nx, Ny), names = c("Popn"),
              Nx = Nx,Ny=Ny, dx = dx,dy=dy,rastermap=rastermap, method = rkMethod("rk45ck"))

# ****************************************************************************************************
# Plot in a 2x2 array
# Note that it will require some work to get the axes limits correct.
# Also, it would be nice to make an animation, but that also requires
# a bit of work https://www.r-bloggers.com/animation-in-r/
par(mfrow = c(2,2)) 


for(i in seq(1,4)){
  # Extracting the data to make the 4 different plots at t = T/4, T/2, ... T.
  # (rounded to an interger value of time)
data = matrix(nrow = Nx, ncol = Ny, unname(out[floor(i*T/4),1:NN+1]))

# Create blank raster of the same resolution and distribution as the rastermap
# The output of the PDE solver will be inserted onto this blank map
rastermapBlank.mask = rastermap
rastermapBlank.mask@data@values=0
# Projecting the PDE solution onto the blank map
rastermapBlank.mask@data@values = data
# Remove expretemly low distribution values - this just ignores making
# some calcuations (below) that aren't super interesting (low density).
 rastermapBlank.mask[rastermapBlank.mask@data@values<(1e-8)]=0

# Combine the rastermap and PDE solution into a raster stack object
layers=c(rastermap,rastermapBlank.mask)
rasterComb = stack(layers)
 
# overlay the two layers of the raster stack as one raster layer
# The coastline is coded as 1 in rastermap and everywhere else as 0.
# The population density is small values > 0. The code below says for
# a point on the coastline assign a value of 1.1, otherwise plot the
# popn density for the PDE. It would be possible to finesse this a bit
# more to get better colors/contrast.
rasterDist = overlay(rasterComb,fun=function(x,y){ifelse(x-y>0,1.1,y)})

# Make the plot. Rescaling the size of your plotting window can help with
# the visualization.
plot(rasterDist,col = topo.colors(20), ylim = c(46.5, 48))
}