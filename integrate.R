library(cubature)
library(fields)
library(plotfunctions)

# generate training grid
nx <- 25
x <- seq(0,1,length.out=nx)
xmat <- expand.grid(x,x,x)

# generate prediction grid
np <- nx - 1
x_pred <- seq(0.01,0.99,length.out=np)
xmat_pred <- expand.grid(x_pred, x_pred, x_pred)


# deterministic output
f <- function(x){(x[1]+1)*cos(pi*x[2])+3}

y <- apply(xmat, 1, f)
ya <- array(y, c(nx, nx))
image.plot(ya)

mean(y)

F_x <- hcubature(f, lowerLimit = c(0,0,0), upperLimit = c(1,1,1), tol = 4e-5)
F_x

# make a function, and sample X and y
f <- function(x){(x[1]+1)*x[2]^0.5+3}
nr <- 5000 # how many random draws to make
x <- matrix(runif(nr*2), nr, 2)
y <- apply(x, 1, f)

# make a color palette
my_palette <- topo.colors(50) # in package grDevices
n_colors <- length(my_palette)

# adjust so legend fits on right
par(mar=c(5,4,4,5)) 

plot(x[,1], x[,2], 
     col=my_palette[ceiling((y-min(y))/
                              (max(y)-min(y))*(n_colors-1))+1],
     pch=16, xlab="x1",ylab="x2"
)

# generate legend from plotfunctions package
gradientLegend(range(y), my_palette, side=4, 
               pos = .5, length = 1, n.seg = 5)

F_x <- hcubature(f, lowerLimit = c(0,0,0), upperLimit = c(1,1,1), tol = 4e-5)
F_x

