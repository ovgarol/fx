# SPDX-FileCopyrightText: 2026 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de>

library(readxl)
library(quantreg)
library(scales)
library(lubridate)
library(deSolve)
library(GA)
library(latex2exp)

####


setwd("~/Documents/articles/stress_generic")
par(mfrow=c(1,1))

# Parameters
r = 0.5    # recovery rate
c = 0.5    # stress increment f(y) 
mu = 0.4   # growth rate for x
k = 0.75  # parameter of mortality

a = 1 
b = 0

# Define the system of equations
dynamics = function(f, x) {
  df.dt = r*f*(1-f) - c*f
  dx.dt = mu*x*(1-x) - k*(1-f)^a * x^b *x
  return(c(df.dt, dx.dt))
}

#dynamics = function(x, y) {c(y, mu * ( 1 - x^2 ) * y - x )}

equilibrium_points = rbind(c(0,0), c(1-c/r,0), c(0,1-k/mu), c(1-c/r,1-k/mu*c/r))

# Create a grid of (sigma, x) values
sigma_seq = seq(0, 1, length.out = 20)
x_seq = seq(0, 1, length.out = 20)
grid = expand.grid(sigma = sigma_seq, x = x_seq)
arrows_x = numeric(nrow(grid))
arrows_y = numeric(nrow(grid))
magnitudes = numeric(nrow(grid))

for (i in 1:nrow(grid)) {
  derivatives = dynamics(grid$sigma[i], grid$x[i])
  arrows_x[i] = derivatives[1]
  arrows_y[i] = derivatives[2]
  magnitudes[i] = sqrt(arrows_x[i]^2 + arrows_y[i]^2)  # Calculate magnitude
}


vector_field = function(
    f,  # Function describing the vector field
    xmin=0, xmax=1, ymin=0, ymax=1,
    width=600, height=600,
    iterations=50,
    epsilon=.01,
    trace=F
) {
  z = matrix(runif(width*height)>0.5,nr=height)
  #z = outer(1:width, 1:height, function(i, j) (i + j) %% 2)

  i_to_x = function(i) xmin + i / width  * (xmax - xmin)
  j_to_y = function(j) ymin + j / height * (ymax - ymin)
  x_to_i = function(x) pmin( width,  pmax( 1, floor( (x-xmin)/(xmax-xmin) * width  ) ) )
  y_to_j = function(y) pmin( height, pmax( 1, floor( (y-ymin)/(ymax-ymin) * height ) ) )
  i = col(z)
  j = row(z)
  x = i_to_x(i)
  y = j_to_y(j)
  res = z
  
  if(iterations>0)for(k in 1:iterations) {
    v = matrix( f(x, y), nc=2 )
    x = x+epsilon*v[,1]
    y = y+epsilon*v[,2]
    i = x_to_i(x)
    j = y_to_j(y)
    res = res + z[cbind(i,j)]
    if(trace) {
      cat(k, "/", iterations, "\n", sep="")
      dev.hold()
      image(t(res))
      dev.flush()
    }
  }
  if(trace) {
    dev.hold()
    image(res>quantile(res,.6), col=0:1)
    dev.flush()
  }
  return((res-min(res))/(max(res)-min(res)))
}

# Sample data

res = vector_field(
  dynamics,
  xmin=0, xmax=1, ymin=0, ymax=1,
  #xmin=-3, xmax=3, ymin=-3, ymax=3,
  width=100, height=100,
  iterations=100,
  epsilon=0.1,
  trace=F
)

par(mai=c(1,1,1,1),mfrow=c(2,2),las=1)
# Plot the phase plane
plot(grid$sigma, grid$x, type = "n", xlab = TeX('fitness, $f$'), 
     ylab = TeX('normalized area coverage, $x$'), 
     #xaxs='i', yaxs='i',
     xlim = c(0, 1), ylim = c(0, 1), main = "Stability plot")
image(t(log1p(res)),col=colorRampPalette(c('orange','gold','cyan3','dodgerblue3'))(100),add=T,useRaster=T)
#image(t(log1p(res)),col=hcl.colors(20,'mako',rev=F),add=T,useRaster=T)


arrow.factor = 0.3
# Add arrows to indicate direction of flow with color based on magnitude
if(F)arrows(grid$sigma, grid$x, grid$sigma + arrows_x * arrow.factor, grid$x + arrows_y * arrow.factor, 
            length = 0.05, col = 'white',lwd=1)

# Add labels for equilibrium points
points(equilibrium_points, col = "black", pch = c('x','O','O','★'))
text(equilibrium_points[, 1], equilibrium_points[, 2], labels = c("E1", "E2", "E3","E4"), pos = c(4,3,4,1),col='black',font=2)



stabi = function(r.c, mu.k) {
  ifelse(r.c < 1 & mu.k < 1, 1,
         ifelse(r.c > 1 & mu.k < (1 / r.c), 2,
                ifelse(r.c < 1 & mu.k > 1, 3,
                       ifelse(r.c > 1 & (1 + mu.k) > (1 / r.c), 4, NA)
                )
         )
  )
}


x = seq(-2,2,by=0.015)
y = seq(-2,2,by=0.015)
z = outer(x,y)*0


for(i in 1:length(x))for(j in 1:length(x))z[i,j] = stabi(exp(x[i]),exp(x[j])) 

image(exp(x),exp(y),z,log='xy',col=c('orange','gold','cyan3','dodgerblue3'),
      xlab='r/c',
      ylab=TeX('$\\mu/k$'),
      main='Stable equilibrium point')
  
text(0.5,0.5,'E1',font=2)
text(1.25,0.5,'E2',font=2)
text(0.5,1.5,'E3',font=2)
text(1.5,1.5,'E4',col='white',font=2)

points(max(min(r/c,5),0.1),max(min(mu/k,5),0.1),pch=19)
lines(r/(c*c(0.05,1)),c(0,0)+mu/k,pch=19)


######
t = seq(0,100,by=0.1)
state=c(f=0.95,x=0.95)
stress.approx = function(t) 1*abs(sin(t/3)) + 0.5*rnorm(1)

system.equations = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    Y = stress.approx(t) 
    df.dt = r*f*(1-f) - c*Y*f
    dx.dt = mu*x*(1-x) - k*(1-f)^a * x^b *x
    list(c(df.dt, dx.dt))
  })
}


#par(mfrow=c(2,1))
plot(t,0*t,pch=21,type='n',main='relative area coverage',ylab='',xlab='',ylim=c(0.001,1),cex=2,log='y')
#points(data0$time,data0$s.real,pch=21,bg='red',col='white',cex=2)
for(i in 1:5){
  cc = ode(y = state, times = t, func = system.equations, parms = c(r=r,mu=mu,c=c,k=k)*(1+0.1*rnorm(4,0.1)),method = 'euler')
  cc = as.data.frame(cc)
  lines(cc$time,cc$x,lwd=0.5*i,lty=1)
  lines(cc$time,cc$f,col='red',lwd=0.5*i,lty=1)
}

plot(0,xlim=c(0,1),ylim=c(0,1))

for(i in 1:5){
  cc = ode(y = state, times = t, func = system.equations, parms = c(r=r,mu=mu,c=c,k=k)*(1+0.1*rnorm(4,0.1)),method = 'euler')
  cc = as.data.frame(cc)
  # lines(cc$time,cc$x,lwd=0.5*i,lty=1)
  # lines(cc$time,cc$f,col='red',lwd=0.5*i,lty=1)
  lines(cc$f,cc$x,col=alpha(i,0.75))
}
#lines(t,stress.approx(t),col=alpha('blue',0.2))
points(equilibrium_points, col = "black", pch = c('x','O','O','★'))
text(equilibrium_points[, 1], equilibrium_points[, 2], labels = c("E1", "E2", "E3","E4"), pos = c(4,3,4,1),col='black',font=2)


