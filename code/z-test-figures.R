#Throwaway file - this generates the graphics for the FDA talk.

# Graphic 1: Power curve for a simple Normal Test with sigma = 1, testing mu > 0, at level .025

png(file="../figures/z-test-1.png")
npoints = 40
mu = seq(-3.25,3.75, by = .5)
# Power of a test with with sigma = 1, alpha .025:
# = P( N(mu, 1) > z_.975 )
# = P( N(0,1) > -mu + z_.975    )
pow = 1 - pnorm(- mu + qnorm(.975))
plot(type = "n",x = mu, y = pow, cex = .2, xlab = "True Value of mu",ylab = "Power", main = "Power of Standard Normal Test")

temp_mu = seq(-3.25,3.75, length.out = npoints*100)
pow_temp_mu = 1 - pnorm(- temp_mu + qnorm(.975))
lines(x = temp_mu, y = pow_temp_mu, lty = 3)
abline(a = .025,b =0)
dev.off()

#Perhaps add a rectangular box here, on the powerpoint, to show that we are zooming in on it!

# Graphic 2a: Type I Error curve for the same normal test (zeroing out at mu = 0). Sparse point set. As if you had infinitely many Monte Carlo samples!
png(file="../figures/z-test-2.png")
npoints = 3
mu = seq(-1,0, length.out = npoints)
stepsize = mu[2] - mu[1]
mu= mu - stepsize/2
pow = 1 - pnorm(- mu + qnorm(.975))
plot(type = "n",x = mu, y = pow, cex = .5, ylim = c(0, .03), xlim = c(-1,0.25),
       xlab = "True Value of mu",ylab = "Type I Error", main = "Type I Error of Standard Normal Test")
abline(a = .025, b = 0)
abline(a = 0, b = 999999)

temp_mu = seq(-1,0, length.out = npoints*100)
pow_temp_mu = 1 - pnorm(- temp_mu + qnorm(.975))
lines(x = temp_mu, y = pow_temp_mu, lty = 3)
abline(a = .025,b =0)
dev.off()

# Graphic 2b: Sparse point set. As if you had infinitely many Monte Carlo samples!
png(file="../figures/z-test-3.png")
npoints = 3
mu = seq(-1,0, length.out = npoints)
stepsize = mu[2] - mu[1]
mu= mu - stepsize/2
pow = 1 - pnorm(- mu + qnorm(.975))
plot(x = mu, y = pow, cex = .5, ylim = c(0, .03), xlim = c(-1,0.25),
     xlab = "True Value of mu",ylab = "Type I Error", main = "Type I Error of Standard Normal Test")
abline(a = .025, b = 0)
abline(a = 0, b = 999999)
dev.off()

# Graphic 3: Using the first derivative, add knowledge of the slope. Represent this as lines, going in each direction, through the relevant points
png(file="../figures/z-test-5.png")
plot(x = mu, y = pow, cex = .5, ylim = c(0, .03), xlim = c(-1,0.25),
     xlab = "True Value of mu",ylab = "Type I Error", main = "Type I Error of Standard Normal Test")
abline(a = .025, b = 0)
abline(a = 0, b = 999999)
epsilon = .001
mu_plus_eps = mu + epsilon
pow_plus_eps = 1 - pnorm(- mu_plus_eps + qnorm(.975))
pow_derivative = (pow_plus_eps - pow)/epsilon

upper = pow + pow_derivative*stepsize/2
lower = pow - pow_derivative*stepsize/2

for(i in 1:npoints){
  lines(x = mu[i] + stepsize*c(-.5,.5), y = c(lower[i],upper[i]))
}
dev.off()

#. But this is an approximation, not a rigorous upper bound!
# Graphic 3b: Add a plot of the true type I Error over top!
png(file='../figures/z-test-6.png')
plot(x = mu, y = pow, cex = .5, ylim = c(0, .03), xlim = c(-1,0.25),
     xlab = "True Value of mu",ylab = "Type I Error", main = "Type I Error of Standard Normal Test")
abline(a = .025, b = 0)
abline(a = 0, b = 999999)
for(i in 1:npoints){
  lines(x = mu[i] + stepsize*c(-.5,.5), y = c(lower[i],upper[i]))
}
temp_mu = seq(-1,0, length.out = npoints*100)
pow_temp_mu = 1 - pnorm(- temp_mu + qnorm(.975))
lines(x = temp_mu, y = pow_temp_mu, lty = 3)
dev.off()


# Graphic 3c: Show that Taylor's theorem describes the error of this approximation
png(file='../figures/z-test-7.png')
plot(x = mu, y = pow, cex = .5, ylim = c(0, .03), xlim = c(-1,0.25),
     xlab = "True Value of mu",ylab = "Type I Error", main = "Type I Error of Standard Normal Test")
abline(a = .025, b = 0)
abline(a = 0, b = 999999)
for(i in 1:npoints){
  lines(x = mu[i] + stepsize*c(-.5,.5), y = c(lower[i],upper[i]))
}
lines(x = temp_mu, y = pow_temp_mu, lty = 3)

index_point<- length(temp_mu)*15/16
x_point = temp_mu[index_point]
y_point = pow_temp_mu[index_point]
points(x_point, y_point)

x_diff = x_point - mu[length(mu)]
y_diff = y_point - (pow[length(mu)] + pow_derivative[length(mu)]*(x_diff))

x_temp = rep(NA, 50)
y_temp = rep(NA, 50)
s = seq(0,x_diff, length.out = 50)
for(i in 1:50){
  x_temp[i] = x_point - x_diff + s[i]
  y_temp[i] = pow[length(mu)] + s[i]*pow_derivative[length(mu)]+ (y_diff)*(s[i]/x_diff)^2
}
lines(x_temp, y_temp)
dev.off()

# Graphic 4: Add a quadratic curve to the linear point! This is now a rigorous upper bound function.
#Note: I believe the upper bound on the second derivative is equal to 1. The taylor penalty is equal to d^2/2
png(file='../figures/z-test-8.png')
plot(x = mu, y = pow, cex = .5, ylim = c(0, .03), xlim = c(-1,0.25),
     xlab = "True Value of mu",ylab = "Type I Error", main = "Type I Error of Standard Normal Test")
abline(a = .025, b = 0)
abline(a = 0, b = 999999)
for(i in 1:npoints){
  lines(x = mu[i] + stepsize*c(-.5,.5), y = c(lower[i],upper[i]))
}
lines(x = temp_mu, y = pow_temp_mu, lty = 3)
points(x_point, y_point)
lines(x_temp, y_temp)

penalty_length = 20
penalized_sequence <- seq(-.5,.5, length.out = penalty_length)
for(i in 1:npoints){
      penalized_bound = pow[i] + pow_derivative[i]*stepsize*penalized_sequence + (1/2)*(stepsize^2)*(penalized_sequence)^2
      lines(x = mu[i] + stepsize*penalized_sequence, y = penalized_bound)
}
dev.off()

# Graphics 5-8: animate how much better this gets as we change the number of points, showing that it rapidly converges!
smooth_out = function(npoints, fig.idx){
    png(file=paste('../figures/z-test-', fig.idx, '.png', sep=''))

    # Copy/pasted the above!
    mu = seq(-1,0, length.out = npoints)
    stepsize = mu[2] - mu[1]
    mu= mu - stepsize/2
    pow = 1 - pnorm(- mu + qnorm(.975))
    plot(x = mu, y = pow, cex = .5, ylim = c(0, .03), xlim = c(-1,0.25),
         xlab = "True Value of mu",ylab = "Type I Error", main = "Type I Error of Standard Normal Test")
    abline(a = .025, b = 0)
    abline(a = 0, b = 999999)

    epsilon = .001
    mu_plus_eps = mu + epsilon
    pow_plus_eps = 1 - pnorm(- mu_plus_eps + qnorm(.975))
    pow_derivative = (pow_plus_eps - pow)/epsilon

    upper = pow + pow_derivative*stepsize/2
    lower = pow - pow_derivative*stepsize/2

    for(i in 1:npoints){
        lines(x = mu[i] + stepsize*c(-.5,.5), y = c(lower[i],upper[i]))
    }

    temp_mu = seq(-1,0, length.out = npoints*100)
    pow_temp_mu = 1 - pnorm(- temp_mu + qnorm(.975))
    lines(x = temp_mu, y = pow_temp_mu, lty = 3)

    penalty_length = 20
    penalized_sequence <- seq(-.5,.5, length.out = penalty_length)

    for(i in 1:npoints){
        penalized_bound = pow[i] + pow_derivative[i]*stepsize*penalized_sequence + (1/2)*(stepsize^2)*(penalized_sequence)^2
        lines(x = mu[i] + stepsize*penalized_sequence, y = penalized_bound)
    }

    dev.off()
}

smooth_out(5, 9)
smooth_out(10, 10)
smooth_out(20, 11)
smooth_out(40, 12)

# Ok, that is how you would do it in 1 dimension!
# Except, wait, we don't actually know this 0th order and 1st order information. What should we do instead?
# We can make confidence intervals for both of them!
# The 0th order, the point estimates here, are what you're used to seeing estimated by Monte Carlo 

# [Insert new graphic - the points, with vertical confidence band]

# The first order, the point estimates here, can actually be estimated with Monte Carlo as well - 
# The estimator we need to use here, is a certain type of score estimator

# It turns out that we can bound these
# [Insert new graphic - the points are corners of triangles, which have shaded-in arcs of the area 
# that can be reached by the score confidence interval]


#Graphic 8 : error bars on Monte Carlo
require(plotrix)

png(file='../figures/z-test-13.png')
#### Copied from previous!
npoints = 3
mu = seq(-1,0, length.out = npoints)
stepsize = mu[2] - mu[1]
mu= mu - stepsize/2
pow = 1 - pnorm(- mu + qnorm(.975))
plot(x = mu, y = pow, cex = 1, ylim = c(0, .03), xlim = c(-1,0.25),
     xlab = "True Value of mu",ylab = "Type I Error", main = "Type I Error of Standard Normal Test")
abline(a = .025, b = 0)
abline(a = 0, b = 999999)
######

uplimit = pow + (1/100)*sqrt(pow * (1 - pow))
lowlimit = pow - (1/100)*sqrt(pow * (1 - pow))
plotCI(mu, pow, ui = uplimit, li=lowlimit, add = TRUE)
dev.off()

#Graphic 9: error bars on the score estimate
png(file='../figures/z-test-14.png')
plot(x = mu, y = pow, cex = 1, ylim = c(0, .03), xlim = c(-1,0.25),
     xlab = "True Value of mu",ylab = "Type I Error", main = "Type I Error of Standard Normal Test")
abline(a = .025, b = 0)
abline(a = 0, b = 999999)

epsilon = .001
mu_plus_eps = mu + epsilon
pow_plus_eps = 1 - pnorm(- mu_plus_eps + qnorm(.975))
pow_derivative = (pow_plus_eps - pow)/epsilon

width = 1/100
vertical = (1/100)*sqrt(pow * (1 - pow))
for(i in 1:npoints){
  x = c(mu[i] , mu[i] , mu[i] + stepsize/2 , mu[i]+stepsize/2)
  y = c(pow[i] - vertical[i], pow[i] + vertical[i],
        pow[i] + (stepsize/2)*pow_derivative[i] + width*stepsize/2 + vertical[i], 
        pow[i] + (stepsize/2)*pow_derivative[i] - width*stepsize/2 - vertical[i])
  polygon(x,y,border=NA,col=blues9[3])
}
for(i in 1:npoints){
  x = c(mu[i], mu[i], mu[i] - stepsize/2, mu[i] - stepsize/2)
  y = c(pow[i] - vertical[i], pow[i] + vertical[i], 
        pow[i] - (stepsize/2)*pow_derivative[i] + width*stepsize/2 + vertical[i], 
        pow[i] - (stepsize/2)*pow_derivative[i] - width*stepsize/2 - vertical[i])
  polygon(x,y,border=NA,col=blues9[3])
}
plotCI(mu, pow, ui = uplimit, li=lowlimit, add = TRUE)
dev.off()

# In order to come up with a final overall bound:
# we can take a 1 - delta/2 confidence bound in Monte Carlo for the 0th order;
# a 1 - delta/2 confidence bound in Monte Carlo for the score;
# and, using this process, compute the worst upper bound over both of those confidence intervals.
