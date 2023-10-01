# Derivative of true regression function
f0_prime = function(x){
    (1-2*x) * sin(2*pi/(x+0.5)) / (2*sqrt(x*(1-x))) - 2*pi*sqrt(x*(1-x)) * cos(2*pi/(x+0.5)) / (x+0.5)^2
}
