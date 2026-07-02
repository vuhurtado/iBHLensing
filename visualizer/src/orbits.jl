# {
#  "cells": [],
#  "metadata": {},
#  "nbformat": 4,
#  "nbformat_minor": 5
# }

module Orbit

export times

"""
make array of times based off period
input:
    Period :: Float64 
"""

function times(n,period)
    ts = collect(0:1e-1:n*period) #grid of times from 0 to n*period
    return ts
end


export mean_anomaly
"""
    mean_anomaly(times, t_p, ,period):
    computes the mean anomaly [radians], sometimes called f 
    inputs:
        times :: #array of times Float64 #[]
        time at periapse t_p :: #Float64 #[]
        period P :: #Float64 #[].
"""

function mean_anomaly(times,tp,period)
    return ((2*π)/period).*(times.-tp)
end

export g_ks
"""
g_ks(eccentric_anomaly,e,mean_anomaly):
function that solves kepler's equations (g)
inputs:
    eccentric_anomaly :: #eccentric anomaly
    e :: #eccentricity
    meanan :: #mean anomaly
"""

function g_ks(eccentric_anomaly,e,mean_anomaly)
    return eccentric_anomaly - e*sin(eccentric_anomaly) - mean_anomaly    
end

export g_prime
"""
g_prime(eccentric_anomaly,e)
derivative of g (kepler solver) with respect to the eccentric anomaly
inputs:
    eccentric_anomaly :: #eccentric anomaly
    e :: #eccentricity
"""

function g_prime(eccentric_anomaly,e)
    return 1 - e*cos(eccentric_anomaly)
end

export newton_solver
"""
newton_solver(e,mean_anomaly):
returns a first guess to start a eccentric anomaly calculation
inputs:
    e :: #eccentricity
    mean_anomaly :: #mean anomaly [rad]
"""

function newton_solver(e,mean_anomaly)
    ExA0 = mean_anomaly + 0.85*e*sign(sin(mean_anomaly))
    while abs(g_ks(ExA0,e,mean_anomaly)/g_prime(ExA0,e)) > 1e-10
        eccentric_anomaly = ExA0 - g_ks(ExA0,e,mean_anomaly)/g_prime(ExA0,e)
        ExA0 = eccentric_anomaly
    end
    return ExA0
end

export true_anomaly
"""
true_anomaly(e,eccentric_anomaly):
calculates the true anomaly in []
inputs:
    eccentricity :: #Float64
    ExA :: #eccentric_anomaly 
"""

function true_anomaly(e,eccentric_anomaly)
    f = 2*(atan((((1+e)/(1-e))^0.5)*(tan(eccentric_anomaly/2))))
    # true_array = true_anomaly.(e,eccentric_anomaly)
    return f
end

export eccentric_anomaly
"""
eccentric_anomaly = newton_solver.(e,mean_anomaly)
"""

function eccentric_anomaly(e,mean_anomaly)
    return newton_solver.(e,mean_anomaly)
end

export rad_to_deg
"""
rad_to_deg(rad_angle)
converts angle in radians to degrees
input:
    rad_angle :: Float64
"""

function rad_to_deg(rad_angle)
    return rad_angle*(180/π)
end

export AU_to_metric
"""
converts astronomical units to cm or M
input:
    distance :: Float64
    unit :: string #"cm" or "m"
"""

function AU_to_metric(distance,unit)
    if unit == "cm"
        distance = distance.*(1.496e13) #1 AU = 1.496e13 cm
    end
    if unit == "m"
        distance = distance.*(1.496e11) #i think?
    end
    return distance
end

end