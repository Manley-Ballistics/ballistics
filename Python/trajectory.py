#
#  TRAJECTORY.PY
#  CREATED BY GREGORY MANLEY ON 2/6/2023.
#  UPDATED BY GREGORY MANLEY ON 2/6/2023.
#
#  THIS SOURCE CODE IS LICENSED UNDER THE LICENSE FOUDN IN LICENSE FILE IN THE
#  ROOT DIRECTORY OF THIS REPO
#  
#  TODO:
#       ADD DEBUG OUTPUT
#
###############################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@GGGGGGGGGGGGGGGGGGGG@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@GGGGGGGGGGGGGGGGGGGG@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@GGG@@@@@@@@@@@@@@GGG@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@GGG@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@GGG@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@GGG@@@@@@@@@@GGGGGGG@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@GGG@@@@@@@@@@GGGGGGG@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@GGG@@@@@@@@@@@@@@GGG@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@GGG@@@@@@@@@@@@@@GGG@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@GGGGGGGGGGGGGGGGGGGG@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@GGGGGGGGGGGGGGGGGGGG@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@MMMM@@@@@@@@@@@@MMMM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@MMMMM@@@@@@@@@@MMMMM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@MMMMMM@@@@@@@@MMMMMM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@MMM@MMM@@@@@@MMM@MMM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@MMM@@MMM@@@@MMM@@MMM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@MMM@@@MMM@@MMM@@@MMM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@MMM@@@@MMMMMM@@@@MMM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@MMM@@@@@@@@@@@@@@MMM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@MMM@@@@@@@@@@@@@@MMM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@MMM@@@@@@@@@@@@@@MMM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@MMM@@@@@@@@@@@@@@MMM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
###############################################################################
from __future__ import annotations
import math

class Trajectory:
    MAYEWSKI_CONSTANT = 246
    PI = 3.14159265359
    impactHeight = 0.0
    
    moa_value = 1.05
    retardCoefficientRate = 0.50
    breakVelocity = 0.0

    def __init__(self, distanceOfInterest, muzzleVelocity, bulletWeight, ballisticCoefficient, windSpeed, windDirection, temperature, altitude, pressure, scopeHeight, zeroRange):
        self.temperature = temperature
        self.altitude = altitude
        self.pressure = pressure
        self.windSpeed = windSpeed
        self.windDirection = windDirection
    
        self.bulletWeight = bulletWeight
        self.ballisticCoefficient = ballisticCoefficient
    
        self.distanceOfInterest = distanceOfInterest
        self.muzzleVelocity = muzzleVelocity
        self.scopeHeight = scopeHeight
        self.zeroRange = zeroRange
        
        p1 = ballisticCoefficient*(460+temperature)
        self.adjustedCoefficient = p1/(519-altitude/280)*math.exp(altitude/31654)*(2-pressure/pressure)
        self.retardCoefficient = ballisticCoefficient*self.MAYEWSKI_CONSTANT*pow(muzzleVelocity,0.45)
        self.adjretardCoefficient = self.retardCoefficient*(460+temperature)/(519-altitude/280)*math.exp(altitude/31654)*(2-pressure/pressure)
        self.dropAtZero = pow(((41.68/muzzleVelocity)/((1/(0+zeroRange))-(1/(self.adjretardCoefficient-(0.75+0.00006*zeroRange)*self.retardCoefficientRate*zeroRange)))),2)
        
        self.impactSpeed = muzzleVelocity*pow((1-3*self.retardCoefficientRate*(distanceOfInterest/self.adjretardCoefficient)),(1/self.retardCoefficientRate))
        self.energy = bulletWeight*pow(self.impactSpeed,2)/450380
        self.dropInches = pow(((41.68/muzzleVelocity)/((1/(0+distanceOfInterest))-(1/(self.adjretardCoefficient-(0.75+0.00006*distanceOfInterest)*self.retardCoefficientRate*distanceOfInterest)))),2)
        self.path = -(self.dropInches+scopeHeight)+(self.dropAtZero+scopeHeight+self.impactHeight)*distanceOfInterest/zeroRange
        self.elevnMOA = -self.path/distanceOfInterest/self.moa_value*100
        self.windageMOA = (79.2*distanceOfInterest*windSpeed/muzzleVelocity/(self.adjretardCoefficient/distanceOfInterest-1-self.retardCoefficientRate))/distanceOfInterest/self.moa_value*100*math.sin(windDirection/12*2*self.PI)
        self.timeOfFlight = (self.adjretardCoefficient/muzzleVelocity)/(1-self.retardCoefficientRate)*(pow((muzzleVelocity/self.impactSpeed), (1-self.retardCoefficientRate))-1)

    def runCalc(self):
        p1 = self.ballisticCoefficient*(460+self.temperature)
        self.adjustedCoefficient = p1/(519-self.altitude/280)*math.exp(self.altitude/31654)*(2-self.pressure/self.pressure)
        self.retardCoefficient = self.ballisticCoefficient*self.MAYEWSKI_CONSTANT*pow(self.muzzleVelocity,0.45)
        self.adjretardCoefficient = self.retardCoefficient*(460+self.temperature)/(519-self.altitude/280)*math.exp(self.altitude/31654)*(2-self.pressure/self.pressure)
        self.dropAtZero = pow(((41.68/self.muzzleVelocity)/((1/(0+self.zeroRange))-(1/(self.adjretardCoefficient-(0.75+0.00006*self.zeroRange)*self.retardCoefficientRate*self.zeroRange)))),2)
        
        self.impactSpeed = self.muzzleVelocity*pow((1-3*self.retardCoefficientRate*(self.distanceOfInterest/self.adjretardCoefficient)),(1/self.retardCoefficientRate))
        self.energy = self.bulletWeight*pow(self.impactSpeed,2)/450380
        self.dropInches = pow(((41.68/self.muzzleVelocity)/((1/(0+self.distanceOfInterest))-(1/(self.adjretardCoefficient-(0.75+0.00006*self.distanceOfInterest)*self.retardCoefficientRate*self.distanceOfInterest)))),2)
        self.path = -(self.dropInches+self.scopeHeight)+(self.dropAtZero+self.scopeHeight+self.impactHeight)*self.distanceOfInterest/self.zeroRange
        self.elevnMOA = -self.path/self.distanceOfInterest/self.moa_value*100
        self.windageMOA = (79.2*self.distanceOfInterest*self.windSpeed/self.muzzleVelocity/(self.adjretardCoefficient/self.distanceOfInterest-1-self.retardCoefficientRate))/self.distanceOfInterest/self.moa_value*100*math.sin(self.windDirection/12*2*self.PI)
        self.timeOfFlight = (self.adjretardCoefficient/self.muzzleVelocity)/(1-self.retardCoefficientRate)*(pow((self.muzzleVelocity/self.impactSpeed), (1-self.retardCoefficientRate))-1)

    #
    # RECALCULATE BSAED ON COMMON VALUE CHANGES
    #
    def recalculate(self, distanceOfInterest):
        self.distanceOfInterest = distanceOfInterest

        self.runCalc()

    def recalculate(self, temperature, altitude, pressure):
        self.temperature = temperature
        self.altitude = altitude
        self.pressure = pressure

        self.runCalc()

    def recalculate(self, scopeHeight, zeroRange):
        self.scopeHeight = scopeHeight
        self.zeroRange = zeroRange

        self.runCalc()

    def recalculate(self, bulletWeight, muzzleVelocity, ballisticsCoefficient):
        self.bulletWeight = bulletWeight
        self.muzzleVelocity = muzzleVelocity
        self.ballisticCoefficient = ballisticsCoefficient

        self.runCalc()

    #
    # COMPLETE RECALCULATION
    #
    def recalculate(self, distanceOfInterest, muzzleVelocity, bulletWeight, ballisticCoefficient, windSpeed, windDirection, temperature, altitude, pressure, scopeHeight, zeroRange):
        self.temperature = temperature
        self.altitude = altitude
        self.pressure = pressure
        self.windSpeed = windSpeed
        self.windDirection = windDirection
    
        self.bulletWeight = bulletWeight
        self.ballisticCoefficient = ballisticCoefficient
    
        self.distanceOfInterest = distanceOfInterest
        self.muzzleVelocity = muzzleVelocity
        self.scopeHeight = scopeHeight
        self.zeroRange = zeroRange
    
        self.runCalc()

# EXAMPLE USE
bullet = Trajectory(500, 2820, 120, 0.297, 10.0, 9.0, 65, 500, 1000, 2.0, 200)
print(bullet.dropInches)
