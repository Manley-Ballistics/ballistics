import Foundation

class Trajectory {
    let MAYEWSKI_CONSTANT: Double = 246
    let PI = 3.14159265359
    let impactHeight = 0.0
    
    let moa_value = 1.05
    let retardCoefficientRate = 0.50
    let breakVelocity = 0.0
    
    private(set) var temperature: Double
    private(set) var altitude: Double
    private(set) var pressure: Double
    private(set) var windSpeed: Double
    private(set) var windDirection: Double
    
    private(set) var bulletWeight: Double
    private(set) var ballisticCoefficient: Double
    
    private(set) var distanceOfInterest: Double
    private(set) var muzzleVelocity: Double
    private(set) var scopeHeight: Double
    private(set) var zeroRange: Double
    
    private(set) var adjustedCoefficient: Double
    private(set) var retardCoefficient: Double
    private(set) var adjretardCoefficient: Double
    private(set) var dropAtZero: Double
    
    private(set) var impactSpeed: Double
    private(set) var energy: Double
    private(set) var dropInches: Double
    private(set) var path: Double // Use this to know actual drop in inches at the distnace
    private(set) var elevnMOA: Double
    private(set) var windageMOA: Double
    private(set) var timeOfFlight: Double
    
    init(distanceOfInterest: Double, muzzleVelocity: Double, bulletWeight: Double, ballisticCoefficient: Double, windSpeed: Double, windDirection: Double, temperature: Double, altitude: Double, pressure: Double, scopeHeight: Double, zeroRange: Double) {
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
        
        let p1: Double = ballisticCoefficient*(460+temperature)
        self.adjustedCoefficient = p1/(519-altitude/280)*exp(altitude/31654)*(2-pressure/pressure)
        self.retardCoefficient = ballisticCoefficient*MAYEWSKI_CONSTANT*pow(muzzleVelocity,0.45)
        self.adjretardCoefficient = retardCoefficient*(460+temperature)/(519-altitude/280)*exp(altitude/31654)*(2-pressure/pressure)
        self.dropAtZero = pow(((41.68/muzzleVelocity)/((1/(0+zeroRange))-(1/(self.adjretardCoefficient-(0.75+0.00006*zeroRange)*self.retardCoefficientRate*zeroRange)))),2)
        
        self.impactSpeed = muzzleVelocity*pow((1-3*self.retardCoefficientRate*(distanceOfInterest/self.adjretardCoefficient)),(1/self.retardCoefficientRate))
        self.energy = bulletWeight*pow(impactSpeed,2)/450380
        self.dropInches = pow(((41.68/muzzleVelocity)/((1/(0+distanceOfInterest))-(1/(adjretardCoefficient-(0.75+0.00006*distanceOfInterest)*self.retardCoefficientRate*distanceOfInterest)))),2)
        self.path = -(dropInches+scopeHeight)+(dropAtZero+scopeHeight+impactHeight)*distanceOfInterest/zeroRange
        self.elevnMOA = -path/distanceOfInterest/moa_value*100
        self.windageMOA = (79.2*distanceOfInterest*windSpeed/muzzleVelocity/(self.adjretardCoefficient/distanceOfInterest-1-self.retardCoefficientRate))/distanceOfInterest/moa_value*100*sin(windDirection/12*2*PI)
        self.timeOfFlight = (self.adjretardCoefficient/muzzleVelocity)/(1-self.retardCoefficientRate)*(pow((muzzleVelocity/self.impactSpeed), (1-self.retardCoefficientRate))-1)
    }
    
    func runCalc() {
        let p1: Double = self.ballisticCoefficient*(460+self.temperature)
        self.adjustedCoefficient = p1/(519-self.altitude/280)*exp(self.altitude/31654)*(2-self.pressure/self.pressure)
        self.retardCoefficient = self.ballisticCoefficient*MAYEWSKI_CONSTANT*pow(self.muzzleVelocity,0.45)
        self.adjretardCoefficient = self.retardCoefficient*(460+self.temperature)/(519-self.altitude/280)*exp(self.altitude/31654)*(2-self.pressure/self.pressure)
        self.dropAtZero = pow(((41.68/self.muzzleVelocity)/((1/(0+self.zeroRange))-(1/(self.adjretardCoefficient-(0.75+0.00006*self.zeroRange)*self.retardCoefficientRate*self.zeroRange)))),2)
        
        self.impactSpeed = self.muzzleVelocity*pow((1-3*self.retardCoefficientRate*(self.distanceOfInterest/self.adjretardCoefficient)),(1/self.retardCoefficientRate))
        self.energy = self.bulletWeight*pow(self.impactSpeed,2)/450380
        self.dropInches = pow(((41.68/self.muzzleVelocity)/((1/(0+self.distanceOfInterest))-(1/(self.adjretardCoefficient-(0.75+0.00006*self.distanceOfInterest)*self.retardCoefficientRate*self.distanceOfInterest)))),2)
        self.path = -(self.dropInches+self.scopeHeight)+(self.dropAtZero+self.scopeHeight+self.impactHeight)*self.distanceOfInterest/self.zeroRange
        self.elevnMOA = -self.path/self.distanceOfInterest/self.moa_value*100
        self.windageMOA = (79.2*self.distanceOfInterest*self.windSpeed/self.muzzleVelocity/(self.adjretardCoefficient/self.distanceOfInterest-1-self.retardCoefficientRate))/self.distanceOfInterest/self.moa_value*100*sin(self.windDirection/12*2*PI)
        self.timeOfFlight = (self.adjretardCoefficient/self.muzzleVelocity)/(1-self.retardCoefficientRate)*(pow((self.muzzleVelocity/self.impactSpeed), (1-self.retardCoefficientRate))-1)
    }
    
    /* 
    RECALCULATE BASED ON COMMON VALUE CHANGES
    */
    func recalculate(distanceOfInterest: Double) {
        self.distanceOfInterest = distanceOfInterest
    
        self.runCalc()
    }
    
    func recalculate(distanceOfInterest: Double, windSpeed: Double, windDirection: Double) {
        self.distanceOfInterest = distanceOfInterest
        self.windSpeed = windSpeed
        self.windDirection = windDirection
    
        self.runCalc()
    }
    
    func recalculate(temperature: Double, altitude: Double, pressure: Double) {
        self.temperature = temperature
        self.altitude = altitude
        self.pressure = pressure
        
        self.runCalc()
    }
    
    func recalculate(scopeHeight: Double, zeroRange: Double) {
        self.scopeHeight = scopeHeight
        self.zeroRange = zeroRange
        
        self.runCalc()
    }
    
    func recalculate(bulletWeight: Double, muzzleVelocity: Double, ballisticCoefficient: Double) {
        self.bulletWeight = bulletWeight
        self.muzzleVelocity = muzzleVelocity
        self.ballisticCoefficient = ballisticCoefficient
        
        self.runCalc()
    }
    
    /*
    COMPLETE RECALCULATION
    */
    func recalculate(distanceOfInterest: Double, muzzleVelocity: Double, bulletWeight: Double, ballisticCoefficient: Double, windSpeed: Double, windDirection: Double, temperature: Double, altitude: Double, pressure: Double, scopeHeight: Double, zeroRange: Double) {
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
    }
}

// Example use
let bullet = Trajectory(distanceOfInterest: 500, muzzleVelocity: 2820, bulletWeight: 120, ballisticCoefficient: 0.297, windSpeed: 10.0, windDirection: 9.0, temperature: 65, altitude: 500, pressure: 1000, scopeHeight: 2.0, zeroRange: 200)
