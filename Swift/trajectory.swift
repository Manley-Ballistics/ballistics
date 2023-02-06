import Foundation

class Trajectory {
    let MAYEWSKI_CONSTANT: Double = 246
    let PI = 3.14159265359
    let impactHeight = 0.0
    
    let moa_value = 1.05
    let retardCoefficientRate = 0.50
    let breakVelocity = 0.0
    
    var adjustedCoefficient: Double
    var retardCoefficient: Double
    var adjretardCoefficient: Double
    var dropAtZero: Double
    
    var impactSpeed: Double
    var energy: Double
    var dropInches: Double
    var path: Double
    var elevnMOA: Double
    var windageMOA: Double
    var timeOfFlight: Double
    
    init(distanceOfInterest: Double, muzzleVelocity: Double, bulletWeight: Double, ballisticCoefficient: Double, windSpeed: Double, windDirection: Double, temperature: Double, altitude: Double, pressure: Double, scopeHeight: Double, zeroRange: Double) {
        let p1: Double = ballisticCoefficient*(460+temperature)
        self.adjustedCoefficient = p1/(519-altitude/280)*exp(altitude/31654)*(2-pressure/pressure)
        self.retardCoefficient = ballisticCoefficient*MAYEWSKI_CONSTANT*pow(muzzleVelocity,0.45)
        self.adjretardCoefficient = retardCoefficient*(460+temperature)/(519-altitude/280)*exp(altitude/31654)*(2-pressure/pressure)
        self.dropAtZero = pow(((41.68/muzzleVelocity)/((1/(0+zeroRange))-(1/(self.adjretardCoefficient-(0.75+0.00006*zeroRange)*retardCoefficientRate*zeroRange)))),2)
        
        self.impactSpeed = muzzleVelocity*pow((1-3*retardCoefficientRate*(distanceOfInterest/adjretardCoefficient)),(1/retardCoefficientRate))
        self.energy = bulletWeight*pow(impactSpeed,2)/450380
        self.dropInches = pow(((41.68/muzzleVelocity)/((1/(0+distanceOfInterest))-(1/(adjretardCoefficient-(0.75+0.00006*distanceOfInterest)*retardCoefficientRate*distanceOfInterest)))),2)
        self.path = -(dropInches+scopeHeight)+(dropAtZero+scopeHeight+impactHeight)*distanceOfInterest/zeroRange
        self.elevnMOA = -path/distanceOfInterest/moa_value*100
        self.windageMOA = (79.2*distanceOfInterest*windSpeed/muzzleVelocity/(adjretardCoefficient/distanceOfInterest-1-retardCoefficientRate))/distanceOfInterest/moa_value*100*sin(windDirection/12*2*PI)
        self.timeOfFlight = (adjretardCoefficient/muzzleVelocity)/(1-retardCoefficientRate)*(pow((muzzleVelocity/impactSpeed), (1-retardCoefficientRate))-1)
    }
}

// Example use
let bullet = Trajectory(distanceOfInterest: 500, muzzleVelocity: 2820, bulletWeight: 120, ballisticCoefficient: 0.297, windSpeed: 10.0, windDirection: 9.0, temperature: 65, altitude: 500, pressure: 1000, scopeHeight: 2.0, zeroRange: 200)
