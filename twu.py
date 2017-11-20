from scipy.optimize import newton
import math


class Twu:
    def __init__(self, boil_temp, specific_gravity):
        self.boilTemp = boil_temp
        self.specificGravity = specific_gravity
        self.alkaneCriticalTemperature = self.calc_alkane_critical_temp()
        self.alpha = self.calc_alpha()
        self.alkaneCriticalVolume = self.calc_alkane_critical_volume()
        self.alkaneSpecificGravity = self.calc_alkane_specific_gravity()
        self.alkaneMolecularWeight = self.calc_alkane_mw()
        self.alkaneCriticalPressure = self.calc_alkane_critical_pressure()
        self.criticalTemperature = self.calc_critical_temperature()
        self.criticalVolume = self.calc_critical_volume()

    def calc_alkane_critical_temp(self):
        tb = self.boilTemp
        tc = tb * (0.533272 + (0.191017e-3 * tb) + (0.779681e-7 * tb ** 2) - (0.284376e-10 * tb ** 3) + (
            0.959468e28 / tb ** 13)) ** (-1)
        return tc

    def calc_alpha(self):
        return 1 - (self.boilTemp / self.alkaneCriticalTemperature)

    def calc_alkane_critical_volume(self):
        vc = (1 - (
             0.419869 - 0.505839 * self.alpha - 1.56436 * self.alpha ** 3 - 9481.70 * self.alpha ** 14)) ** (-8)
        return vc

    def calc_alkane_specific_gravity(self):
        sg = 0.843593 - 0.128624 * self.alpha - 3.36159 * self.alpha ** 3 - 13749.5 * self.alpha ** 12
        return sg

    def objective_function(self, theta):
        a = math.exp((5.71419 + 2.71579*theta - 0.286590*theta**2 - 39.8544/theta - 0.122488/theta**2))
        return a - 24.7522*theta + 35.3155*theta**2 - self.boilTemp

    def calc_alkane_mw(self):
        mw_guess = self.boilTemp/(10.44 - 0.0052*self.boilTemp)
        theta_guess = math.log(mw_guess)
        theta = newton(self.objective_function, theta_guess)
        return math.exp(theta)

    def calc_alkane_critical_pressure(self):
        a = 3.83354
        b = 1.19629
        c = 34.8888
        d = 36.1952
        e = 104.193
        return (a + b*self.alpha**0.5 + c*self.alpha + d*self.alpha**2 + e*self.alpha**4)**2

    def calc_critical_temperature(self):
        a = 0.362456
        b = 0.0398285
        c = 0.948125
        delta_sg_t = math.exp(5*(self.alkaneSpecificGravity - self.specificGravity))-1
        ft = delta_sg_t*(-a/self.boilTemp**0.5 + (b - c/self.boilTemp**0.5)*delta_sg_t)
        tc = self.alkaneCriticalTemperature*((1 + 2*ft)/(1 - 2*ft))**2
        return tc

    def calc_critical_volume(self):
        a = 0.466590
        b = 0.182421
        c = 3.01721
        delta_sg_v = math.exp(4*(self.alkaneSpecificGravity**2 - self.specificGravity**2)) - 1
        fv = delta_sg_v*(a/self.boilTemp**0.5 + (-b + c/self.boilTemp**0.5)*delta_sg_v)
        vc = self.alkaneCriticalVolume*((1 + 2*fv)/(1 - 2*fv))**2
        return vc


if __name__ == '__main__':
    var = Twu(919.34, 1.097)

    print("The Boiling point is {:0.2f} R".format(var.boilTemp))
    print("The specific gravity is {:0.3f}".format(var.specificGravity))
    print("The alkane specific gravity is {:0.3f}".format(var.alkaneSpecificGravity))
    print("The alkane critical temperature is {:0.2f} R".format(var.alkaneCriticalTemperature))
    print("The alkane critical volume {:0.2f} ft3 lb-1 mol-1".format(var.alkaneCriticalVolume))
    print("The alkane critical pressure is {:0.2f} psia".format(var.alkaneCriticalPressure))
    print("The component critical temperature is {:0.2f} R".format(var.criticalTemperature))
    print("The component critical volume is {:0.3f} ft3 lb-1 mol-1".format(var.criticalVolume))
