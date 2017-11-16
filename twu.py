from scipy.optimize import brentq

class Twu:
    def __init__(self, boil_temp, specific_gravity):
        self.boilTemp = boil_temp
        self.specificGravity = specific_gravity
        self.__alkaneCriticalTemperature = self.calc_critical_temp()
        self.__alpha = self.calc_alpha()
        self.__alkaneCriticalVolume = self.calc_critical_volume()
        self.__alkaneSpecificGravity = self.calc_specific_gravity()

    def calc_critical_temp(self):
        tb = self.boilTemp
        tc = tb*(0.533272 + (0.191017e-3*tb) + (0.779681e-7*tb**2) - (0.284376e-10*tb**3) + (0.959468e28/tb**13))**(-1)
        return tc

    def calc_alpha(self):
        return 1 - (self.boilTemp / self.__alkaneCriticalTemperature)

    def calc_critical_volume(self):
        vc = (1 - (0.419869 - 0.505839 * self.__alpha - 1.56436 * self.__alpha ** 3 - 9481.70 * self.__alpha ** 14)) ** (-8)
        return vc

    def calc_specific_gravity(self):
        sg = 0.843593 - 0.128624 * self.__alpha - 3.36159*self.__alpha**3 - 13749.5*self.__alpha**12
        return sg




if __name__ == '__main__':

    var = Twu(919.34, 1.097)