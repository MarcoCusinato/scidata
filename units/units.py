class units:
    """
    Class with the various conversion units and constants in cgs.
    """
    def __init__(self):
        self.solar_mass = 1.98847e33 #g
        self.speed_light = 2.99792458e10 #cm/s
        self.MeV = 1.360218e-6 #erg
        self.km = 1e5 #cm
        self.ms = 1e-3 #s

    def convert_to_solar_masses(self, quantity):
        return quantity / self.solar_mass
    
    def convert_to_grams(self, quantity):
        return quantity * self.solar_mass
    
    def convert_to_grams(self, quantity):
        return quantity * self.solar_mass

    def convert_in_c_units(self, quantity):
        return quantity / self.speed_light

    def convert_to_MeV(self, quantity):
        return quantity / self.MeV

    def convert_to_erg(self, quantity):
        return quantity * self.MeV

    def convert_to_km(self, quantity):
        return quantity / self.km

    def convert_to_cm(self, quantity):
        return quantity * self.km

    def convert_to_ms(self, quantity):
        return quantity / self.ms

    def convert_to_s(self, quantity):
        return quantity * self.ms
