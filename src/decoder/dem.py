import stim
import sinter

def dem(circuit):
    dem = circuit.detector_error_model()
    return dem

class Boomerang:
    def __init__(self, dem):
        self.dem = dem

    def throw(self):
        dem = self.dem
        for line in repr(dem):
            find 
        print("Boomerang thrown!")
        return self.data

    def catch(self):
        print("Boomerang caught!")
        return self.data