from interfaces.IPhysicalParametrized import IPhysicalParametrized


class Sharik(IPhysicalParametrized):
    def getAllPhysicalParameters(self):
        print('')

    def __init__(self, nachalnayaTemperatura, plavleniyaTemperatura, plotnostKapla, teploEmkostKapla, radiusKapla, dolaPoverhnostiKapli, teplotaPlavleniya, teploEmkostKristall):
        self.NachalnayaTemperatura = nachalnayaTemperatura
        self.PlavleniyaTemperatura = plavleniyaTemperatura
        self.PlotnostKapla = plotnostKapla
        self.TeploEmkostKapla = teploEmkostKapla
        self.RadiusKapla = radiusKapla
        self.DolaPoverhnostiKapli = dolaPoverhnostiKapli
        self.TeplotaPlavleniya = teplotaPlavleniya
        self.TeploEmkostKristall = teploEmkostKristall