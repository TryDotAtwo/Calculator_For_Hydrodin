from interfaces.IPhysicalParametrized import IPhysicalParametrized


class Jija(IPhysicalParametrized):
    def getAllPhysicalParameters(self):
        paramsDict = {'Leidenfrost temperature': self.LeidenfrostaTemperatura,
                      'Jija temperature': self.JijiTemperatura,
                      'Plotnost para': self.PlotnostPara,
                      'Koef teplo provodimosti': self.KoefTeploprovodnosti,
                      'teplota paroobrazovania': self.TeplotaParoobrazovaniya,
                      'Poverhn natiagenie': self.PoverhNat
                      }


    def __init__(self, leidenfrostaTemperatura, jijiTemperatura, plotnostPara, koefTeploprovodnosti, teplotaParoobrazovaniya, poverhNat, dinamicVyazkost, plotnostJija, TeploEmkostPara):
        self.LeidenfrostaTemperatura = leidenfrostaTemperatura
        self.JijiTemperatura = jijiTemperatura
        self.PlotnostPara = plotnostPara
        self.KoefTeploprovodnosti = koefTeploprovodnosti
        self.TeplotaParoobrazovaniya = teplotaParoobrazovaniya
        self.PoverhNat = poverhNat
        self.DinamicVyazkost = dinamicVyazkost
        self.PlotnostJija = plotnostJija
        self.TeploEmkostPara = TeploEmkostPara
