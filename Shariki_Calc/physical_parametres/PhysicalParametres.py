from algorithms import Calculations as calcAlg

class PhysicalParametrFactory:
    def __init__(self, large_name, small_name, units, values):
        self.large_name = large_name
        self.small_name = small_name
        self.units = units
        self.values = values

    def add_new_method(self, Name, func):
        setattr(self, Name, func)
    
    def SetValue(self, values):
        self.values = values
    
    def GetValue(self):
        return self.values
        
# Начинаем создавать физичекие параметры
InitialMassive = [] #Заполняем параметр значение

# Создаём объект сила, в котором существует метод для расчета силы тяжести шарика и силы давления паровой подушки
Force = PhysicalParametrFactory("Force", "F", "µN", InitialMassive)
Force.add_new_method('Gravity', calcAlg.ForceGravity)
Force.add_new_method('Pressure', calcAlg.ForcePressure)
Force.add_new_method('PressureHconst', calcAlg.ForcePressureHconst)

# Создаём объект время, в котором существует метод для расчета длительности остывания шарика, кристаллизации и остывания шарика, а также их суммы
Time = PhysicalParametrFactory("Time", "t", "s", InitialMassive)
Time.add_new_method('1Stage', calcAlg.CalcDlitelnostOhlagdeniya)
Time.add_new_method('2Stage', calcAlg.CalcDlitelnostKristall)
Time.add_new_method('3Stage', calcAlg.CalcDlitelnostOstivaniya)
Time.add_new_method('Summary', calcAlg.SummTime)

# Создаём объект давление, в котором существует метод для расчета давления пара на шарик
Pressure = PhysicalParametrFactory("Pressure", "Pa", "mPa", InitialMassive)
Pressure.add_new_method('Calc', calcAlg.Pressure)
Pressure.add_new_method('Hconst', calcAlg.Pressure)

# Создаём объект длина, в котором существует метод для расчета длительности остывания шарика, кристаллизации и остывания шарика, а также их суммы
Height = PhysicalParametrFactory("Height", "H", "µm", InitialMassive)
Height.add_new_method('1Stage', calcAlg.CalcH)
Height.add_new_method('2Stage', calcAlg.CalcHWave)
Height.add_new_method('3Stage', calcAlg.CalcHTipa)
Height.add_new_method('Summary', calcAlg.CalcHV2)

# Создаём объект радиус, в котором существует метод для расчета максимального радиуса шарика
Radius = PhysicalParametrFactory("MaxRadius", "R", "mm", InitialMassive)
Radius.add_new_method('Calc', calcAlg.CalcMaxRadius)

# Создаём объект площадь, в котором существует метод для расчета максимального радиуса шарика
Square = PhysicalParametrFactory("Square", "S", "mm²", InitialMassive)
Square.add_new_method('Ball', calcAlg.SquareBall)
Square.add_new_method('BallInJij', calcAlg.SquareBallInJija)

# Создаём объект Radius, в котором существует метод для расчета максимального радиуса шарика
Volume = PhysicalParametrFactory("Volume", "V", "mm³", InitialMassive)
Volume.add_new_method('Calc', calcAlg.Volume)

# Создаём объект температура
Temperature = PhysicalParametrFactory("Temperature", "T", "K", InitialMassive)