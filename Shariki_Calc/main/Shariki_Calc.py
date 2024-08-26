from __future__ import annotations
import numpy as np

from algorithms import Calculations as calcAlg
from algorithms import Modeler as modeler
from plotters.Plotter import Plotter as graphPlotter

from physical_objects.Jija import Jija
from physical_objects.Sharik import Sharik

from physical_parametres import PhysicalParamentres as PhysicalParamentres

G = 9.81
Azot = Jija(127, 77.4, 4.5, 7.26 / 1000, 200 * 1000, 8.85 / 1000, 55.2 / 10000000, 811.67, 1151)
Argon = Jija(127, 87.29, 5.707, 6 / 1000, 161.7 * 1000, 11 / 1000, 240 / 10000000, 1392, 1151)
Neon = Jija(100, 27.102, 11.5, 7.2 / 1000, 86 * 1000, 4.76 / 1000, 124 / 10000000, 1207, 1151)
Helium = Jija(20, 4.215, 15.5, 8.98 / 1000, 20.3 * 1000, 0.1 / 1000, 3 / 10000000, 125, 1151)
Hidrogen = Jija(100, 20.38, 1.19, 16.7 / 1000, 452 * 1000, 2.10 / 1000, 13.5 / 10000000, 70.75, 1151)

Mezytilen = Sharik(300, 228.43, 861, 1.75 * 1000, 4 / 2 / 1000, 0.5, 80.12 * 1000, 1.21 * 1000)
Metan = Sharik(111.66, 90.67, 494, 3.37 * 1000, 3.7 / 2 / 1000, 0.5, 57.87 * 1000, 2.71 * 1000)

JijaList = {'N': Azot, 'Ar': Argon, 'Ne': Neon, 'He': Helium, 'H': Hidrogen}

T0Array = []
TimeArray = []
modeler.promodelirovatParashu_1(T0Array, TimeArray)
graphPlotter.plotGraphic_TotT0(T0Array, TimeArray)

PhysicalParamentres.Force.SetValue(0)
print(PhysicalParamentres.Force.GetValue)

ExperimentalTime = [36, 31, 22, 20, 20, 26, 25, 29, 25, 23, 18, 25, 27, 17, 20, 23, 29, 28, 31, 31, 17, 23, 33, 23, 19,
                    23, 25]
Kek228 = []
modeler.promodelirovatPoVremeni(ExperimentalTime, Kek228)
graphPlotter.plotGraphic_TotR(np.linspace(0, len(Kek228), len(Kek228)), Kek228)


ExperimentalTime = [19, 17, 29, 12, 14, 14, 26, 20, 13, 15, 28, 20, 12, 20, 23, 24, 18, 17, 19, 25, 23.5, 25, 17, 14,
                    27, 19, 13]
print(np.average(ExperimentalTime))
Kek228 = []
modeler.promodelirovatPoVremeni(ExperimentalTime, Kek228)
graphPlotter.plotGraphic_TotR(np.linspace(0, len(Kek228), len(Kek228)), Kek228)



CountsOfX = pow(10, 4)

Xot = 3.6 / 1000 / 2
Xdo = 4.2 / 1000 / 2
mu = (Xdo + Xot) / 2
sigma = (Xdo - mu) / 2
RadiusArray = np.random.normal(mu, sigma, CountsOfX)
Mezytilen.RadiusKapla = RadiusArray
ConstantDrop = calcAlg.CalcConstantDrop(Azot, Mezytilen)
t1 = calcAlg.CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen)
t2 = calcAlg.CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen)
t3 = calcAlg.CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen)
t = t1 + t2 + t3
Mezytilen.RadiusKapla = mu
calcAlg.CalcConstantDrop(Azot, Mezytilen)
tKek = (calcAlg.CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen)
        + calcAlg.CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen)
        + calcAlg.CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen))
Average = np.mean(t)

Kek228 = [t[0]]
Padla = 0
for i in range(1, CountsOfX + 1):
    Padla = Padla + t[i - 1]
    Kek228.append(Padla / i)

print(Average)
print(Kek228[CountsOfX])
Kek = Average - Kek228[CountsOfX]
print(Kek)
graphPlotter.plotGraphic_TotR(np.linspace(0, CountsOfX, len(Kek228)), Kek228)

Result = []
Result1 = []
modeler.promodelirovatRadiusKapli(Mezytilen, Azot, Result, Result1)
graphPlotter.plotGraphic_FotR(Result1, Result)

Result = []
Result1 = []
modeler.promodelirovatGavno_1(Mezytilen, Azot, Result, Result1)
graphPlotter.plotGraphic_FotDT(Result1, Result)

# MaxRadius = 2.5/1000
# MinRadius = 1/1000
# Result = []
# Result1 = []
# Result2 = []
# Result3 = []
# DT = 200
# k = DT
# Mezytilen.RadiusKapla = 2.5/1000
# Kek = Mezytilen.RadiusKapla
# for j in range(1, round((MaxRadius - MinRadius)*1000)):
#     for i  in range(1, 200):
#        F2 = 4 / 3 * 3.14 * pow(Mezytilen.RadiusKapla, 3) * Mezytilen.PlotnostKapla * G
#        KekLol = 2 * 3.14 * 24 * Azot.DinamicVyazkost * pow(Mezytilen.RadiusKapla, 4) * Azot.KoefTeploprovodnosti*DT/(Azot.PlotnostPara*Azot.TeplotaParoobrazovaniya*pow(50/pow(10,6),4))*0.125
#        F4 = KekLol - F2
#        Result.append(F4*1000)
#        Result1.append(DT)
#        DT = k - i     
#     Result2.append(Mezytilen.RadiusKapla*1000)
#     plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
#     Mezytilen.RadiusKapla = Kek - j/1000
#     Result1 = []
#     Result = []
# plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
# plt.plot(Result2, Result, marker='o', color='b', linestyle='-')
# plt.xlabel('DT, K')
# plt.ylabel('F, mN')
# plt.title('F ot DT')
# plt.grid(True)
# plt.show()


for i in range(1, 1000):
    P = 24 * 27 * Azot.DinamicVyazkost * calcAlg.CalcConstantDrop(Azot, Mezytilen) / Mezytilen.RadiusKapla * 0.125
    S = 2 * 3.14 * pow(Mezytilen.RadiusKapla, 2)
    F = P * S
    F2 = 4 / 3 * 3.14 * pow(Mezytilen.RadiusKapla, 3) * Mezytilen.PlotnostKapla * G
    F3 = F - F2
    Mezytilen.RadiusKapla = Mezytilen.RadiusKapla - i / 1000000

# CalcNormalGause(mu, sigma, Color, ColorArea, NameOfPlot, NameOfY, NameOfX, NameOfArea, Xot, Xdo, CountsOfX, linewidth)
# annotations = []
# Xot = 3.3
# Xdo = 4.1
# CountsOfX = 1000
# linewidth = 3
# NameOfY = 'Couts of ball, Conventional units'
# NameOfX = 'Diameter, mm'
# NameOfPlot = ''
# CalcNormalGause(3.78, 0.1,'b-' , 'b', NameOfPlot, NameOfY, NameOfX, '1 ball in 7  second', Xot, Xdo, CountsOfX, linewidth)
# CalcNormalGause(3.77, 0.1,'r-' , 'r', NameOfPlot, NameOfY, NameOfX, '1 ball in 10 second', Xot, Xdo, CountsOfX, linewidth)
# CalcNormalGause(3.7, 0.1,'g-' , 'g', NameOfPlot, NameOfY, NameOfX, '1 ball in 30 second', Xot, Xdo, CountsOfX, linewidth)
# plt.show()

# for symbol, Jija in JijaList.items():

#     RadiusMax = CalcMaxRadius(Jija, Mezytilen)
#     # Mezytilen.RadiusKapla = RadiusMax*2
#     ConstantDrop = CalcConstantDrop(Jija, Mezytilen)
#     t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Jija, Mezytilen)
#     t2 = CalcDlitelnostKristall(ConstantDrop, Jija, Mezytilen)
#     t3 = CalcDlitelnostOstivaniya(ConstantDrop, Jija, Mezytilen)
#     t = t1 + t2 + t3
#     print("Jija:", symbol)
#     print("ConstantDrop:", ConstantDrop)
#     print("Dmax, mm:", RadiusMax*2000)
#     print("t, s:", t1)
#     print("t, s:", t2)
#     print("t, s:", t3)
#     print("t, s:", t)

# for symbol, Jija in JijaList.items():

#     RadiusMax = CalcMaxRadius(Jija, Metan)
#     Metan.RadiusKapla = RadiusMax*2
#     ConstantDrop = CalcConstantDrop(Jija, Metan)
#     t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Jija, Metan)
#     t2 = CalcDlitelnostKristall(ConstantDrop, Jija, Metan)
#     t3 = CalcDlitelnostOstivaniya(ConstantDrop, Jija, Metan)
#     t = t1 + t2 + t3
#     print("Jija:", symbol)
#     print("ConstantDrop:", ConstantDrop)
#     print("Dmax, mm:", RadiusMax*2000)
#     print("t, s:", t1)
#     print("t, s:", t2)
#     print("t, s:", t3)
#     print("t, s:", t)


# # # Azot = PropertiesJija(126, 77.4, 4.5, 7.26 / 1000, 200 * 1000, 8.85/1000, 55.2 / 10000000, 811.67, 1151)
# # # Argon = PropertiesJija(136, 87.29, 5.707,  6 / 1000, 161.7 * 1000, 12.50/1000, 240 / 10000000, 1392, 1151)
# # # Neon = PropertiesJija(60, 27.102, 11.5, 7.2 / 1000, 86 * 1000, 4.76/1000, 124 / 10000000, 1207, 1151)
# # # Helium = PropertiesJija(20, 4.215, 15.5, 8.98 / 1000, 20.3 * 1000, 0.09/1000, 3 / 10000000, 125, 1151)
# # # Hidrogen = PropertiesJija(40, 20.38, 1.19, 16.7 / 1000, 452 * 1000, 2.10/1000, 13.5 / 10000000, 70.75, 1151)

# # # Mezytilen = PropertiesSharik(300, 228.43, 861, 1.75 * 1000, 3.5 / 2 / 1000, 0.5, 80.12 * 1000, 1.21 * 1000)
# # # Metan = PropertiesSharik(100.66, 90.67, 494, 3.37 * 1000, 3.5 / 2 / 1000, 0.5, 57.87 * 1000, 2.71 * 1000)

# # # Result = []
# # # Result1 = []

# # # h = CalcH(Azot, Mezytilen)
# # # tV2 = CalcVremyaLeidenV2(h, Azot, Mezytilen)
# # # ConstantDrop = CalcConstantDrop(Azot, Mezytilen)
# # # t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen)
# # # t2 = CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen) 
# # # t3 = CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen)
# # # t = t1 + t2 + t3

# # # print(tV2)
# # # print(t)


Result = []
Result1 = []
modeler.promodelirovatGavno_2(Mezytilen, Azot, Result, Result1)
graphPlotter.plotGraphic_T1otT0(Result1, Result)

Result = []
Result1 = []
modeler.promodelirovatGavno_3(Mezytilen, Azot, Result, Result1)
graphPlotter.plotGraphic_T1otR(Result1, Result)

# Result = []
# Result1 = []
# Azot.LeidenfrostaTemperatura = 149
# k = Azot.LeidenfrostaTemperatura
# for i in range(1, 7000):
#     RMet = CalcMaxRadius(Azot, Mezytilen)
#     Mezytilen.RadiusKapla = RMet*2
#     ConstantDrop = CalcConstantDrop(Azot, Mezytilen)
#     t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen)
#     t2 = CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen)
#     t3 = CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen)
#     t = t1 + t2 + t3
#     Result.append(t)
#     Result1.append(Azot.LeidenfrostaTemperatura)
#     Azot.LeidenfrostaTemperatura = k - i/100
#     graphPlotter.plotGraphic_TotTL_mezitilenInAzot(Result1, Result)


# Result = []
# Result1 = []
# k = Neon.LeidenfrostaTemperatura
# for i in range(1, 7000):
#     RMet = CalcMaxRadius(Neon, Metan)
#     Metan.RadiusKapla = RMet*2
#     ConstantDrop = CalcConstantDrop(Neon, Metan)
#     t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Neon, Metan)
#     t2 = CalcDlitelnostKristall(ConstantDrop, Neon, Metan)
#     t3 = CalcDlitelnostOstivaniya(ConstantDrop, Neon, Metan)
#     t = t1 + t2 + t3
#     Result.append(t)
#     Result1.append(Neon.LeidenfrostaTemperatura)
#     Neon.LeidenfrostaTemperatura = k - i/100
#     graphPlotter.plotGraphic_TotTL_metanInNeon(Result1, Result)
#

# Result = []
# Result1 = []
# k = Hidrogen.LeidenfrostaTemperatura
# for i in range(1, 8000):
#     RMet = CalcMaxRadius(Hidrogen, Metan)
#     Metan.RadiusKapla = RMet*2
#     ConstantDrop = CalcConstantDrop(Hidrogen, Metan)
#     t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Hidrogen, Metan)
#     t2 = CalcDlitelnostKristall(ConstantDrop, Hidrogen, Metan)
#     t3 = CalcDlitelnostOstivaniya(ConstantDrop, Hidrogen, Metan)
#     t = t1 + t2 + t3
#     Result.append(t)
#     Result1.append(Hidrogen.LeidenfrostaTemperatura)
#     Hidrogen.LeidenfrostaTemperatura = k - i/100
#     graphPlotter.plotGraphic_TotTL(Result1, Result)


# k = 3
# for g in range(1, k):
#     PropertiesSharik.NachalnayaTemperatura = 230
#     for i in range(1, 101):

#         ConstantDrop = CalcConstantDrop(PropertiesSharik.PlotnostKapla, PropertiesSharik.TeploEmkostKapla, PropertiesSharik.RadiusKapla, PropertiesSharik.DolaPoverhnostiKapli, PropertiesJija.DinamicVyazkost, PropertiesJija.PlotnostPara, PropertiesJija.KoefTeploprovodnosti, PropertiesJija.TeplotaParoobrazovaniya)
#         t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, PropertiesSharik.NachalnayaTemperatura, PropertiesSharik.PlavleniyaTemperatura, PropertiesJija.JijiTemperatura)
#         t2 = CalcDlitelnostKristall(ConstantDrop, PropertiesSharik.TeploEmkostKapla, PropertiesSharik.TeplotaPlavleniya, PropertiesSharik.PlavleniyaTemperatura, PropertiesJija.JijiTemperatura)
#         t3 = CalcDlitelnostOstivaniya(ConstantDrop, PropertiesSharik.TeploEmkostKapla, PropertiesSharik.TeploEmkostKristall, PropertiesJija.LeidenfrostaTemperatura, PropertiesSharik.PlavleniyaTemperatura, PropertiesJija.JijiTemperatura)

#         t = t1 + t2 + t3

#         Result.append(t)

#         PropertiesSharik.NachalnayaTemperatura = 230 + i/100


#     PropertiesSharik.RadiusKapla =  (3.5+g/1000) / 2 / 1000

#     # if Result[i-1] < 0:
#     #     print(i)
#     Result1 = []
#     Result2 = []
# for i in range (1, 101):
#     Result1.append(Result[i-1])
#     g = i + 99
#     Result2.append(Result[g])

# NachalnayaTemperatura = np.linspace(230, 300, num=len(Result1))
# graphPlotter.plotDependency(Result1, NachalnayaTemperatura)
# graphPlotter.plotTwoGraphics_TotT0(NachalnayaTemperatura, Result1, Result2)


# Пример данных для построения графика
