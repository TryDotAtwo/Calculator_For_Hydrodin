import numpy as np

from algorithms import Calculations as calcAlg

# Define global parameters
G = 9.81
PressureIntegral = 0.25
MaxRadius = 2.5 / 1000
MinRadius = 1 / 1000


def promodelirovatParashu_1(arg1, arg2):
    T0 = 300
    KekW = T0
    TL = 127
    SpeedTemp = 1500
    SecInMin = 60
    for i in range(1, T0 - TL):
        Gopa = (T0 - TL) / (1500 / SecInMin)
        arg1.append(T0)
        arg2.append(Gopa)
        T0 = KekW - i


def promodelirovatPoVremeni(arg1, arg2):
    Padla = 0
    for i in range(1, len(arg1)):
        Padla = Padla + arg1[i - 1]
        arg2.append(Padla / i)


def promodelirovatRadiusKapli(Mezytilen, Azot, arg1, arg2):
    Mezytilen.RadiusKapla = MaxRadius
    k = Mezytilen.RadiusKapla

    for i in range(1, round((MaxRadius - MinRadius) * 1000000)):
        F2 = 4 / 3 * 3.14 * pow(Mezytilen.RadiusKapla, 3) * Mezytilen.PlotnostKapla * G
        KekLol = 3.14 * 24 * Azot.DinamicVyazkost * pow(Mezytilen.RadiusKapla, 4) * Azot.KoefTeploprovodnosti * 200 / (
                Azot.PlotnostPara * Azot.TeplotaParoobrazovaniya * pow(50 / pow(10, 6), 4)) * PressureIntegral
        F4 = KekLol - F2
        arg1.append(F4 * 1000)
        arg2.append(Mezytilen.RadiusKapla * 1000)
        Mezytilen.RadiusKapla = k - i / 1000 / 1000


def promodelirovatGavno_1(Mezytilen, Azot, arg1, arg2):
    Mezytilen.RadiusKapla = MaxRadius
    DT = 200
    k = DT
    for i in range(1, 200):
        F2 = 4 / 3 * 3.14 * pow(Mezytilen.RadiusKapla, 3) * Mezytilen.PlotnostKapla * G
        KekLol = 3.14 * 24 * Azot.DinamicVyazkost * pow(Mezytilen.RadiusKapla, 4) * Azot.KoefTeploprovodnosti * DT / (
                Azot.PlotnostPara * Azot.TeplotaParoobrazovaniya * pow(50 / pow(10, 6), 4)) * PressureIntegral
        F4 = KekLol - F2
        arg1.append(F4 * 1000)
        arg2.append(DT)
        DT = k - i


def promodelirovatGavno_2(Mezytilen, Azot, arg1, arg2):
    Mezytilen.RadiusKapla = MaxRadius
    k = Mezytilen.NachalnayaTemperatura
    for i in range(1, 100 * round(Mezytilen.NachalnayaTemperatura - Mezytilen.PlavleniyaTemperatura - 1)):
        RMez = calcAlg.CalcMaxRadius(Azot, Mezytilen)
        ConstantDrop = calcAlg.CalcConstantDrop(Azot, Mezytilen)
        t1 = calcAlg.CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen)
        t2 = calcAlg.CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen)
        t3 = calcAlg.CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen)
        t = t1 + t2 + t3
        arg1.append(t1)
        arg2.append(Mezytilen.NachalnayaTemperatura)
        Mezytilen.NachalnayaTemperatura = k - i / 100


def promodelirovatGavno_3(Mezytilen, Azot, arg1, arg2):
    Mezytilen.RadiusKapla = MaxRadius
    k = Mezytilen.RadiusKapla
    for i in range(1, round((MaxRadius - MinRadius) * 1000000)):
        ConstantDrop = calcAlg.CalcConstantDrop(Azot, Mezytilen)
        t1 = calcAlg.CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen)
        t2 = calcAlg.CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen)
        t3 = calcAlg.CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen)
        t = t1 + t2 + t3
        arg1.append(t1)
        arg2.append(Mezytilen.RadiusKapla * 1000)
        Mezytilen.RadiusKapla = k - i / 1000 / 1000


def promodelirovatGavno_4(Mezytilen, Azot, CountsOfX, Kek228, t):
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
    Padla = 0
    for i in range(1, CountsOfX + 1):
        Padla = Padla + t[i - 1]
        Kek228.append(Padla / i)


