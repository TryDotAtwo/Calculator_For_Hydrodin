import numpy as np
import math
from scipy.stats import norm
from scipy.integrate import quad
from physical_objects import Jija as jija
from physical_objects import  Sharik as Kapla

G = 9.81
PressureIntegral = 0.25
H = 50 / pow(10, 6)
ConstantStef_Bolth = 5.670367 / pow(10, 8)


def CalcH(Jija, Kapla):
    e = Jija.KoefTeploprovodnosti / Kapla.RadiusKapla
    k = 2 * Kapla.PlotnostKapla * Jija.PlotnostJija * G * pow(Jija.KoefTeploprovodnosti, 3) * Kapla.TeplotaPlavleniya
    ds = 9 * Jija.DinamicVyazkost * Kapla.RadiusKapla * (Kapla.NachalnayaTemperarura - Jija.JijiTemperatura)
    ks = k / ds
    Kek = 1 + Jija.TeploEmkostPara * (Kapla.NachalnayaTemperarura - Jija.JijiTemperatura) / 2 / Kapla.TeplotaPlavleniya
    f = np.power((ks * Kek), 0.25)
    h = e + f
    return h


def CalcVremyaLeidenV2(h, Jija, Kapla):
    a = 2 * Kapla.RadiusKapla * Kapla.PlotnostKapla / 3
    b = (Kapla.TeploEmkostKapla * (
            Kapla.NachalnayaTemperarura - Kapla.PlavleniyaTemperatura) + Kapla.TeplotaPlavleniya) / (
                Kapla.PlavleniyaTemperatura - Jija.JijiTemperatura)
    c = 1.1 * Kapla.TeploEmkostKristall * 0.5 * math.log(abs((Kapla.PlavleniyaTemperatura - Jija.JijiTemperatura) / (
            Jija.LeidenfrostaTemperarura - Jija.JijiTemperatura)))
    t = a * (b + c) / h
    return t


def CalcMaxRadius(Jija, Kapla):
    Lc = np.power(Jija.PoverhNat / (Jija.PlotnostJija * G), 0.5)
    R = Lc * np.power((3 * Jija.PlotnostJija) / (2 * Kapla.PlotnostKapla), 0.5)
    return R


def CalcConstantDrop(Jija, Kapla):
    a = 4 * Kapla.PlotnostKapla * Kapla.TeploEmkostKapla * Kapla.RadiusKapla / (3 * Kapla.DolaPoverhnostiKapli)
    b = pow(9 * Jija.DinamicVyazkost * Kapla.RadiusKapla / (
            2 * Kapla.PlotnostKapla * Jija.PlotnostPara * G * pow(Jija.KoefTeploprovodnosti,
                                                                  3) * Jija.TeplotaParoobrazovaniya), 0.25)
    ConstantDrop = a * b
    return ConstantDrop


def CalcDlitelnostOhlagdeniya(Jija, Kapla):
    a = pow((Kapla.NachalnayaTemperatura - Jija.JijiTemperatura), 0.25)
    b = pow((Kapla.PlavleniyaTemperatura - Jija.JijiTemperatura), 0.25)
    DlitelnostOhlagdeniya = CalcConstantDrop(Jija, Kapla) * (a - b)
    return DlitelnostOhlagdeniya


def CalcDlitelnostKristall(Jija, Kapla):
    a = pow((Kapla.PlavleniyaTemperatura - Jija.JijiTemperatura), -0.75)
    DlitelnostKristall = CalcConstantDrop(Jija, Kapla) * Kapla.TeplotaPlavleniya * a / Kapla.TeploEmkostKapla
    return DlitelnostKristall


def CalcDlitelnostOstivaniya(Jija, Kapla):
    a = pow((Kapla.PlavleniyaTemperatura - Jija.JijiTemperatura), 0.25)
    b = pow((Jija.LeidenfrostaTemperatura - Jija.JijiTemperatura), 0.25)
    c = Kapla.TeploEmkostKristall / Kapla.TeploEmkostKapla
    DlitelnostOstivaniya = CalcConstantDrop(Jija, Kapla) * c * (a - b)
    return DlitelnostOstivaniya

def SummTime(Jija, Kapla):
    a = CalcDlitelnostOhlagdeniya(Jija, Kapla)
    b = CalcDlitelnostKristall(Jija, Kapla)
    c = CalcDlitelnostOstivaniya(Jija, Kapla)
    T = a + b + c
    return T

def CalcHV2(Jija, Kapla):
    return 0

def CalcHWave(Jija, Kapla):
    H = pow(Jija.DinamicVyazkost*Kapla.RadiusKapla*ConstantStef_Bolth*(pow(Kapla.NachalnayaTemperarura, 4)-pow(Jija.JijiTemperarura, 4))/(Jija.PlotnostJija*Kapla.PlotnostKapla*G*Jija.TeplotaParoobrazovaniya), 1/3)
    return H

def CalcHTipa(Jija, Kapla):
    H = pow(Jija.DinamicVyazkost * Jija.KoefTeploprovodnosti * Kapla.RadiusKapla * (Kapla.NachalnayaTemperarura - Jija.JijiTemperarura)/(Jija.PlotnostJija*Kapla.PlotnostKapla*G*Jija.TeplotaParoobrazovaniya), 0.25)
    return H

def SquareBall(Jija, Kapla):
    Square = 4 * np.pi * pow(Kapla.RadiusKapla, 2)
    return Square

def SquareBallInJija(Jija, Kapla):
    Square = 4 * np.pi * pow(Kapla.RadiusKapla, 2)/2
    return Square

def Volume(Jija, Kapla):
    Volume = 4/3 * np.pi * pow(Kapla.RadiusKapla, 3)
    return Volume

def ForceGravity(Jija, Kapla):
    ForceGravity = Volume(Jija, Kapla) * Kapla.PlotnostKapla * G
    return ForceGravity

def ForcePressure(Jija, Kapla):
    ForcePressure =  Pressure(Jija, Kapla)  * SquareBallInJija(Jija, Kapla) * PressureIntegral
    return ForcePressure

def ForcePressureHconst(Jija, Kapla, Hconst):
    ForcePressure = PressureHconst(Jija, Kapla, Hconst) * PressureIntegral
    return ForcePressure

def Pressure(Jija, Kapla):
    Pressure = 24 * Jija.DinamicVyazkost * pow(Kapla.RadiusKapla, 2) * Jija.KoefTeploprovodnosti * (Kapla.NachalnayaTemperatura - Jija.JijiTemperatura) / (
                Jija.PlotnostPara * Jija.TeplotaParoobrazovaniya * pow(CalcHV2(Jija, Kapla), 4))
    return Pressure

def PressureHconst(Jija, Kapla, Hconst):
    Pressure = 24 * Jija.DinamicVyazkost * pow(Kapla.RadiusKapla, 2) * Jija.KoefTeploprovodnosti * (Kapla.NachalnayaTemperatura - Jija.JijiTemperatura) / (
                Jija.PlotnostPara * Jija.TeplotaParoobrazovaniya * pow(Hconst, 4))
    return Pressure

def CalcNormalGause(mu, sigma, Color, ColorArea, NameOfPlot, NameOfY, NameOfX, NameOfArea, Xot, Xdo, CountsOfX, linewidth):
    # Генерируем данные для графика
    z = np.linspace(Xot, Xdo, CountsOfX)
    k = norm.pdf(z, mu, sigma)
    NormalCoef = np.max(k)
    x = np.linspace(Xot, Xdo, CountsOfX)
    y = norm.pdf(x, mu, sigma) / NormalCoef

    # Задаем интервалы для закраски
    x_fill_left = np.linspace(Xot, 3.5, 100)
    y_fill_left = norm.pdf(x_fill_left, mu, sigma) / NormalCoef

    x_fill_right = np.linspace(3.9, Xdo, 100)
    y_fill_right = norm.pdf(x_fill_right, mu, sigma) / NormalCoef

    # Вычисляем площадь закрашенных областей в процентах
    area_left, _ = quad(lambda x: norm.pdf(x, mu, sigma), -np.inf, 3.5)
    area_right, _ = quad(lambda x: norm.pdf(x, mu, sigma), 3.9, np.inf)
    percentage_left = area_left * 100
    percentage_right = area_right * 100
    percentage = percentage_left + percentage_right

    # annotations.append(f"{NameOfArea}: {percentage:.2f}%")
