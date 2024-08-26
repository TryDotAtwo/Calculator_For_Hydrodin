import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad
from scipy.stats import norm

G = 9.81
ConstantStef_Bolth = 5.670367 / pow(10, 8)

# def plot_dependency(Result, NachalnayaTemperarura):
#     plt.plot(NachalnayaTemperarura, Result, marker='o', color='b', linestyle='-')
#     plt.xlabel('T0')
#     plt.ylabel('t')
#     plt.title('t ot T0')
#     plt.grid(True)
#     # plt.show()

def CalcNormalGause(mu, sigma, Color, ColorArea, NameOfPlot, NameOfY, NameOfX, NameOfArea, Xot, Xdo, CountsOfX, linewidth):
    # Генерируем данные для графика 
    z = np.linspace(Xot, Xdo, CountsOfX)
    k = norm.pdf(z, mu, sigma)
    NormalCoef = np.max(k)
    x = np.linspace(Xot, Xdo, CountsOfX)
    y = norm.pdf(x, mu, sigma)/NormalCoef

    # Создаем график нормального распределения
    plt.plot(x, y, Color, label=NameOfPlot, linewidth=linewidth)
    
    # Задаем интервалы для закраски
    x_fill_left = np.linspace(Xot, 3.5, 100)
    y_fill_left = norm.pdf(x_fill_left, mu, sigma)/NormalCoef

    x_fill_right = np.linspace(3.9, Xdo, 100)
    y_fill_right = norm.pdf(x_fill_right, mu, sigma)/NormalCoef

    # Закрашиваем указанные интервалы
    plt.fill_between(x_fill_left, y_fill_left, color=ColorArea)
    plt.fill_between(x_fill_right, y_fill_right, color=ColorArea)

    # Вычисляем площадь закрашенных областей в процентах
    area_left, _ = quad(lambda x: norm.pdf(x, mu,sigma), -np.inf, 3.5)
    area_right, _ = quad(lambda x: norm.pdf(x, mu,sigma), 3.9, np.inf)
    percentage_left = area_left * 100
    percentage_right = area_right * 100
    percentage = percentage_left + percentage_right
    
    
    # annotations.append(f"{NameOfArea}: {percentage:.2f}%")
    
    # Выводим аннотации справа вверху окна
    for i, text in enumerate(annotations):
        plt.annotate(text, xy=(1, 1 - 0.08*i), xycoords='axes fraction', fontsize=14, color=ColorArea, ha='right', va='top')
        
    plt.title(NameOfPlot)
    plt.xlabel(NameOfX)
    plt.ylabel(NameOfY)
    plt.grid(True)

    

def CalcH(Jija, Kapla):
    e = Jija.KoefTeploprovodnosti/Kapla.RadiusKapla
    k = 2*Kapla.PlotnostKapla*Jija.PlotnostJija*G*pow(Jija.KoefTeploprovodnosti,3)*Kapla.TeplotaPlavleniya
    ds = 9*Jija.DinamicVyazkost*Kapla.RadiusKapla*(Kapla.NachalnayaTemperarura-Jija.JijiTemperarura)
    ks = k/ds
    Kek = 1 + Jija.TeploEmkostPara*(Kapla.NachalnayaTemperarura-Jija.JijiTemperarura)/2/Kapla.TeplotaPlavleniya
    f = np.power((ks*Kek),0.25)
    h = e + f   
    return h

def CalcVremyaLeidenV2(h, Jija, Kapla):
    a = 2 * Kapla.RadiusKapla * Kapla.PlotnostKapla/3
    b = (Kapla.TeploEmkostKapla * (Kapla.NachalnayaTemperarura - Kapla.PlavleniyaTemperarura) + Kapla.TeplotaPlavleniya)/(Kapla.PlavleniyaTemperarura - Jija.JijiTemperarura)
    c = 1.1*Kapla.TeploEmkostKristall*0.5*math.log(abs((Kapla.PlavleniyaTemperarura-Jija.JijiTemperarura)/(Jija.LeidenfrostaTemperarura-Jija.JijiTemperarura)))
    t = a * (b + c)/h
    return t

def CalcMaxRadius(Jija, Kapla):
    Lc = np.power(Jija.PoverhNat / (Jija.PlotnostJija * G), 0.5)
    R = Lc * np.power((3 * Jija.PlotnostJija) / (2 * Kapla.PlotnostKapla), 0.5)
    return R

def CalcMaxRadiusV2(Jija, Kapla):
    R = np.power(3*Jija.PoverhNat /(G*(2*Kapla.PlotnostKapla-Jija.PlotnostJija)), 0.5)
    return R


def CalcMaxRadiusSferoid(Jija, Kapla, dR):
    R = dR + pow(pow(dR,2)+12*Jija.PoverhNat /(G*(2*Kapla.PlotnostKapla)), 0.5)
    return R/2

def CalcMaxRadiusSferoidV2(Jija, Kapla, dR):
    R = dR + pow(pow(dR,2)+12*Jija.PoverhNat /(G*(2*Kapla.PlotnostKapla-Jija.PlotnostJija)), 0.5)
    return R/2

def CalcConstantDrop(Jija, Kapla):
    a = 4*Kapla.PlotnostKapla*Kapla.TeploEmkostKapla*Kapla.RadiusKapla/(3*Kapla.DolaPoverhnostiKapli)
    b = pow(9*Jija.DinamicVyazkost*Kapla.RadiusKapla/(2*Kapla.PlotnostKapla*Jija.PlotnostPara*G*pow(Jija.KoefTeploprovodnosti, 3)*Jija.TeplotaParoobrazovaniya), 0.25)
    ConstantDrop = a*b
    return ConstantDrop

def CalcDlitelnostOhlagdeniya(ConstantDrop, Jija, Kapla):
    a = pow((Kapla.NachalnayaTemperarura-Jija.JijiTemperarura), 0.25)
    b = pow((Kapla.PlavleniyaTemperarura-Jija.JijiTemperarura), 0.25)
    DlitelnostOhlagdeniya = ConstantDrop * (a-b)
    return DlitelnostOhlagdeniya

def CalcDlitelnostKristall(ConstantDrop, Jija, Kapla):
    a = pow((Kapla.PlavleniyaTemperarura-Jija.JijiTemperarura), -0.75)
    DlitelnostKristall = ConstantDrop * Kapla.TeplotaPlavleniya * a / Kapla.TeploEmkostKapla
    return DlitelnostKristall

def CalcDlitelnostOstivaniya(ConstantDrop, Jija, Kapla):
    a = pow((Kapla.PlavleniyaTemperarura-Jija.JijiTemperarura), 0.25)
    b = pow((Jija.LeidenfrostaTemperarura-Jija.JijiTemperarura), 0.25)
    c = Kapla.TeploEmkostKristall/Kapla.TeploEmkostKapla
    DlitelnostOstivaniya = ConstantDrop * c * (a-b)
    return DlitelnostOstivaniya



def Diffur(Jija, Kapla, N, dt, H0):
    h = np.zeros(N+2)
    T = np.zeros(N+2)
    t = np.zeros(N+2)
    t[0] = 0
    t[1] = 0
    h[0] = H0
    h[1] = H0
    T[0] = Kapla.NachalnayaTemperarura - Jija.JijiTemperarura
    T[1] = T[0]
    for i in range(2, N+2):
        T[i] = T[i-1] - (3 * Jija.KoefTeploprovodnosti * T[i-1])/(2 * Kapla.RadiusKapla * h[i-1])*dt
        a = 9*Jija.DinamicVyazkost*Kapla.RadiusKapla*Jija.KoefTeploprovodnosti*T[i]
        b = Jija.PlotnostPara*Kapla.PlotnostKapla*Jija.TeplotaParoobrazovaniya*np.power(h[i-1],4)
        c = 9*Jija.DinamicVyazkost*Kapla.RadiusKapla
        d = Kapla.PlotnostKapla*np.power(h[i-1],2)
        h[i] = (2*h[i-1]-h[i-2]+(a/b)*np.power(dt,2)-(c*d)*dt-G*np.power(dt,2))/(1-c/(d*h[i-1]))
        t[i] = t[i-1] + dt
        if (T[i] < Jija.LeidenfrostaTemperarura):
            break
        if (h[i] < pow(10, -9)) or (T[i] == pow(10, -9)):
            print('i = ' + str(i) + ' h = ' + str(h[i]) + ' T = ' + str(T[i]))
            break
        
    return h, T, t

def Diffur2(Jija, Kapla, N, dt, H0):
    h = np.zeros(N+2)
    T = np.zeros(N+2)
    t = np.zeros(N+2)
    t[0] = 0
    t[1] = 0
    h[0] = H0
    h[1] = H0
    T[0] = Kapla.NachalnayaTemperarura - Jija.JijiTemperarura
    T[1] = T[0]
    for i in range(2, N+2):
        T[i] = T[i-1] - (3 * Jija.KoefTeploprovodnosti * T[i-1])/(2 * Kapla.RadiusKapla * h[i-1]*Kapla.PlotnostKapla*Kapla.TeploEmkostKapla)*dt
        a = 9*Jija.DinamicVyazkost*Kapla.RadiusKapla*Jija.KoefTeploprovodnosti*T[i-1]
        b = Jija.PlotnostPara*Kapla.PlotnostKapla*Jija.TeplotaParoobrazovaniya*np.power(h[i-1],4)
        h[i] = 2*h[i-1]-h[i-2]+((a/b)-G)*np.power(dt,2)
        t[i] = t[i-1] + dt
        if (T[i] < Jija.LeidenfrostaTemperarura):
            break
        if (h[i] < pow(10, -12)) or (T[i] < pow(10, -12)):
            print('i = ' + str(i) + ' h = ' + str(h[i]) + ' T = ' + str(T[i]))
            break
        
    return h, T, t


def Diffur3(Jija, Kapla, N, dt, H0, v0):
    h = np.zeros(N+2)
    T = np.zeros(N+2)
    t = np.zeros(N+2)
    a = np.zeros(N+2)
    v = np.zeros(N+2)
    t[0] = 0
    t[1] = 0
    a[0] = 0
    a[1] = 0
    v[0] = v0
    v[1] = v0
    h[0] = H0
    h[1] = H0
    T[0] = Kapla.NachalnayaTemperarura - Jija.JijiTemperarura
    T[1] = T[0]
    for i in range(2, N+2):
        T[i] = T[i-1] - (3 * Jija.KoefTeploprovodnosti * T[i-1])/(2 * Kapla.RadiusKapla * h[i-1]*Kapla.PlotnostKapla*Kapla.TeploEmkostKapla)*dt
        a[i] = 0 - G + (9*Jija.DinamicVyazkost*Kapla.RadiusKapla*Jija.KoefTeploprovodnosti*T[i-1])/(Jija.PlotnostPara*Kapla.PlotnostKapla*Jija.TeplotaParoobrazovaniya*np.power(h[i-1],4))
        v[i] = v[i-1] + 1.5 * a[i-1] * dt - 0.5 * a[i-2] * dt
        h[i] = h[i-1] + 1.5 * v[i-1] * dt - 0.5 * v[i-2] * dt
        t[i] = t[i-1] + dt
        if (T[i] < Jija.LeidenfrostaTemperarura):
            break
        if (h[i] < pow(10, -12)) or (T[i] < pow(10, -12)):
            print('i = ' + str(i) + ' h = ' + str(h[i]) + ' T = ' + str(T[i]))
            break
        
    return h, T, t, v, a

def DiffurHoin(Jija, Kapla, N, dt, H0, v0):
    h = np.zeros(N+2)
    T = np.zeros(N+2)
    t = np.zeros(N+2)
    ac = np.zeros(N+2)
    v = np.zeros(N+2)
    t[0] = 0
    t[1] = 0
    ac[0] = 0
    ac[1] = 0
    v[0] = v0
    v[1] = v0
    h[0] = H0
    h[1] = H0
    T[0] = Kapla.NachalnayaTemperarura - Jija.JijiTemperarura
    T[1] = T[0]
    
    def a(Kek, Zalupek):
         return 0 - G + (9*Jija.DinamicVyazkost*Kapla.RadiusKapla*Jija.KoefTeploprovodnosti*Zalupek)/(Jija.PlotnostPara*Kapla.PlotnostKapla*Jija.TeplotaParoobrazovaniya*np.power(Kek,4))
   
    for i in range(2, N+2):
        T[i] = T[i-1] - (3 * Jija.KoefTeploprovodnosti * T[i-1])/(2 * Kapla.RadiusKapla * h[i-1]*Kapla.PlotnostKapla*Kapla.TeploEmkostKapla)*dt
 
        ac[i] = 0 - G + (9*Jija.DinamicVyazkost*Kapla.RadiusKapla*Jija.KoefTeploprovodnosti*T[i-1])/(Jija.PlotnostPara*Kapla.PlotnostKapla*Jija.TeplotaParoobrazovaniya*np.power(h[i-1],4))
        v[i] = v[i-1] + 0.5 * dt * (a(h[i-1], T[i-1]) + a(h[i-1] + dt * a(h[i-1], T[i-1]), T[i-1]))
        h[i] = h[i-1] + v[i-1] * dt

        t[i] = t[i-1] + dt
        if (T[i] < (Kapla.PlavleniyaTemperarura-Jija.JijiTemperarura)):
            break
        if (h[i] < pow(10, -12)) or (T[i] < pow(10, -12)):
            print('i = ' + str(i) + ' h = ' + str(h[i]) + ' T = ' + str(T[i]))
            break
        
    return h, T, t, v, ac

class PropertiesSharik:
    def __init__(self, nachalnaya_temperatura, plavleniya_temperatura, plotnost_kapla, teplo_emkost_kapla, radius_kapla, dola_poverhnosti_kapli, teplota_plavleniya, teplo_emkost_kristall, teplo_provodnost):
        self.NachalnayaTemperarura = nachalnaya_temperatura
        self.PlavleniyaTemperarura = plavleniya_temperatura
        self.PlotnostKapla = plotnost_kapla
        self.TeploEmkostKapla = teplo_emkost_kapla
        self.RadiusKapla = radius_kapla
        self.DolaPoverhnostiKapli = dola_poverhnosti_kapli
        self.TeplotaPlavleniya = teplota_plavleniya
        self.TeploEmkostKristall = teplo_emkost_kristall
        self.TeploProvodnost = teplo_provodnost

    


class PropertiesJija:
    def __init__(self, LeidenfrostaTemperarura, JijiTemperarura, PlotnostPara, KoefTeploprovodnosti, TeplotaParoobrazovaniya, PoverhNat, DinamicVyazkost, PlotnostJija, TeploEmkostPara):
        self.LeidenfrostaTemperarura = LeidenfrostaTemperarura
        self.JijiTemperarura = JijiTemperarura
        self.PlotnostPara = PlotnostPara
        self.KoefTeploprovodnosti = KoefTeploprovodnosti
        self.TeplotaParoobrazovaniya = TeplotaParoobrazovaniya
        self.PoverhNat = PoverhNat
        self.DinamicVyazkost = DinamicVyazkost
        self.PlotnostJija = PlotnostJija
        self.TeploEmkostPara = TeploEmkostPara


    


Azot = PropertiesJija(127, 77.4, 4.5, 7.26 / 1000, 200 * 1000, 8.85/1000, 55.2 / 10000000, 811.67, 1151)
Argon = PropertiesJija(127, 87.29, 5.707,  6 / 1000, 161.7 * 1000, 11/1000, 240 / 10000000, 1392, 1151)
Neon = PropertiesJija(100, 27.102, 11.5, 7.2 / 1000, 86 * 1000, 4.76/1000, 124 / 10000000, 1207, 1151)
Helium = PropertiesJija(20, 4.215, 15.5, 8.98 / 1000, 20.3 * 1000, 0.1/1000, 3 / 10000000, 125, 1151)
Hidrogen = PropertiesJija(100, 20.38, 1.19, 16.7 / 1000, 452 * 1000, 2.10/1000, 13.5 / 10000000, 70.75, 1151)

Mezytilen = PropertiesSharik(300, 228.43, 861, 1.75 * 1000, 3.75 / 2 / 1000, 0.5, 80.12 * 1000, 1.21 * 1000, 14.7)
Metan = PropertiesSharik(111.66, 90.67, 494, 3.37 * 1000, 3.7 / 2 / 1000, 0.5, 57.87 * 1000, 2.71 * 1000, 1)

JijaList = {'N': Azot, 'Ar': Argon, 'Ne': Neon, 'He': Helium, 'H': Hidrogen}

KalcAkk = pow(Mezytilen.TeploEmkostKapla*Mezytilen.PlotnostKapla*Mezytilen.TeploProvodnost, 0.5)/pow(10,4)
print(KalcAkk)
KalcAkk = pow(3800*200*8900, 0.5)/pow(10,4)
print(KalcAkk)
KalcAkk = pow(90*160*7800, 0.5)/pow(10,4)
print(KalcAkk)
KalcAkk = 0.32
print(KalcAkk)
Mezytilen.RadiusKapla = (3.9+3.6)/4/1000
ConstantDrop = CalcConstantDrop(Azot, Mezytilen)
t2 = CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen)
print(t2)

Result = []
Result1 = []
MaxRadius = 2.5/1000
MinRadius = 1/1000
Mezytilen.RadiusKapla = MaxRadius
k = Mezytilen.RadiusKapla
for i in range(1, round((MaxRadius - MinRadius)*1000000)):
    ConstantDrop = CalcConstantDrop(Azot, Mezytilen)
    t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen)
    t2 = CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen)
    t3 = CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen)
    t = t1 + t2 + t3
    Result.append(t2)
    Result1.append(Mezytilen.RadiusKapla*1000)
    Mezytilen.RadiusKapla = k - i/1000/1000
plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
plt.xlabel('R, mm')
plt.ylabel('t2, s')
plt.title('t2 ot R')
plt.grid(True)
plt.show()







dR = 0.25/1000
K = CalcMaxRadius(Azot, Mezytilen)*1000
S = CalcMaxRadiusV2(Azot, Mezytilen)*1000
K1 = CalcMaxRadiusSferoid(Azot, Mezytilen, dR)*1000
S1 = CalcMaxRadiusSferoidV2(Azot, Mezytilen, dR)*1000
print(K)
print(S)
print(K1)
print(S1)


t = CalcDlitelnostOhlagdeniya(CalcConstantDrop(Azot, Mezytilen), Azot, Mezytilen)
print(t)
print()
# Diffur(Jija, Kapla, N, dt, H0)
v0 = 1
DT = Mezytilen.NachalnayaTemperarura - Azot.JijiTemperarura
S = 1000*pow(10,-3)
Suka = pow(9*Azot.DinamicVyazkost*Azot.KoefTeploprovodnosti*Mezytilen.RadiusKapla*DT/(Mezytilen.PlotnostKapla*G*Azot.PlotnostPara*Azot.TeplotaParoobrazovaniya), 0.25)
# h, T, t, v, a = Diffur3(Azot, Mezytilen, pow(10, 6), pow(10,-6), Suka, v0)
h, T, t, v, a = DiffurHoin(Azot, Mezytilen, 3*pow(10, 6), pow(10,-9), Suka, v0)
# plt.plot(np.arange(0, len(h)), h, marker='o', color='g', linestyle='-')
# plt.xlabel('i')
# plt.ylabel('h, m')
# plt.title('h ot i')
# plt.grid(True)
# plt.show()

screen_width = 1980
screen_height = 1080

Window_width = 5
Window_height = 4

x1 = 0
x2 = int((screen_width - 700) / 2)
x3 = screen_width - 700
y = int((screen_height - 700) / 2)

# plt.plot(t, T, marker='o', color='g', linestyle='-')
# plt.xlabel('t, s')
# plt.ylabel('T, K')
# plt.title('T ot t')
# plt.grid(True)
# plt.show()


fig1 = plt.figure(figsize=(Window_width, Window_height))
mngr1 = plt.get_current_fig_manager()
mngr1.window.wm_geometry("+{}+{}".format(x1, y))

plt.plot(t, h, marker='o', color='g', linestyle='-')
plt.xlabel('t, s')
plt.ylabel('h, m')
plt.title('h ot t')
plt.grid(True)

# fig2 = plt.figure(figsize=(Window_width, Window_height))
# mngr2 = plt.get_current_fig_manager()
# mngr2.window.wm_geometry("+{}+{}".format(x2, y))

# plt.plot(t, v, marker='o', color='g', linestyle='-')
# plt.xlabel('t, s')
# plt.ylabel('v, m/s')
# plt.title('v ot t')


# fig3 = plt.figure(figsize=(Window_width, Window_height))
# mngr3 = plt.get_current_fig_manager()
# mngr3.window.wm_geometry("+{}+{}".format(x3, y))

# plt.plot(t, a, marker='o', color='g', linestyle='-')
# plt.xlabel('t, s')
# plt.ylabel('a, m/s/s')
# plt.title('a ot t')
# plt.grid(True)


plt.show()












K = CalcMaxRadius(Azot, Mezytilen)*1000
S = CalcMaxRadiusV2(Azot, Mezytilen)*1000
print(K)
print(S)
print(K-S)
dR = 1/1000

Result = []
Result1 = []
k = dR
for i in range(1, 1000):
      Kek = CalcMaxRadiusSferoid(Azot, Mezytilen, dR)*1000
      Result.append(Kek)
      Result1.append(dR*1000)
      dR = k - i/1000/1000
plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
plt.xlabel('dR, mm')
plt.ylabel('R, mm')
plt.title('R ot dR')
plt.grid(True)
plt.show()

Result = []
Result1 = []
DT = Mezytilen.NachalnayaTemperarura - Azot.JijiTemperarura
MaxR = Mezytilen.RadiusKapla
MinR = Mezytilen.RadiusKapla/4
k = MaxR
for i in range(1, round((MaxR - MinR)*1000*1000)):
      Suka = pow(9*Azot.DinamicVyazkost*Azot.KoefTeploprovodnosti*Mezytilen.RadiusKapla*DT/(Mezytilen.PlotnostKapla*G*Azot.PlotnostPara*Azot.TeplotaParoobrazovaniya), 0.25)
      Result.append(Suka*pow(10,6))
      Result1.append(Mezytilen.RadiusKapla*1000)
      Mezytilen.RadiusKapla = k - i/1000/1000
plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
plt.xlabel('R, mm')
plt.ylabel('h, mkm')
plt.title('h ot R')
plt.grid(True)
plt.show()

DT = Mezytilen.NachalnayaTemperarura - Azot.JijiTemperarura
Result = []
Result1 = []
k = DT
for i in range(1, round(Mezytilen.NachalnayaTemperarura - Azot.JijiTemperarura)*90):
   Suka = pow(9*Azot.DinamicVyazkost*Azot.KoefTeploprovodnosti*Mezytilen.RadiusKapla*DT/(Mezytilen.PlotnostKapla*G*Azot.PlotnostPara*Azot.TeplotaParoobrazovaniya), 0.25)
   Result.append(Suka*pow(10,6))
   Result1.append(DT)
   DT = k - i/100
plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
plt.xlabel('DT, K')
plt.ylabel('h, mkm')
plt.title('h ot DT')
plt.grid(True)
plt.show()

Suka = Mezytilen.PlotnostKapla*Mezytilen.TeploEmkostKapla*30/pow(10,6)*Mezytilen.RadiusKapla*2/3/Azot.KoefTeploprovodnosti*np.log(222-48)
print(Suka)

Result = []
Result1 = []
DT = 200
k = DT
for i in range(1, 200):
   F2 = 4 / 3 * 3.14 * pow(Mezytilen.RadiusKapla, 3) * Mezytilen.PlotnostKapla * G
   KekLol = 2 * 3.14 * 24 * Azot.DinamicVyazkost * pow(Mezytilen.RadiusKapla, 4) * Azot.KoefTeploprovodnosti*DT/(Azot.PlotnostPara*Azot.TeplotaParoobrazovaniya*pow(38/pow(10,6),4))*0.25
   F4 = KekLol - F2
   Result.append(F4*1000)
   Result1.append(DT)
   DT = k - i
plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
plt.xlabel('DT, K')
plt.ylabel('F, mN')
plt.title('F ot DT')
plt.grid(True)
plt.show()

Result = []
Result1 = []
MaxH = 40/pow(10,6)
MinH = 38/pow(10,6)
Mezytilen.NachalnayaTemperarura = 300
k = MaxH
F2 = 4 / 3 * 3.14 * pow(Mezytilen.RadiusKapla, 3) * Mezytilen.PlotnostKapla * G
Dt = 48
Heigh = MaxH
F4 = 0
for i in range(1, round((MaxH - MinH)*pow(10,6)*1000)):
    KekLol = 3.14 * 24 * Azot.DinamicVyazkost * pow(Mezytilen.RadiusKapla, 4) * Azot.KoefTeploprovodnosti*Dt/(Azot.PlotnostPara*Azot.TeplotaParoobrazovaniya*pow(Heigh,4))*0.25-F2
    Result.append(KekLol*1000)
    Result1.append(Heigh*pow(10,6))
    Heigh = k - i/pow(10,6)/1000
plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
plt.xlabel('H, mkm')
plt.ylabel('F, mN')
plt.title('F ot dT')
plt.grid(True)
plt.show()



Result = []
Result1 = []
MaxRadius = 2.5/1000
MinRadius = 1/1000
Mezytilen.RadiusKapla = MaxRadius
k = Mezytilen.RadiusKapla
for i in range(1, round((MaxRadius - MinRadius)*1000000)):
    F2 = 4 / 3 * 3.14 * pow(Mezytilen.RadiusKapla, 3) * Mezytilen.PlotnostKapla * G
    KekLol = 3.14 * 24 * Azot.DinamicVyazkost * pow(Mezytilen.RadiusKapla, 4) * Azot.KoefTeploprovodnosti*200/(Azot.PlotnostPara*Azot.TeplotaParoobrazovaniya*pow(50/pow(10,6),4))*0.25
    F4 = KekLol - F2
    Result.append(F4*1000)
    Result1.append(Mezytilen.RadiusKapla*1000)
    Mezytilen.RadiusKapla = k - i/1000/1000
plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
plt.xlabel('R, mm')
plt.ylabel('F, mN')
plt.title('F ot R')
plt.grid(True)
plt.show()




def CalcHWave(Jija, Kapla):
    H = pow(Jija.DinamicVyazkost*Kapla.RadiusKapla*ConstantStef_Bolth*(pow(Kapla.NachalnayaTemperarura, 4)-pow(Jija.JijiTemperarura, 4))/(Jija.PlotnostJija*Kapla.PlotnostKapla*G*Jija.TeplotaParoobrazovaniya), 1/3)
    return H

def CalcHTipa(Jija, Kapla):
    H = pow(Jija.DinamicVyazkost * Jija.KoefTeploprovodnosti * Kapla.RadiusKapla * (Kapla.NachalnayaTemperarura - Jija.JijiTemperarura)/(Jija.PlotnostJija*Kapla.PlotnostKapla*G*Jija.TeplotaParoobrazovaniya), 0.25)
    return H

print(CalcHWave(Azot, Mezytilen)*1000000*5)
print(' mkm')
print()
print(CalcHTipa(Azot, Mezytilen)*1000000*5)
print(' mkm')
print()

def CalcHV2(Jija, Kapla):
    return 0

def ForceGravity(Kapla):
    ForceGravity = 4 / 3 * np.pi * pow(Kapla.RadiusKapla, 3) * Kapla.PlotnostKapla * G
    return ForceGravity

def ForcePressure(Jija, Kapla):
    ForcePressure = np.pi * 24 * Jija.DinamicVyazkost * pow(Kapla.RadiusKapla, 4) * Jija.KoefTeploprovodnosti * (Kapla.NachalnayaTemperatura - Jija.JijiTemperatura) / (
                Jija.PlotnostPara * Jija.TeplotaParoobrazovaniya * pow(CalcHV2(Jija, Kapla), 4)) * PressureIntegral
    return ForcePressure

class PhysicalParametr:
    def __init__(self, large_name, small_name, values):
        self.large_name = large_name
        self.small_name = small_name
        self.values = values

    def add_new_method(self, Name, func):
        setattr(self, Name, func)
    
    def SetValue(self, values):
        self.values = values
    
    def GetValue(self):
        return self.values
        
# Начинаем создавать физичекие параметры
InitialMassive = [] 
Force = PhysicalParametr("Force", "N", InitialMassive)
Force.add_new_method('Gravity', ForceGravity)
Force.add_new_method('Pressure', ForcePressure)
K = Force.Gravity(Mezytilen)
Force.SetValue(K)
print (Force.GetValue())
# class Power:
#     def __init__(self, Name, Result):
#         self.Name = Name
#         self.Result = Result

# def plotGraphic_TotR(arguments1, arguments2):
#     plt.plot(arguments1, arguments2, marker='o', color='g', linestyle='-')
#     plt.xlabel(arguments1.Name)
#     plt.ylabel(arguments2.Name)
#     plt.title(arguments1.Name + "ot" + arguments2.Name)
#     plt.grid(True)
#     plt.show()

Result = []
Result1 = []
MaxRadius = 2.5/1000
MinRadius = 1/1000
Mezytilen.RadiusKapla = MaxRadius
k = Mezytilen.RadiusKapla
for i in range(1, round((MaxRadius - MinRadius)*1000000)):
    F2 = 4 / 3 * 3.14 * pow(Mezytilen.RadiusKapla, 3) * Mezytilen.PlotnostKapla * G
    KekLol = 3.14 * 24 * Azot.DinamicVyazkost * pow(Mezytilen.RadiusKapla, 4) * Azot.KoefTeploprovodnosti*200/(Azot.PlotnostPara*Azot.TeplotaParoobrazovaniya*pow(10/pow(10,6),4))*0.25
    F4 = KekLol - F2
    Result.append(F4*1000)
    Result1.append(Mezytilen.RadiusKapla*1000)
    Mezytilen.RadiusKapla = k - i/1000/1000
plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
plt.xlabel('T0, K')
plt.ylabel('t, s')
plt.title('t ot T0')
plt.grid(True)
plt.show()



K = Force.Gravity(Mezytilen)
G = Time.Ohlagdeeniya(Mezytilen, Azot)
Force.SetValue(K)
Time.SetValue(G)
Plotter(Force, Time)



T0Array = []
TimeArray = []
T0 = 300
KekW = T0
TL = 127
SpeedTemp = 1500
SecInMin = 60
for i in range (1, T0 - TL):
    Gopa = (T0 - TL)/(1500/SecInMin)
    T0Array.append(T0)
    TimeArray.append(Gopa)
    T0 = KekW - i
    plt.plotKav(Gavno1, Gavno2)
plt.plot(T0Array, TimeArray, marker='o', color='g', linestyle='-')
plt.xlabel('T0, K')
plt.ylabel('t, s')
plt.title('t ot T0')
plt.grid(True)
plt.show()

ExperimentalTime = [36, 31, 22, 20, 20, 26, 25, 29, 25, 23, 18, 25, 27, 17, 20, 23, 29, 28, 31, 31, 17, 23, 33, 23, 19, 23, 25]
Kek228 = []
Padla = 0
for i in range(1, len(ExperimentalTime)):
   Padla = Padla + ExperimentalTime[i-1]
   Kek228.append(Padla/i)
plt.plot(np.linspace(0,len(Kek228),len(Kek228)), Kek228, marker='o', color='g', linestyle='-')
plt.xlabel('R, mm')
plt.ylabel('t, s')
plt.title('t ot R')
plt.grid(True)
plt.show()


ExperimentalTime = [19, 17, 29, 12, 14, 14, 26, 20, 13, 15, 28, 20, 12, 20, 23, 24, 18, 17, 19, 25, 23.5, 25, 17, 14, 27, 19, 13]
print(np.average(ExperimentalTime))
Kek228 = []
Padla = 0
for i in range(1, len(ExperimentalTime)):
   Padla = Padla + ExperimentalTime[i-1]
   Kek228.append(Padla/i)
plt.plot(np.linspace(0,len(Kek228),len(Kek228)), Kek228, marker='o', color='g', linestyle='-')
plt.xlabel('R, mm')
plt.ylabel('t, s')
plt.title('t ot R')
plt.grid(True)
plt.show()

Xot = 3.6/1000/2
Xdo = 4.2/1000/2
CountsOfX = pow(10,8)
mu = (Xdo + Xot)/2
sigma = (Xdo - mu)/2
RadiusArray = np.random.normal(mu, sigma, CountsOfX)
Mezytilen.RadiusKapla = RadiusArray
ConstantDrop = CalcConstantDrop(Azot, Mezytilen)
t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen)
t2 = CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen)
t3 = CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen)
t = t1 + t2 + t3
Mezytilen.RadiusKapla = mu
CalcConstantDrop(Azot, Mezytilen)
tKek = CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen) + CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen) + CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen)
Average = np.mean(t)
Kek228 = []
Kek228.append(t[0])
Padla = 0
for i in range(1, CountsOfX+1):
   Padla = Padla + t[i-1]
   Kek228.append(Padla/i)

print(Average)
print(Kek228[CountsOfX])
Kek = Average - Kek228[CountsOfX]
print (Kek)

plt.plot(np.linspace(0,CountsOfX,len(Kek228)), Kek228, marker='o', color='g', linestyle='-')
plt.xlabel('R, mm')
plt.ylabel('t, s')
plt.title('t ot R')
plt.grid(True)
plt.show()

PressureIntegral = 0.25



Result = []
Result1 = []
MaxRadius = 2.5/1000
MinRadius = 1/1000
Mezytilen.RadiusKapla = MaxRadius
k = Mezytilen.RadiusKapla
for i in range(1, round((MaxRadius - MinRadius)*1000000)):
    F2 = 4 / 3 * 3.14 * pow(Mezytilen.RadiusKapla, 3) * Mezytilen.PlotnostKapla * G
    KekLol = 3.14 * 24 * Azot.DinamicVyazkost * pow(Mezytilen.RadiusKapla, 4) * Azot.KoefTeploprovodnosti*200/(Azot.PlotnostPara*Azot.TeplotaParoobrazovaniya*pow(50/pow(10,6),4))*PressureIntegral
    F4 = KekLol - F2
    Result.append(F4*1000)
    Result1.append(Mezytilen.RadiusKapla*1000)
    Mezytilen.RadiusKapla = k - i/1000/1000
k = Power("Н", Result1)
g = Power("s", Result)
plotGraphic_TotR(Result1, Result)

Result = []
Result1 = []
Mezytilen.RadiusKapla = MaxRadius
DT = 200
k = DT
for i in range(1, 200):
   F2 = 4 / 3 * 3.14 * pow(Mezytilen.RadiusKapla, 3) * Mezytilen.PlotnostKapla * G
   KekLol = 3.14 * 24 * Azot.DinamicVyazkost * pow(Mezytilen.RadiusKapla, 4) * Azot.KoefTeploprovodnosti*DT/(Azot.PlotnostPara*Azot.TeplotaParoobrazovaniya*pow(50/pow(10,6),4))*PressureIntegral
   F4 = KekLol - F2
   Result.append(F4*1000)
   Result1.append(DT)
   DT = k - i
plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
plt.xlabel('DT, K')
plt.ylabel('F, mN')
plt.title('F ot DT')
plt.grid(True)
plt.show()

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
    P = 24 * 27 * Azot.DinamicVyazkost * CalcConstantDrop(Azot, Mezytilen)/Mezytilen.RadiusKapla*0.125
    S = 2 * 3.14 * pow(Mezytilen.RadiusKapla, 2)
    F = P * S
    F2 = 4 / 3 * 3.14 * pow(Mezytilen.RadiusKapla, 3) * Mezytilen.PlotnostKapla * G
    F3 = F - F2
    Mezytilen.RadiusKapla = Mezytilen.RadiusKapla - i/1000000

 #CalcNormalGause(mu, sigma, Color, ColorArea, NameOfPlot, NameOfY, NameOfX, NameOfArea, Xot, Xdo, CountsOfX, linewidth)
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


MaxRadius = 2.5/1000
MinRadius = 1/1000
Mezytilen.RadiusKapla = MaxRadius
Result = []
Result1 = []
k = Mezytilen.NachalnayaTemperarura
for i in range(1, 100*round(Mezytilen.NachalnayaTemperarura - Mezytilen.PlavleniyaTemperarura-1)):
    RMez = CalcMaxRadius(Azot, Mezytilen)
    ConstantDrop = CalcConstantDrop(Azot, Mezytilen)
    t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen)
    t2 = CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen)
    t3 = CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen)
    t = t1 + t2 + t3
    Result.append(t1)
    Result1.append(Mezytilen.NachalnayaTemperarura)
    Mezytilen.NachalnayaTemperarura = k - i/100
plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
plt.xlabel('T0, K')
plt.ylabel('t1, s')
plt.title('t1 ot T0')
plt.grid(True)
plt.show()

Result = []
Result1 = []
MaxRadius = 2.5/1000
MinRadius = 1/1000
Mezytilen.RadiusKapla = MaxRadius
k = Mezytilen.RadiusKapla
for i in range(1, round((MaxRadius - MinRadius)*1000000)):
    ConstantDrop = CalcConstantDrop(Azot, Mezytilen)
    t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen)
    t2 = CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen)
    t3 = CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen)
    t = t1 + t2 + t3
    Result.append(t1)
    Result1.append(Mezytilen.RadiusKapla*1000)
    Mezytilen.RadiusKapla = k - i/1000/1000
plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
plt.xlabel('R, mm')
plt.ylabel('t1, s')
plt.title('t1 ot R')
plt.grid(True)
plt.show()

# Result = []
# Result1 = []
# Azot.LeidenfrostaTemperarura = 149
# k = Azot.LeidenfrostaTemperarura
# for i in range(1, 7000):
#     RMet = CalcMaxRadius(Azot, Mezytilen)
#     Mezytilen.RadiusKapla = RMet*2
#     ConstantDrop = CalcConstantDrop(Azot, Mezytilen)
#     t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Azot, Mezytilen)
#     t2 = CalcDlitelnostKristall(ConstantDrop, Azot, Mezytilen)
#     t3 = CalcDlitelnostOstivaniya(ConstantDrop, Azot, Mezytilen)
#     t = t1 + t2 + t3
#     Result.append(t)
#     Result1.append(Azot.LeidenfrostaTemperarura)
#     Azot.LeidenfrostaTemperarura = k - i/100
# plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
# plt.xlabel('TL, K')
# plt.ylabel('t, s')
# plt.title('t ot TL, Mezytilen in Azot')
# plt.grid(True)
# plt.show()


# Result = []
# Result1 = []
# k = Neon.LeidenfrostaTemperarura
# for i in range(1, 7000):
#     RMet = CalcMaxRadius(Neon, Metan)
#     Metan.RadiusKapla = RMet*2
#     ConstantDrop = CalcConstantDrop(Neon, Metan)
#     t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Neon, Metan)
#     t2 = CalcDlitelnostKristall(ConstantDrop, Neon, Metan)
#     t3 = CalcDlitelnostOstivaniya(ConstantDrop, Neon, Metan)
#     t = t1 + t2 + t3
#     Result.append(t)
#     Result1.append(Neon.LeidenfrostaTemperarura)
#     Neon.LeidenfrostaTemperarura = k - i/100
# plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
# plt.xlabel('TL, K')
# plt.ylabel('t, s')
# plt.title('t ot TL, Metan in Neon')
# plt.grid(True)
# plt.show()

# Result = []
# Result1 = []
# k = Hidrogen.LeidenfrostaTemperarura
# for i in range(1, 8000):
#     RMet = CalcMaxRadius(Hidrogen, Metan)
#     Metan.RadiusKapla = RMet*2
#     ConstantDrop = CalcConstantDrop(Hidrogen, Metan)
#     t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, Hidrogen, Metan)
#     t2 = CalcDlitelnostKristall(ConstantDrop, Hidrogen, Metan)
#     t3 = CalcDlitelnostOstivaniya(ConstantDrop, Hidrogen, Metan)
#     t = t1 + t2 + t3
#     Result.append(t)
#     Result1.append(Hidrogen.LeidenfrostaTemperarura)
#     Hidrogen.LeidenfrostaTemperarura = k - i/100
# plt.plot(Result1, Result, marker='o', color='g', linestyle='-')
# plt.xlabel('TL, K')
# plt.ylabel('t, s')
# plt.title('t ot TL, Metan in Hydrogen')
# plt.grid(True)
# plt.show()


# k = 3
# for g in range(1, k):
#     PropertiesSharik.NachalnayaTemperarura = 230
#     for i in range(1, 101):
    
#         ConstantDrop = CalcConstantDrop(PropertiesSharik.PlotnostKapla, PropertiesSharik.TeploEmkostKapla, PropertiesSharik.RadiusKapla, PropertiesSharik.DolaPoverhnostiKapli, PropertiesJija.DinamicVyazkost, PropertiesJija.PlotnostPara, PropertiesJija.KoefTeploprovodnosti, PropertiesJija.TeplotaParoobrazovaniya)
#         t1 = CalcDlitelnostOhlagdeniya(ConstantDrop, PropertiesSharik.NachalnayaTemperarura, PropertiesSharik.PlavleniyaTemperarura, PropertiesJija.JijiTemperarura)
#         t2 = CalcDlitelnostKristall(ConstantDrop, PropertiesSharik.TeploEmkostKapla, PropertiesSharik.TeplotaPlavleniya, PropertiesSharik.PlavleniyaTemperarura, PropertiesJija.JijiTemperarura)
#         t3 = CalcDlitelnostOstivaniya(ConstantDrop, PropertiesSharik.TeploEmkostKapla, PropertiesSharik.TeploEmkostKristall, PropertiesJija.LeidenfrostaTemperarura, PropertiesSharik.PlavleniyaTemperarura, PropertiesJija.JijiTemperarura)

#         t = t1 + t2 + t3
        
#         Result.append(t)
    
#         PropertiesSharik.NachalnayaTemperarura = 230 + i/100
    
    
      
#     PropertiesSharik.RadiusKapla =  (3.5+g/1000) / 2 / 1000
   
#     # if Result[i-1] < 0:
#     #     print(i)
#     Result1 = []
#     Result2 = []
# for i in range (1, 101):
#     Result1.append(Result[i-1])
#     g = i + 99
#     Result2.append(Result[g])

# NachalnayaTemperarura = np.linspace(230, 300, num=len(Result1))
# # plot_dependency(Result1, NachalnayaTemperarura) 
# plt.plot(NachalnayaTemperarura, Result1, marker='o', color='b', linestyle='-')
# plt.plot(NachalnayaTemperarura, Result2, marker='o', color='r', linestyle='-')
# plt.xlabel('T0')
# plt.ylabel('t')
# plt.title('t ot T0')
# plt.grid(True)
# plt.show()
   

# Пример данных для построения графика

