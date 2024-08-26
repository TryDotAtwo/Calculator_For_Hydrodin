import matplotlib.pyplot as plt


def plotGraphic(arguments1, arguments2):
    plt.plot(arguments1.values, arguments2.values, marker='o', color='g', linestyle='-')
    plt.xlabel(arguments1.small_name + ", " + arguments1.units)
    plt.ylabel(arguments2.small_name + ", " + arguments2.units)
    plt.title(arguments1.large_name + ' ot ' + arguments2.large_name)
    plt.grid(True)
    plt.show()