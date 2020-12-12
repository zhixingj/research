import pandas as pd
import matplotlib.pyplot as plt
import math as math
import numpy as np
from astropy.modeling import models
from matplotlib import pyplot as plt, pylab
from matplotlib.pyplot import *
from scipy.interpolate import interpolate
from numpy import ma
from matplotlib import cbook
from matplotlib.colors import Normalize

import matplotlib.colors as colors
from numpy import zeros

np.set_printoptions(threshold=np.inf)

pd.set_option("display.precision", 4)

Hill_2019 = pd.read_csv('Hill_2019_2.txt', sep="\t",
                        names=['[Fe/H]', '[O/Fe]', '[Na/Fe]', '[Mg/Fe]', '[Si/Fe]', '[Ca/Fe]', '[Sc/Fe]', '[Ti1/Fe]',
                               '[Ti2/Fe]', '[Cr/Fe]', '[Fe2/Fe]', '[Co/Fe]', '[Ni/Fe]', '[Zn/Fe]', '[Ba/Fe]', '[La/Fe]',
                               '[Nd/Fe]', '[Eu/Fe]', '[Sr/Fe]', '[Al1/Fe]', '[Si1/Fe]', '[Sc2/Fe]', '[Mn/Fe]',
                               '[Sr2/Fe]', '[Y2/Fe]', '[Ba2/Fe]', '[Eu2/Fe]', '[C/Fe]'])
Hill_2019_Uncertainty = pd.read_csv('Hill_2019_Uncertainty.txt', sep="\t",
                        names=['e[Fe/H]', 'e[O/Fe]', 'e[Na/Fe]', 'e[Mg/Fe]', 'e[Si/Fe]', 'e[Ca/Fe]', 'e[Sc/Fe]', 'e[Ti1/Fe]',
                               'e[Ti2/Fe]', 'e[Cr/Fe]', 'e[Fe2/Fe]', 'e[Co/Fe]', 'e[Ni/Fe]', 'e[Zn/Fe]', 'e[Ba/Fe]', 'e[La/Fe]',
                               'e[Nd/Fe]', 'e[Eu/Fe]', 'e[Sr/Fe]', 'e[Al1/Fe]', 'e[Si1/Fe]', 'e[Sc2/Fe]', 'e[Mn/Fe]',
                               'e[Sr2/Fe]', 'e[Y2/Fe]', 'e[Ba2/Fe]', 'e[Eu2/Fe]', 'e[C/Fe]'])
Starkenburg_2013 = pd.read_csv('Starkenburg_2013_2.txt', sep="\t",
                               names=['[Fe/H]', '[O/Fe]', '[Na/Fe]', '[Mg/Fe]', '[Si/Fe]', '[Ca/Fe]', '[Sc/Fe]',
                                      '[Ti1/Fe]', '[Ti2/Fe]', '[Cr/Fe]', '[Fe2/Fe]', '[Co/Fe]', '[Ni/Fe]', '[Zn/Fe]',
                                      '[Ba/Fe]', '[La/Fe]', '[Nd/Fe]', '[Eu/Fe]', '[Sr/Fe]', '[Al1/Fe]', '[Si1/Fe]',
                                      '[Sc2/Fe]', '[Mn/Fe]', '[Sr2/Fe]', '[Y2/Fe]', '[Ba2/Fe]', '[Eu2/Fe]', '[C/Fe]'])
Starkenburg_2013_Uncertainty = pd.read_csv('Starkenburg_2013_Uncertainty.txt', sep="\t",
                        names=['e[Fe/H]', 'e[O/Fe]', 'e[Na/Fe]', 'e[Mg/Fe]', 'e[Si/Fe]', 'e[Ca/Fe]', 'e[Sc/Fe]', 'e[Ti1/Fe]',
                               'e[Ti2/Fe]', 'e[Cr/Fe]', 'e[Fe2/Fe]', 'e[Co/Fe]', 'e[Ni/Fe]', 'e[Zn/Fe]', 'e[Ba/Fe]', 'e[La/Fe]',
                               'e[Nd/Fe]', 'e[Eu/Fe]', 'e[Sr/Fe]', 'e[Al1/Fe]', 'e[Si1/Fe]', 'e[Sc2/Fe]', 'e[Mn/Fe]',
                               'e[Sr2/Fe]', 'e[Y2/Fe]', 'e[Ba2/Fe]', 'e[Eu2/Fe]', 'e[C/Fe]'])
Tafelmeyer_2018 = pd.read_csv('Tafelmeyer_2018_2.txt', sep="\t",
                              names=['[Fe/H]', '[O/Fe]', '[Na/Fe]', '[Mg/Fe]', '[Si/Fe]', '[Ca/Fe]', '[Sc/Fe]',
                                     '[Ti1/Fe]', '[Ti2/Fe]', '[Cr/Fe]', '[Fe2/Fe]', '[Co/Fe]', '[Ni/Fe]', '[Zn/Fe]',
                                     '[Ba/Fe]', '[La/Fe]', '[Nd/Fe]', '[Eu/Fe]', '[Sr/Fe]', '[Al1/Fe]', '[Si1/Fe]',
                                     '[Sc2/Fe]', '[Mn/Fe]', '[Sr2/Fe]', '[Y2/Fe]', '[Ba2/Fe]', '[Eu2/Fe]', '[C/Fe]'])
Tafelmeyer_2018_Uncertainty = pd.read_csv('Tafelmeyer_2018_Uncertainty.txt', sep="\t",
                        names=['e[Fe/H]', 'e[O/Fe]', 'e[Na/Fe]', 'e[Mg/Fe]', 'e[Si/Fe]', 'e[Ca/Fe]', 'e[Sc/Fe]', 'e[Ti1/Fe]',
                               'e[Ti2/Fe]', 'e[Cr/Fe]', 'e[Fe2/Fe]', 'e[Co/Fe]', 'e[Ni/Fe]', 'e[Zn/Fe]', 'e[Ba/Fe]', 'e[La/Fe]',
                               'e[Nd/Fe]', 'e[Eu/Fe]', 'e[Sr/Fe]', 'e[Al1/Fe]', 'e[Si1/Fe]', 'e[Sc2/Fe]', 'e[Mn/Fe]',
                               'e[Sr2/Fe]', 'e[Y2/Fe]', 'e[Ba2/Fe]', 'e[Eu2/Fe]', 'e[C/Fe]'])
Kirby_2009 = pd.read_csv('Kirby_2009_2.txt', sep="\t",
                         names=['[Fe/H]', '[O/Fe]', '[Na/Fe]', '[Mg/Fe]', '[Si/Fe]', '[Ca/Fe]', '[Sc/Fe]', '[Ti1/Fe]',
                                '[Ti2/Fe]', '[Cr/Fe]', '[Fe2/Fe]', '[Co/Fe]', '[Ni/Fe]', '[Zn/Fe]', '[Ba/Fe]',
                                '[La/Fe]', '[Nd/Fe]', '[Eu/Fe]', '[Sr/Fe]', '[Al1/Fe]', '[Si1/Fe]', '[Sc2/Fe]',
                                '[Mn/Fe]', '[Sr2/Fe]', '[Y2/Fe]', '[Ba2/Fe]', '[Eu2/Fe]', '[C/Fe]'])
Kirby_2009_Uncertainty = pd.read_csv('Kirby_2009_Uncertainty.txt', sep="\t",
                        names=['e[Fe/H]', 'e[O/Fe]', 'e[Na/Fe]', 'e[Mg/Fe]', 'e[Si/Fe]', 'e[Ca/Fe]', 'e[Sc/Fe]', 'e[Ti1/Fe]',
                               'e[Ti2/Fe]', 'e[Cr/Fe]', 'e[Fe2/Fe]', 'e[Co/Fe]', 'e[Ni/Fe]', 'e[Zn/Fe]', 'e[Ba/Fe]', 'e[La/Fe]',
                               'e[Nd/Fe]', 'e[Eu/Fe]', 'e[Sr/Fe]', 'e[Al1/Fe]', 'e[Si1/Fe]', 'e[Sc2/Fe]', 'e[Mn/Fe]',
                               'e[Sr2/Fe]', 'e[Y2/Fe]', 'e[Ba2/Fe]', 'e[Eu2/Fe]', 'e[C/Fe]'])
Kirby_2018 = pd.read_csv('Kirby_2018.txt', sep="\t",
                         names=['[Fe/H]', '[O/Fe]', '[Na/Fe]', '[Mg/Fe]', '[Si/Fe]', '[Ca/Fe]', '[Sc/Fe]', '[Ti1/Fe]',
                                '[Ti2/Fe]', '[Cr/Fe]', '[Fe2/Fe]', '[Co/Fe]', '[Ni/Fe]', '[Zn/Fe]', '[Ba/Fe]',
                                '[La/Fe]', '[Nd/Fe]', '[Eu/Fe]', '[Sr/Fe]', '[Al1/Fe]', '[Si1/Fe]', '[Sc2/Fe]',
                                '[Mn/Fe]', '[Sr2/Fe]', '[Y2/Fe]', '[Ba2/Fe]', '[Eu2/Fe]', '[C/Fe]'])
Kirby_2018_Uncertainty = pd.read_csv('Kirby_2018_Uncertainty.txt', sep="\t",
                        names=['e[Fe/H]', 'e[O/Fe]', 'e[Na/Fe]', 'e[Mg/Fe]', 'e[Si/Fe]', 'e[Ca/Fe]', 'e[Sc/Fe]', 'e[Ti1/Fe]',
                               'e[Ti2/Fe]', 'e[Cr/Fe]', 'e[Fe2/Fe]', 'e[Co/Fe]', 'e[Ni/Fe]', 'e[Zn/Fe]', 'e[Ba/Fe]', 'e[La/Fe]',
                               'e[Nd/Fe]', 'e[Eu/Fe]', 'e[Sr/Fe]', 'e[Al1/Fe]', 'e[Si1/Fe]', 'e[Sc2/Fe]', 'e[Mn/Fe]',
                               'e[Sr2/Fe]', 'e[Y2/Fe]', 'e[Ba2/Fe]', 'e[Eu2/Fe]', 'e[C/Fe]'])
Duggan_2018 = pd.read_csv('Duggan_2018.txt', sep="\t",
                         names=['[Fe/H]', '[O/Fe]', '[Na/Fe]', '[Mg/Fe]', '[Si/Fe]', '[Ca/Fe]', '[Sc/Fe]', '[Ti1/Fe]',
                                '[Ti2/Fe]', '[Cr/Fe]', '[Fe2/Fe]', '[Co/Fe]', '[Ni/Fe]', '[Zn/Fe]', '[Ba/Fe]',
                                '[La/Fe]', '[Nd/Fe]', '[Eu/Fe]', '[Sr/Fe]', '[Al1/Fe]', '[Si1/Fe]', '[Sc2/Fe]',
                                '[Mn/Fe]', '[Sr2/Fe]', '[Y2/Fe]', '[Ba2/Fe]', '[Eu2/Fe]', '[C/Fe]'])
Duggan_2018_Uncertainty = pd.read_csv('Duggan_2018_Uncertainty.txt', sep="\t",
                        names=['e[Fe/H]', 'e[O/Fe]', 'e[Na/Fe]', 'e[Mg/Fe]', 'e[Si/Fe]', 'e[Ca/Fe]', 'e[Sc/Fe]', 'e[Ti1/Fe]',
                               'e[Ti2/Fe]', 'e[Cr/Fe]', 'e[Fe2/Fe]', 'e[Co/Fe]', 'e[Ni/Fe]', 'e[Zn/Fe]', 'e[Ba/Fe]', 'e[La/Fe]',
                               'e[Nd/Fe]', 'e[Eu/Fe]', 'e[Sr/Fe]', 'e[Al1/Fe]', 'e[Si1/Fe]', 'e[Sc2/Fe]', 'e[Mn/Fe]',
                               'e[Sr2/Fe]', 'e[Y2/Fe]', 'e[Ba2/Fe]', 'e[Eu2/Fe]', 'e[C/Fe]'])
Skuladottir_2017 = pd.read_csv('Skuladottir_2017.txt', sep="\t",
                         names=['[Fe/H]', '[O/Fe]', '[Na/Fe]', '[Mg/Fe]', '[Si/Fe]', '[Ca/Fe]', '[Sc/Fe]', '[Ti1/Fe]',
                                '[Ti2/Fe]', '[Cr/Fe]', '[Fe2/Fe]', '[Co/Fe]', '[Ni/Fe]', '[Zn/Fe]', '[Ba/Fe]',
                                '[La/Fe]', '[Nd/Fe]', '[Eu/Fe]', '[Sr/Fe]', '[Al1/Fe]', '[Si1/Fe]', '[Sc2/Fe]',
                                '[Mn/Fe]', '[Sr2/Fe]', '[Y2/Fe]', '[Ba2/Fe]', '[Eu2/Fe]', '[C/Fe]'])
Skuladottir_2017_Uncertainty = pd.read_csv('Skuladottir_2017_Uncertainty.txt', sep="\t",
                        names=['e[Fe/H]', 'e[O/Fe]', 'e[Na/Fe]', 'e[Mg/Fe]', 'e[Si/Fe]', 'e[Ca/Fe]', 'e[Sc/Fe]', 'e[Ti1/Fe]',
                               'e[Ti2/Fe]', 'e[Cr/Fe]', 'e[Fe2/Fe]', 'e[Co/Fe]', 'e[Ni/Fe]', 'e[Zn/Fe]', 'e[Ba/Fe]', 'e[La/Fe]',
                               'e[Nd/Fe]', 'e[Eu/Fe]', 'e[Sr/Fe]', 'e[Al1/Fe]', 'e[Si1/Fe]', 'e[Sc2/Fe]', 'e[Mn/Fe]',
                               'e[Sr2/Fe]', 'e[Y2/Fe]', 'e[Ba2/Fe]', 'e[Eu2/Fe]', 'e[C/Fe]'])

Lodders_2020 = pd.read_csv('Lodders_2020.txt', sep="\t", names=['Abundance'])
AG_1989 = pd.read_csv('A&G_1989.txt', sep="\t", names=['Abundance'])


# method takes textfile name and datafile
# method iterates over every column in the datafile
# for every value above, method goes to corresponding isotope solar abundance in Lodders text file and given textfile
# method recalculates the value above by adding log of solar abundance from given textfile and subtracting solar abundance from Lodders text file
def callibrateSolarAbundances(df, textfile):

    df.loc["[Fe/H]"] += textfile.at["Fe56", "Abundance"] - textfile.at["H1", "Abundance"] - Lodders_2020.at[
        "Fe56", "Abundance"] + Lodders_2020.at["H1", "Abundance"]
    df.loc["[O/Fe]"] += textfile.at["O16", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "O16", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Na/Fe]"] += textfile.at["Na23", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Na23", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Mg/Fe]"] += textfile.at["Mg24", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Mg24", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Si/Fe]"] += textfile.at["Si28", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Si28", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Ca/Fe]"] += textfile.at["Ca40", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Ca40", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Sc/Fe]"] += textfile.at["Sc45", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Sc45", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Ti1/Fe]"] += textfile.at["Ti48", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Ti48", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Ti2/Fe]"] += textfile.at["Ti46", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Ti46", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Cr/Fe]"] += textfile.at["Cr52", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Cr52", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Fe2/Fe]"] += textfile.at["Fe54", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Fe54", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Co/Fe]"] += textfile.at["Co59", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Co59", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Ni/Fe]"] += textfile.at["Ni58", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Ni58", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Zn/Fe]"] += textfile.at["Zn64", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Zn64", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Ba/Fe]"] += textfile.at["Ba138", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Ba138", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[La/Fe]"] += textfile.at["La139", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "La139", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Nd/Fe]"] += textfile.at["Nd142", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Nd142", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Eu/Fe]"] += textfile.at["Eu153", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Eu153", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Sr/Fe]"] += textfile.at["Sr86", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Sr86", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Al1/Fe]"] += textfile.at["Al27", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Al27", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Si1/Fe]"] += textfile.at["Si28", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Si28", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]

    # no values for Sc2 in A&G

    df.loc["[Mn/Fe]"] += textfile.at["Mn55", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Mn55", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Sr2/Fe]"] += textfile.at["Sr87", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Sr87", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Y2/Fe]"] += textfile.at["Y89", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Y89", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Ba2/Fe]"] += textfile.at["Ba137", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Ba137", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[Eu2/Fe]"] += textfile.at["Eu151", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "Eu151", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]
    df.loc["[C/Fe]"] += textfile.at["C12", "Abundance"] - textfile.at["Fe56", "Abundance"] - Lodders_2020.at[
        "C12", "Abundance"] + Lodders_2020.at["Fe56", "Abundance"]

def getMaxBound(element, Data):
    for i in Data.columns:
        if i == element:
            element = Data.columns.get_loc(i)
            return Data.iloc[:, element].max()

def getMinBound(element, Data):
    for i in Data.columns:
        if i == element:
            element = Data.columns.get_loc(i)
            return Data.iloc[:, element].min()

def getMaxBoundFe(element, Data):
    max = -100
    FeH = Data.columns.get_loc("[Fe/H]")
    for i in Data.columns:
        if i == element:
            elementPos = Data.columns.get_loc(i)
            for j in range(0, len(Data.index)):
                if pd.notna(Data.iloc[j, elementPos]):
                    x0 = Data.iloc[j, FeH]
                    if x0 > max:
                        max = x0
    return max

def getMinBoundFe(element, Data):
    min = 100
    FeH = Data.columns.get_loc("[Fe/H]")
    for i in Data.columns:
        if i == element:
            elementPos = Data.columns.get_loc(i)
            for j in range(0, len(Data.index)):
                if pd.notna(Data.iloc[j, elementPos]):
                    x0 = Data.iloc[j, FeH]
                    if x0 < min:
                        min = x0
    return min


# Creates plot for all isotopic abundance vs metallicity graphs, and can also take an isotope as parameter and plot Isotopic Abundance vs Metallicity for that element
def makePlot(df, element=None):
    if element is None:
        for i in df.columns:
            if df.columns.get_loc(i) > 0:
                if i != '[Fe/H]':
                    df.plot.scatter(x='[Fe/H]', y=df.columns.get_loc(i), color='red', marker= ".")
                    plt.show(block=True)
    else:
        for i in df.columns:
            if i == element:
                df.plot.scatter(x='[Fe/H]', y=df.columns.get_loc(i), color='red', marker = ".")
                plt.show(block=True)


# Creates plot for all isotopic abundance vs metallicity graphs, and can also take an isotope as parameter and plot Isotopic Abundance vs Metallicity for that element. Also shows source of data, and seperates data into colors
def createPlotWithSources(df, element=None):
    if element is None:
        for i in df.columns:
            print i
            if df.columns.get_loc(i) > 0:
                if i != '[Fe/H]':
                    hill_plot = Hill_2019.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i), color='blue',
                                                       marker='.',
                                                       label="Hill_2019", extent=(getMinBound("[Fe/H]", df), getMaxBound("[Fe/H]", df),
                                                       getMinBound(element, df), getMaxBound(element, df)))
                    kirby_plot = Kirby_2009.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i), color='green',
                                                         marker="*",
                                                         label="Kirby_2009", ax=hill_plot)
                    starkenburg_plot = Starkenburg_2013.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                                     color='chocolate', marker="+",
                                                                     label="Starkenburg_2013",
                                                                     ax=hill_plot)
                    tafelmeyer_plot = Tafelmeyer_2018.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                                   color='red', marker="x", label="Tafelmeyer_2018",
                                                                   ax=hill_plot)
                    kirby2_plot = Kirby_2018.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                          color='purple', marker="1", label="Kirby_2018",
                                                          ax=hill_plot)
                    Duggan_plot = Duggan_2018.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                          color='orange', marker="2", label="Duggan_2018",
                                                          ax=hill_plot)
                    Skuladottir_plot = Skuladottir_2017.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                           color='darkviolet', marker="3", label=r'Sk$\mathrm{\acute{u}}$lad$\mathrm{\acute{o}}$ttir',
                                                           ax=hill_plot)
                    plt.xlim(getMinBound("[Fe/H]", df) - 1, getMaxBound("[Fe/H]", df) + 1)
                    plt.ylim(getMinBound(element, df) - 1, getMaxBound(element, df) + 1)
                    return hill_plot
    else:
        for i in df.columns:
            if i == element:
                hill_plot = Hill_2019.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i), color='blue',
                                                   marker='.',
                                                   label="Hill_2019")
                kirby_plot = Kirby_2009.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i), color='green',
                                                     marker="*",
                                                     label="Kirby_2009", ax=hill_plot)
                starkenburg_plot = Starkenburg_2013.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                                 color='chocolate', marker="+",
                                                                 label="Starkenburg_2013",
                                                                 ax=hill_plot)
                tafelmeyer_plot = Tafelmeyer_2018.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                               color='red', marker="x", label="Tafelmeyer_2018",
                                                               ax=hill_plot)
                kirby2_plot = Kirby_2018.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                      color='purple', marker="1", label="Kirby_2018",
                                                      ax=hill_plot)
                Duggan_plot = Duggan_2018.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                       color='orange', marker="2", label="Duggan_2018",
                                                       ax=hill_plot)
                Skuladottir_plot = Skuladottir_2017.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                                 color='darkviolet', marker="3", label= r'Sk$\mathrm{\acute{u}}$lad$\mathrm{\acute{o}}$ttir',
                                                                 ax=hill_plot)
                plt.xlim(getMinBound("[Fe/H]", df) - 1, getMaxBound("[Fe/H]", df) + 1)
                plt.ylim(getMinBound(element, df) - 1, getMaxBound(element, df) + 1)
                return hill_plot


def showData(df, element=None):
    if element is None:
        print (df.to_string());
    else:
        for i in df.columns:
            if i == element:
                print(df.loc[:, element])

# Pass an array and a given value, and return the largest element in array less than value
def findCorrectBin(array, value):
    if value < array[0]:
        return 0

    elif value > array[len(array) - 1]:
        return (len(array) - 1)

    else:
        minValue = array[0]

        for i in array:
            if minValue < value:
                break;
            minValue = i

        for j in array:
            if j >= value:
                continue;
            if (j > minValue):
                minValue = j

        position = 0
        for k in range(0, len(array)):
            if array[k] == minValue:
                position = k

        return position

def calculateAverages(df):
    for column in df.columns:
        if column == "Mean":
            df.loc[:, "Mean"] = df.loc[:, "Sum"] / df.loc[:, "Number"]
        if column == "ErrorMean":
            df.loc[:, "ErrorMean"] = df.loc[:, "ErrorSum"] / df.loc[:, "Number"]
    return df

'''
# creates a 2d array by applying a gaussian spread function using metallicity abundance values and corresponding error values
def createGaussianArray(element, AbundanceData, ErrorData, z_max_FeH, z_max_Mg, logz_FeH, logz_Mg):
    rows = z_max_FeH
    columns = z_max_Mg
    FeH = AbundanceData.columns.get_loc("[Fe/H]")

    IsotopicAbundances = [[0 for x in range(columns)] for x in range(rows)]

    for i in AbundanceData.columns:
        if i == element:
            element = AbundanceData.columns.get_loc(i)
            for j in range(0, len(AbundanceData.index)):
                if pd.notna(AbundanceData.iloc[j, element]):
                    x0 = AbundanceData.iloc[j, FeH]
                    y0 = AbundanceData.iloc[j, element]

                    sigx = ErrorData.iloc[j, FeH]
                    sigy = ErrorData.iloc[j, element]

                    correctBinFeH = findCorrectBin(logz_FeH, AbundanceData.iloc[j, FeH])
                    correctBinElement = findCorrectBin(logz_Mg, AbundanceData.iloc[j, element])

                    for x in range(rows):
                        for y in range(columns):
                            newx = logz_FeH[x]
                            newy = logz_Mg[y]

                            fx = -0.5 * (((newx - x0) / sigx) ** 2)
                            fy = -0.5 * (((newy - y0) / sigy) ** 2)
                            fxresult = math.exp(fx)
                            fyresult = math.exp(fy)
                            ftot = fxresult * fyresult

                            IsotopicAbundances[correctBinFeH][correctBinElement] += ftot

    return IsotopicAbundances
'''

# creates a 2d plot by applying a gaussian spread function using metallicity abundance values and corresponding error values
def createGaussianPlot(element, AbundanceData, ErrorData, columns, rows, logz_FeH, logz_Element):
    FeH = AbundanceData.columns.get_loc("[Fe/H]")

    z = np.zeros((rows, columns))
    x, y = np.meshgrid(np.arange(0, rows), np.arange(0, columns), indexing="xy")

    for i in AbundanceData.columns:
        if i == element:
            element = AbundanceData.columns.get_loc(i)
            for j in range(0, len(AbundanceData.index)): #len(AbundanceData.index)
                if pd.notna(AbundanceData.iloc[j, element]):
                    x0 = AbundanceData.iloc[j, FeH]
                    y0 = AbundanceData.iloc[j, element]
                    #print str(x0) + ", " + str(y0)
                    sigx = 1.5 #ErrorData.iloc[j, FeH]
                    if (pd.isna(sigx)) or sigx > 0.5:
                        sigx = 1.5
                    sigy = 1.5 #ErrorData.iloc[j, element]
                    if pd.isna(sigy) or sigy > 0.5:
                        sigy = 1.5
                    correctBinFeH = findCorrectBin(logz_FeH, x0)
                    correctBinElement = findCorrectBin(logz_Element, y0)
                    #print correctBinFeH
                    #print correctBinElement
                    g2d = models.Gaussian2D(amplitude=1, #0.075
                                            x_mean=correctBinFeH,
                                            y_mean=correctBinElement,
                                            x_stddev=sigx,
                                            y_stddev=sigy)
                    z += g2d(x, y)

    #max = np.max(z)
    #z = np.divide(z, max)
    total = np.sum(z)
    z = np.divide(z, total)
    for ix, iy in np.ndindex(z.shape):
        if z[ix, iy] < 10e-9:
            z[ix, iy] = 0
    return z

def displayPlot(plot1, element, AbundanceData, title, ylabel, xlabel):
    plt.imshow(plot1,
               extent=(getMinBound("[Fe/H]", AbundanceData), getMaxBound("[Fe/H]", AbundanceData),
                       getMinBound(element, AbundanceData), getMaxBound(element, AbundanceData)),
               origin="lower", interpolation="gaussian", cmap = "Purples")
    plt.colorbar(spacing="proportional")
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)
    plt.show()

def createAverageScatterPlot(plot):
    rows = range(0, len(plot[0, :]))
    rowDimension = len(rows)
    averagePlot = np.zeros((len(rows), len(rows)))

    for i in range(0, rowDimension):
        if plot[:, i].sum() == 0:
            average = 0
        else:
            average = np.average(rows, weights=plot[:, i])
        average = int(np.round(average, 0))
        averagePlot[average][i] = 1

    return averagePlot

def createAverageLinePlot(plot, logz_Element):
    rows = range(0, len(plot[0, :]))
    rowDimension = len(rows)
    linePlot = np.zeros(len(rows))
    for i in range(0, rowDimension):
        if plot[:, i].sum() == 0:
            average = np.NAN
        else:
            average = np.average(rows, weights=plot[:, i])
            average = int(np.round(average, 0))
            average = logz_Element[average]
        linePlot[i] = average
    linePlot = smooth(linePlot, 10)
    ok = ~np.isnan(linePlot)
    xp = ok.ravel().nonzero()[0]
    fp = linePlot[~np.isnan(linePlot)]
    x = np.isnan(linePlot).ravel().nonzero()[0]
    linePlot[np.isnan(linePlot)] = np.interp(x, xp, fp)
    return linePlot

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def createStandardDeviationLinePlot(plot, logz_Element):
    rows = range(0, len(plot[0, :]))
    rowDimension = len(rows)
    std = np.zeros(len(rows))
    for i in range(0, rowDimension):
        stdValue = weighted_std(rows, plot[:, i])
        stdValue = int(np.round(stdValue, 0))
        stdValue = logz_Element[stdValue]
        std[i] = stdValue
    std = smooth(std, 30)
    ok = ~np.isnan(std)
    xp = ok.ravel().nonzero()[0]
    fp = std[~np.isnan(std)]
    x = np.isnan(std).ravel().nonzero()[0]
    std[np.isnan(std)] = np.interp(x, xp, fp)
    return std

def displayArray(array):
    array2 = np.matrix(array)
    plt.contourf(array2)
    plt.colorbar()
    plt.show()

def makeSmoothCurve(ylinePlot, xlinePlot, XaxisBins):
    x_new = np.linspace(0, len(xlinePlot), XaxisBins)
    a_BSpline = interpolate.make_interp_spline(xlinePlot, ylinePlot)
    y_new = a_BSpline(x_new)
    return y_new

def weighted_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    if weights.sum() == 0:
        average = np.zeros(len(values))
        variance = np.average((values - average) ** 2)
    else:
        average = np.average(values, weights=weights)
        variance = np.average((values - average) ** 2, weights=weights)

    return (math.sqrt(variance))

def showAllPlots(element, elementBrackets, Data, ErrorData):
    XaxisBins = 150
    '''x_minBound = -5.0
    x_maxBound = 0.5
    y_minBound = -1.5
    y_maxBound = 1.5'''

    x_minBound = getMinBoundFe(elementBrackets, Data) - 0.25
    x_maxBound = getMaxBoundFe(elementBrackets, Data) + 0.25
    y_minBound = getMinBound(elementBrackets, Data) - 0.25
    y_maxBound = getMaxBound(elementBrackets, Data) + 0.25

    #logz_FeH = np.linspace(getMinBound("[Fe/H]", Data), getMaxBound("[Fe/H]", Data), num=z_max_FeH)
    xBinNumber = 100
    xBins = np.linspace(x_minBound, x_maxBound, num=xBinNumber)
    #z_FeH = np.power(10, logz_FeH)
    #zero_FeH = (np.abs(z_FeH - 1)).argmin()
    #z_FeH[zero_FeH] = 1.0
    FeHAxis = np.linspace(x_minBound, x_maxBound, num=XaxisBins)

    yBinNumber = 100
    yBins = np.linspace(y_minBound, y_maxBound, num=yBinNumber)
    #z_Element = np.power(10, logz_Element)
    #zero_Element = (np.abs(z_Element - 1)).argmin()
    #z_Element[zero_Element] = 1.0


    ElementScatter = createPlotWithSources(Data, elementBrackets)
    ElementGaussian = createGaussianPlot(elementBrackets, Data, ErrorData, xBinNumber, yBinNumber, xBins, yBins)
    ElementAverageLinePlot = createAverageLinePlot(ElementGaussian, yBins)
    ElementAverageLinePlot = makeSmoothCurve(ElementAverageLinePlot, range(0, yBinNumber), XaxisBins)
    ElementSTDLinePlot = createStandardDeviationLinePlot(ElementGaussian, yBins)
    ElementSTDLinePlot = makeSmoothCurve(ElementSTDLinePlot, range(0, yBinNumber), XaxisBins)


    #norm = MidPointNorm(midpoint=0.0008) #0.0003
    im = plt.imshow(ElementGaussian,
               extent=(x_minBound, x_maxBound,
                       y_minBound, y_maxBound),
               origin="lower", interpolation="gaussian", cmap="Greys", norm = colors.PowerNorm(gamma = 0.45, vmin = 0.000001, vmax = ElementGaussian.max()))#cmap = bone #norm = norm
    plt.plot(FeHAxis, ElementAverageLinePlot, color = "red", linewidth = 1.5)
    plt.plot(FeHAxis, ElementAverageLinePlot - ElementSTDLinePlot, color = "red", linestyle = "dashed", linewidth = 0.5)
    plt.plot(FeHAxis, ElementAverageLinePlot + ElementSTDLinePlot, color = "red", linestyle = "dashed", linewidth = 0.5)
    plt.ylabel("[" + str(element) + "/Fe]")
    plt.xlabel('[Fe/H]')
    plt.title(elementBrackets)
    plt.legend(["Average", "STD", "STD", "Hill2019", "Kirby2009", "Starkenburg2013", "Tafelmeyer2018", "Kirby2018", "Duggan_2018", "Skuladottir_2017"], loc = "right", fancybox = True)
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [1, 2, 6, 3, 4, 5, 0]
    plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc = "upper right", fancybox = True)
    plt.xlim(x_minBound, x_maxBound)
    plt.ylim(y_minBound, y_maxBound)
    cbar = plt.colorbar()
    #plt.clim(10e-4, 0)
    cbar.set_label('Relative Contribution')
    #cbar.formatter = LogFormatterExponent(base=10, labelOnlyBase=True, minor_thresholds=(-1, -0.9))  # 10 is the default
    #cbar.ticker = LogLocator(base=10.0, subs=(1.0, ), numdecs=4, numticks=None)
    #cbar.ax.set_yticklabels(['-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '0'])
    cbar.update_ticks()
    plt.show()


class MidPointNorm(Normalize):
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint
            resdat[resdat>0] /= abs(vmax - midpoint)
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)

        if is_scalar:
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val-0.5)
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (value - 0.5)
            if val < 0:
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint

'''
# Go to correct column in Dataframe
# Iterate through column and go to corresponding [Fe/H] value
# Put the [Fe/H] value into correct "bin"
# After entire column has been iterated, average all the values in each "bin"
def createBins(AbundanceData, ErrorData, array, element, logz):
    Bins = {"Bin": logz, "Sum": 0, "Number": 0, "Mean": 0, "ErrorSum": 0, "ErrorMean": 0}
    BinsDF = pd.DataFrame(Bins)
    for i in AbundanceData.columns:
        if i == element:
            element = AbundanceData.columns.get_loc(i)
            FeH = AbundanceData.columns.get_loc("[Fe/H]")
            for j in range(0, len(AbundanceData.index)):
                if pd.notna(AbundanceData.iloc[j, element]):
                    correctBin = findCorrectBin(array, AbundanceData.iloc[j, FeH])
                    BinsDF.iloc[correctBin, 5] += AbundanceData.iloc[j, element]
                    BinsDF.iloc[correctBin, 4] += 1
                    BinsDF.iloc[correctBin, 2] += ErrorData.iloc[j, element]

            BinsDF = calculateAverages(BinsDF)
            return BinsDF
'''

'''
# Creates plot for all isotopic abundance vs metallicity graphs, and can also take an isotope as
# parameter and plot Isotopic Abundance vs Metallicity for that element. Also shows source of data, seperates data
# into colors, and shows a mean line
def makePlotWithSourcesAndMeanLine(Data, ErrorData, logz_FeH, element=None, ElementBin=None):
    if element is None:
        for i in Data.columns:
            if Data.columns.get_loc(i) > 0:
                if i != '[Fe/H]':
                    Bins = createBins(Data, ErrorData, logz_FeH, i)

                    hill_plot = Hill_2019.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i), color='blue',
                                                       label="Hill_2019", figsize = (15, 8))
                    kirby_plot = Kirby_2009.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                         color='green',
                                                         label="Kirby_2009", ax=hill_plot, figsize = (15, 8))
                    starkenburg_plot = Starkenburg_2013.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                                     color='black', label="Starkenburg_2013",
                                                                     ax=hill_plot, figsize = (15, 8))
                    tafelmeyer_plot = Tafelmeyer_2018.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                                   color='red', label="Tafelmeyer_2018",
                                                                   ax=hill_plot, figsize = (15, 8))
                    AveragePlot = Bins.plot(x='Bin', y='Mean', color='purple', ax=hill_plot, figsize = (15, 8), linewidth = 2.0)

                    plt.show(block=True)
    else:
        for i in Data.columns:
            if i == element:
                ElementBin = createBins(Data, ErrorData, logz_FeH, i)
                hill_plot = Hill_2019.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i), color='blue',
                                                   label="Hill_2019", figsize = (15, 8))
                kirby_plot = Kirby_2009.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i), color='green',
                                                     label="Kirby_2009", ax=hill_plot, figsize = (15, 8))
                starkenburg_plot = Starkenburg_2013.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                                 color='black', label="Starkenburg_2013",
                                                                 ax=hill_plot, figsize = (15, 8))
                tafelmeyer_plot = Tafelmeyer_2018.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i),
                                                               color='red', label="Tafelmeyer_2018",
                                                               ax=hill_plot, figsize = (15, 8))
                AveragePlot = ElementBin.plot(x='Bin', y='Mean', color='purple', ax=hill_plot, figsize = (15, 8), linewidth = 2.0, label = "Mean")
                ErrorPlot = ElementBin.plot(x='Bin', y='ErrorMean', color='pink', ax=hill_plot, figsize = (15, 8), linewidth = 2.0, label = "Error")
                plt.show(block=True)    
'''