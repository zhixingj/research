import pandas as pd
import matplotlib.pyplot as plt
import math as math
import numpy as np
from astropy.modeling import models
from matplotlib import pyplot as plt, pylab
from matplotlib.pyplot import *
from scipy.interpolate import interpolate
from datetime import date
from gce_data import *
from GCEJuly17 import *

# np.set_printoptions(threshold=np.inf)
#
# pd.set_option("display.precision", 4)

source = pd.read_csv(test(), sep = ",", names=["[Fe/H]", "[Mg/Fe]", "mgfe_sig", "feh_sig"])
#source = pd.read_csv('aa_.csv', sep = ",", names=['[Mg/Fe]', "g", "f"])
modelFile = pd.read_csv(gce(), sep = ",", names=['[Fe/H]', "[Mg/Fe]"])

# print(modelFile.at[0, "[Mg/Fe]"])

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
# makePlot(modelFile, "[Mg/Fe]")
#Creates plot for all isotopic abundance vs metallicity graphs, and can also take an isotope as parameter and plot Isotopic Abundance vs Metallicity for that element. Also shows source of data, and seperates data into colors
def createPlotWithSources(df, element=None):
    if element is None:
        for i in df.columns:
            # print i
            if df.columns.get_loc(i) > 0:
                if i != '[Fe/H]':
                    source_plot = source.plot.scatter(x='[Fe/H]', y=source.columns.get_loc(i), color='blue',
                                                       marker='.',
                                                       label="source", extent=(getMinBound("[Fe/H]", df), getMaxBound("[Fe/H]", df),
                                                       getMinBound(element, df), getMaxBound(element, df)))
                    # modelFile_plot = Kirby_2009.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i), color='green',
                    #                                      marker="*",
                    #                                      label="Kirby_2009", ax=hill_plot)

                    plt.xlim(getMinBound("[Fe/H]", df) - 1, getMaxBound("[Fe/H]", df) + 1)
                    plt.ylim(getMinBound(element, df) - 1, getMaxBound(element, df) + 1)
                    return source_plot
    else:
        for i in df.columns:
            if i == element:
                source_plot = source.plot.scatter(x='[Fe/H]', y=source.columns.get_loc(i), color='blue',
                                                   marker='.',
                                                   label="source")
                # kirby_plot = Kirby_2009.plot.scatter(x='[Fe/H]', y=Hill_2019.columns.get_loc(i), color='green',
                #                                      marker="*",
                #                                      label="Kirby_2009", ax=hill_plot)

                plt.xlim(getMinBound("[Fe/H]", df) - 1, getMaxBound("[Fe/H]", df) + 1)
                plt.ylim(getMinBound(element, df) - 1, getMaxBound(element, df) + 1)
                return source_plot


# def showData(df, element=None):
#     if element is None:
#         print (df.to_string());
#     else:
#         for i in df.columns:
#             if i == element:
#                 print(df.loc[:, element])
#
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

# def calculateAverages(df):
#     for column in df.columns:
#         if column == "Mean":
#             df.loc[:, "Mean"] = df.loc[:, "Sum"] / df.loc[:, "Number"]
#         if column == "ErrorMean":
#             df.loc[:, "ErrorMean"] = df.loc[:, "ErrorSum"] / df.loc[:, "Number"]
#     return df
#
# creates a 2d plot by applying a gaussian spread function using metallicity abundance values and corresponding error values
# # def createGaussianPlot(element, AbundanceData, ErrorData, z_max_FeH, z_max_Element, logz_FeH, logz_Element):




def createGaussianPlot(element, AbundanceData, z_max_FeH, z_max_Element, logz_FeH, logz_Element):

    FeH = AbundanceData.columns.get_loc("[Fe/H]")
    mgfe_sig = AbundanceData.columns.get_loc("mgfe_sig")
    feh_sig = AbundanceData.columns.get_loc("feh_sig")
    rows = z_max_Element
    columns = z_max_FeH


    z = np.zeros((rows, columns))
    x, y = np.meshgrid(np.arange(0, rows), np.arange(0, columns), indexing="xy")

    for i in AbundanceData.columns:
        if i == element:
            element = AbundanceData.columns.get_loc(i)
            # for j in range (0,2):
            for j in range(0, len(AbundanceData.index)):
                if pd.notna(AbundanceData.iloc[j, element]):
                    x0 = AbundanceData.iloc[j, FeH]
                    y0 = AbundanceData.iloc[j, element]
                    sigx = AbundanceData.iloc[j, feh_sig]
                    if (pd.isna(sigx)) or sigx > 0.5:
                        sigx = 1.5
                    sigy = AbundanceData.iloc[j, mgfe_sig]
                    if pd.isna(sigy) or sigy > 0.5:
                        sigy = 1.5
                    correctBinFeH = findCorrectBin(logz_FeH, x0) - 1
                    correctBinElement = findCorrectBin(logz_Element, y0)
                    g2d = models.Gaussian2D(amplitude= 1,
                                            x_mean=correctBinFeH,
                                            y_mean=correctBinElement,
                                            x_stddev=sigx,
                                            y_stddev=sigy)
                    # print(g2d(x, y)+z)
                    z = z + g2d(x, y)


                    # z += g2d(x, y)

    return z

XaxisBins = 300
z_max_FeH = 6
Data = source
elementBrackets = "[Mg/Fe]"
logz_FeH = np.linspace(getMinBound("[Fe/H]", Data), getMaxBound("[Fe/H]", Data), num=z_max_FeH)
z_FeH = np.power(10, logz_FeH)
zero_FeH = (np.abs(z_FeH - 1)).argmin()
z_FeH[zero_FeH] = 1.0
FeHAxis = np.linspace(getMinBound("[Fe/H]", Data), getMaxBound("[Fe/H]", Data), num=XaxisBins)

z_max_Element = 6
logz_Element = np.linspace(getMinBound(elementBrackets, Data), getMaxBound(elementBrackets, Data), num=z_max_Element)
z_Element = np.power(10, logz_Element)
zero_Element = (np.abs(z_Element - 1)).argmin()
z_Element[zero_Element] = 1.0

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
#
bb = createGaussianPlot("[Mg/Fe]", source, z_max_FeH, z_max_Element, logz_FeH, logz_Element)
# displayPlot(bb, "[Mg/Fe]", source, "ff", "Mg", "H")
# makePlot(modelFile, "[Mg/Fe]")
makePlot(source, "[Mg/Fe]")
# def createAverageScatterPlot(plot):
#     rows = range(0, len(plot[0, :]))
#     rowDimension = len(rows)
#     averagePlot = np.zeros((len(rows), len(rows)))
#
#     for i in range(0, rowDimension):
#         if plot[:, i].sum() == 0:
#             average = 0
#         else:
#             average = np.average(rows, weights=plot[:, i])
#         average = int(np.round(average, 0))
#         averagePlot[average][i] = 1
#
#     return averagePlot

def createAverageLinePlot(plot, logz_Element):
    rows = range(0, len(plot[0, :]))
    rowDimension = len(rows)
    linePlot = np.zeros(len(rows))
    for i in range(0, rowDimension):
        if plot[:, i].sum() == 0:
            average = 0
        else:
            average = np.average(rows, weights=plot[:, i])
        average = int(np.round(average, 0))
        average = logz_Element[average]
        linePlot[i] = average
    return linePlot

def createStandardDeviationLinePlot(plot, logz_Element):
    rows = range(0, len(plot[0, :]))
    rowDimension = len(rows)
    std = np.zeros(len(rows))
    for i in range(0, rowDimension):
        stdValue = weighted_std(rows, plot[:, i])
        stdValue = int(np.round(stdValue, 0))
        stdValue = logz_Element[stdValue]
        std[i] = stdValue
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

# def showAllPlots(element, elementBrackets, Data, ErrorData):
# def showAllPlots(element, elementBrackets, Data):
#     XaxisBins = 300
#     z_max_FeH = 6
#     logz_FeH = np.linspace(getMinBound("[Fe/H]", Data), getMaxBound("[Fe/H]", Data), num=z_max_FeH)
#     z_FeH = np.power(10, logz_FeH)
#     zero_FeH = (np.abs(z_FeH - 1)).argmin()
#     z_FeH[zero_FeH] = 1.0
#     FeHAxis = np.linspace(getMinBound("[Fe/H]", Data), getMaxBound("[Fe/H]", Data), num=XaxisBins)
#
#     z_max_Element = 6
#     logz_Element = np.linspace(getMinBound(elementBrackets, Data), getMaxBound(elementBrackets, Data), num=z_max_Element)
#     z_Element = np.power(10, logz_Element)
#     zero_Element = (np.abs(z_Element - 1)).argmin()
#     z_Element[zero_Element] = 1.0
#
#     ElementScatter = createPlotWithSources(Data, elementBrackets)
#     ElementGaussian = createGaussianPlot(elementBrackets, Data, z_max_FeH, z_max_Element, logz_FeH, logz_Element)
#     ElementAverageLinePlot = createAverageLinePlot(ElementGaussian, logz_Element)
#     ElementAverageLinePlot = makeSmoothCurve(ElementAverageLinePlot, range(0, z_max_FeH), XaxisBins)
#     ElementSTDLinePlot = createStandardDeviationLinePlot(ElementGaussian, logz_Element)
#     ElementSTDLinePlot = makeSmoothCurve(ElementSTDLinePlot, range(0, z_max_FeH), XaxisBins)
#
#     '''plt.imshow(ElementGaussian,
#                extent=(getMinBound("[Fe/H]", Data), getMaxBound("[Fe/H]", Data),
#                        getMinBound(elementBrackets, Data), getMaxBound(elementBrackets, Data)),
#                origin="lower", interpolation="gaussian", cmap="Greys")'''
#     plt.imshow(ElementGaussian,
#                extent=(getMinBound("[Fe/H]", Data) - 1, getMaxBound("[Fe/H]", Data) + 1,
#                        getMinBound(elementBrackets, Data) - 1, getMaxBound(elementBrackets, Data) + 1),
#                origin="lower", interpolation="gaussian", cmap="Greys")
#     plt.plot(FeHAxis, ElementAverageLinePlot, color = "red", linewidth = 1.5)
#     plt.plot(FeHAxis, ElementAverageLinePlot - ElementSTDLinePlot, color = "red", linestyle = "dashed", linewidth = 0.5)
#     plt.plot(FeHAxis, ElementAverageLinePlot + ElementSTDLinePlot, color = "red", linestyle = "dashed", linewidth = 0.5)
#     plt.ylabel("[" + str(element) + "/Fe]")
#     plt.xlabel('[Fe/H]')
#     plt.title(elementBrackets)
#     #plt.legend(["Average", "STD", "STD", "Hill2019", "Kirby2009", "Starkenburg2013", "Tafelmeyer2018", "Kirby2018"], loc = "right", fancybox = True)
#     handles, labels = plt.gca().get_legend_handles_labels()
#     order = [1, 2, 3, 4, 0]
#
#
#     #plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc = "right", fancybox = True)
#     plt.xlim(getMinBound("[Fe/H]", Data) - 1, getMaxBound("[Fe/H]", Data) + 1.5)
#     plt.ylim(getMinBound(elementBrackets, Data) - 1, getMaxBound(elementBrackets, Data) + 1)
#     cbar = plt.colorbar()
#     cbar.ax.set_yticklabels(['-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '0'])
#     cbar.set_label('log (Relative Contribution)')
#     plt.show()
#
# showAllPlots("[Mg/Fe]", "[Mg/Fe]" , source)
