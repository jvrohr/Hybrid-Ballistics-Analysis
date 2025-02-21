import matplotlib.pyplot as plt

class PlotResults:
    def __init__(self):
        self.resultsDict = {}
        # Initialize results dict structure
        self.resultsDict["Time"] = []
        self.resultsDict["Temperature"] = []
        self.resultsDict["Quantity Gas"] = []
        self.resultsDict["Quantity Liquid"] = []
        self.resultsDict["Thrust"] = []
        self.resultsDict["Pressure Chamber"] = []
        self.resultsDict["Pressure Tank"] = []
        self.resultsDict["Oxidizer Mass Flow"] = []

    def PlotResultsBurn(self):
        title = "Hybrid Motor Burn Result"
        time = self.resultsDict["Time"]
        self.PlotData(time, self.resultsDict["Thrust"], "Thrust [N]", title)
        self.PlotData(time, [i/1e5 for i in self.resultsDict["Pressure Chamber"]], "Pressure Chamber [bar]", title)
        self.PlotData(time, [i/1e5 for i in self.resultsDict["Pressure Tank"]], "Tank Pressure [bar]", title)

    def PlotResultsBlowdown(self):
        title = "Hybrid Motor Blowdown Result"
        time = self.resultsDict["Time"]
        self.PlotData(time, self.resultsDict["Temperature"], "Temperature [K]", title)
        self.PlotData(time, self.resultsDict["Quantity Gas"], "Quantity Gas [mols]", title)
        self.PlotData(time, self.resultsDict["Quantity Liquid"], "Quantity Liquid [mols]", title)
        self.PlotData(time, [i/1e5 for i in self.resultsDict["Pressure Tank"]], "Tank Pressure [bar]", title)
        self.PlotData(time, [i/1e5 for i in self.resultsDict["Pressure Chamber"]], "Pressure Chamber [bar]", title)
        self.PlotData(time, self.resultsDict["Oxidizer Mass Flow"], "Oxidizer Mass Flow [kg/s]", title)

    def PlotData(self, X, Y, ylabel, title, xlabel="Time [s]", show=True):
        plt.plot(X, Y)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid()
        plt.title(title)
        if show:
            plt.show()