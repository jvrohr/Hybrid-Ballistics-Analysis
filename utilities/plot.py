import matplotlib.pyplot as plt

class PlotResults:
    def __init__(self):
        self.results_dict = {}
        # Initialize results dict structure
        self.results_dict["Time"] = []
        self.results_dict["Temperature"] = []
        self.results_dict["Quantity Gas"] = []
        self.results_dict["Quantity Liquid"] = []
        self.results_dict["Thrust"] = []
        self.results_dict["Pressure Chamber"] = []
        self.results_dict["Pressure Tank"] = []
        self.results_dict["Oxidizer Mass Flow"] = []

    def plot_results_burn(self):
        title = "Hybrid Motor Burn Result"
        time = self.results_dict["Time"]
        self.plot_data(time, self.results_dict["Thrust"], "Thrust [N]", title)
        self.plot_data(time, [i/1e5 for i in self.results_dict["Pressure Chamber"]], "Pressure Chamber [bar]", title)
        self.plot_data(time, [i/1e5 for i in self.results_dict["Pressure Tank"]], "Tank Pressure [bar]", title)

    def plot_results_blowdown(self):
        title = "Hybrid Motor Blowdown Result"
        time = self.results_dict["Time"]
        self.plot_data(time, self.results_dict["Temperature"], "Temperature [K]", title)
        self.plot_data(time, self.results_dict["Quantity Gas"], "Quantity Gas [mols]", title)
        self.plot_data(time, self.results_dict["Quantity Liquid"], "Quantity Liquid [mols]", title)
        self.plot_data(time, [i/1e5 for i in self.results_dict["Pressure Tank"]], "Tank Pressure [bar]", title)
        self.plot_data(time, [i/1e5 for i in self.results_dict["Pressure Chamber"]], "Pressure Chamber [bar]", title)
        self.plot_data(time, self.results_dict["Oxidizer Mass Flow"], "Oxidizer Mass Flow [kg/s]", title)

    def plot_data(self, X, Y, ylabel, title, xlabel="Time [s]", show=True):
        plt.plot(X, Y)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid()
        plt.title(title)
        if show:
            plt.show()