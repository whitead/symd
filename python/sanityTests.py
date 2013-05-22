from SimpleMD import SimpleMD
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import probplot
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

class Test:

    def print_summary(self):
        if(self.passed):
            print "%s:Passed" % self.name
        else:
            print "%s" % self.name
            print self.summary
            for p in self.plots:
                print "![%s](%s \"%s\")" % (p.caption, p.path, p.name)
                print p.caption
                print "\n"


class Plot:
    def __init__(self, name, path, caption="", size=[8,6]):
        self.name = name
        self.path = path
        self.caption = caption
        self.size = size


class SanityTester:

    def __init__(self, md, name):
        self.tests = []
        self.add(EnergyConservationTest(md, name))
        self.add(EnergyConservationTest(md, name))

    def add(self, test):
        self.tests.append(test)

    def run(self):
        for t in self.tests:
            success = t.run()
            t.print_summary()
            if(not success):
                return

class EnergyConservationTest(Test):

    def __init__(self, md, name):
        self.short_name = "Energy_%s" % name
        self.name = "Energy Conservation [%s]" % name
        self.md = md
        self.passed = False

        self.dtol = 1e-10
        self.ptol = 0.01
        self.equil_frac = 0.3

    def run(self):

        if(not md.executed):
            md.execute()

        energy = md.te[round(len(md.te) * self.equil_frac):]
        mean = np.mean(energy)
        dev = np.std(energy)
        
        if(dev < self.dtol):
            self.passed = True
            self.summary = "Energy fluctation is %g" % dev
            return self.passed

        norm_energy = (energy - mean) / dev 
        cor = pearsonr(range(len(energy)), norm_energy)

        if(abs(cor[1]) < self.ptol):
            self.passed = True

        self.summary = "Time Correlation = %g (p-value: %g)" % cor

        self.plots = [Plot("%s Q-Q Norm Plot" % self.name, "%s_energy_qq.png" % self.short_name)]

        fig = plt.figure(figsize=(12, 6))

        ax = fig.add_subplot(121)
        n, bins, patches = ax.hist(norm_energy,  50, normed=1, facecolor='yellow', alpha=0.75)

        ax.set_xlabel("Normalized Energy")
        ax.set_ylabel("Density")

        ax.grid(True)

        ax = fig.add_subplot(122)

        ((osm, osr), (slope, intercept, r)) = probplot(norm_energy, fit=True)

        ax.plot(osm, osr, 'o', osm, slope*osm + intercept)
        ax.set_title('Probability Plot')
        ax.set_xlabel('Quantiles')
        ax.set_ylabel('Ordered Values')

        xmin = np.amin(osm)
        xmax = np.amax(osm)
        ymin = np.amin(norm_energy)
        ymax = np.amax(norm_energy)
        posx = xmin + 0.70 * (xmax - xmin)
        posy = ymin + 0.01 * (ymax - ymin)

        plt.savefig(self.plots[0].path)

        return self.passed

                  

md = SimpleMD(1000, 3, temperature=0.7, exeDir="test")
md.setup_positions(0.7)
md.setup_masses(1)

runner = SanityTester(md, "LJ_NVE")
runner.run()
