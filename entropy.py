import math
import scipy
import matplotlib.pyplot as plt
from scipy.stats import *

O = {'A': 119.0, 'B': 66.0, 'C': 720.0}

D = {'A': 112.0, 'B': 51.0, 'C': 742.0}

fee = {('A', 'A'): 3.0, ('A', 'B'): 4.0, ('A', 'C'): 9.0,
       ('B', 'A'): 4.0, ('B', 'B'): 3.0, ('B', 'C'): 6.0,
       ('C', 'A'): 9.0, ('C', 'B'): 6.0, ('C', 'C'): 3.0}

data = {('A', 'A'): 50.0, ('A', 'B'): 20.0, ('A', 'C'): 100.0,
       ('B', 'A'): 20.0, ('B', 'B'): 7.0, ('B', 'C'): 6.0,
       ('C', 'A'): 1.0, ('C', 'B'): 8.0, ('C', 'C'): 1000.0}

class Test:
    def __init__(self, O, D, fee, data):
        self.O = O
        self.D = D
        self.fee = fee
        self.beta = 0.948
        self.data = data

    def calPQ(self):
        count = 0
        P = {}
        Q = {}
        for i in self.O.keys():
            Q[i] = 1
        while count < 100:
            for i in self.O.keys():
                divider = 0
                for j in self.D.keys():
                    divider += Q[j] * self.D[j] * math.exp(- self.beta * self.fee[(i, j)])
                P[i] = 1.0 / divider
            for j in self.D.keys():
                divider = 0
                for i in self.O.keys():
                    divider += P[i] * self.O[i] * math.exp(- self.beta * self.fee[(i, j)])
                Q[j] = 1.0 / divider
            count += 1
        return P, Q

    def calT(self):
        T = {}
        P, Q = self.calPQ()
        for i in self.O.keys():
            for j in self.D.keys():
                T[(i, j)] = P[i] * Q[j] * self.O[i] * self.D[j] * math.exp(- self.beta * self.fee[(i, j)])
        return T

    def normalize(self, data):
        for i in data.keys():
            list = []
            divider = sum(data[i])
            for j in data[i]:
                list.append(j / divider)
            data[i] = list
        return data

    def poisson(self, update, count):
        T = self.calT()
        distribution_x = {}
        distribution_y = {}
        for i in T.keys():
            upper_limit = min(self.O[i[0]], self.D[i[1]])
            distribution_x[i] = list(range(0, int(upper_limit) + 1))
        if update == False:
            for i in distribution_x.keys():
                probability_list = []
                for j in distribution_x[i]:
                    probability = distributions.poisson.pmf(j, T[i])
                    probability_list.append(probability)
                distribution_y[i] = probability_list
        else:
            for i in distribution_x.keys():
                probability_list = []
                for j in distribution_x[i]:
                    probability = distributions.poisson.pmf(j, T[i]) * distributions.poisson.pmf(self.data[i], j)**count
                    probability_list.append(probability)
                distribution_y[i] = probability_list

        return distribution_x, self.normalize(distribution_y)
        

    def makeGraph(self, update, count = 1):
        distribution_x, distribution_y = self.poisson(update, count)
        fig = plt.figure()
        ax1 = fig.add_subplot(3, 3, 1)
        ax2 = fig.add_subplot(3, 3, 2)
        ax3 = fig.add_subplot(3, 3, 3)
        ax4 = fig.add_subplot(3, 3, 4)
        ax5 = fig.add_subplot(3, 3, 5)
        ax6 = fig.add_subplot(3, 3, 6)
        ax7 = fig.add_subplot(3, 3, 7)
        ax8 = fig.add_subplot(3, 3, 8)
        ax9 = fig.add_subplot(3, 3, 9)

        ax1.bar(distribution_x[('A', 'A')], distribution_y[('A', 'A')])
        ax2.bar(distribution_x[('A', 'B')], distribution_y[('A', 'B')])
        ax3.bar(distribution_x[('A', 'C')], distribution_y[('A', 'C')])
        ax4.bar(distribution_x[('B', 'A')], distribution_y[('B', 'A')])
        ax5.bar(distribution_x[('B', 'B')], distribution_y[('B', 'B')])
        ax6.bar(distribution_x[('B', 'C')], distribution_y[('B', 'C')])
        ax7.bar(distribution_x[('C', 'A')], distribution_y[('C', 'A')])
        ax8.bar(distribution_x[('C', 'B')], distribution_y[('C', 'B')])
        ax9.bar(distribution_x[('C', 'C')], distribution_y[('C', 'C')])
        plt.show()


                

if __name__ == '__main__':            
    test = Test(O, D, fee, data)
    test.makeGraph(False)
    test.makeGraph(True)
    test.makeGraph(True, 2)


            
            

