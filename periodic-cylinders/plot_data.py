import csv

import numpy as np
import matplotlib.pyplot as plt


class Table():
    def __init__(self, filename):
        self.load_data(filename)

    def load_data(self, filename):
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            radii = []
            phases = []
            for i, row in enumerate(reader):
                if i == 0:
                    gaps = row[1:]
                elif i == 1:
                    continue
                else:
                    radii.append(row[0])
                    phases.append(row[1:])
            self.radii = np.array(radii, dtype=np.int64)
            self.gaps = np.array(gaps, dtype=np.int64)
            self.phases = np.array(phases, dtype=np.float64)

    def get_phases(self, radii, gaps):
        rows = [np.where(self.radii == r)[0][0] for r in radii]
        rows = np.array([[i] for i in rows])
        cols = [np.where(self.gaps == g)[0][0] for g in gaps]
        cols = np.array(cols)
        return self.phases[rows, cols]


def main():
    table = Table('t155-lambda633-sim-results-Phase.csv')

    radii = list(range(130, 181))
    gaps = [46]
    phases_std = table.get_phases(radii, gaps)
    phases = []
    with open('phi.txt', 'r') as f:
        for line in f:
            phases.append(float(line.rstrip()))
    phases = (np.array(phases)/np.pi) % 2 - 1.0


    plt.plot(radii, phases_std, 'r', label='theirs')
    plt.plot(radii, phases, 'b', label='ours')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
