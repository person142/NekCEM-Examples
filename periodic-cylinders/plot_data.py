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

    radii = [149, 150, 151, 152, 153, 154, 155, 156, 157, 158]
    gaps = [46]
    phases_std = table.get_phases(radii, gaps)
    phases = [0.9467198592747457, 0.9722228421559183,
              1.00281088750689, 1.0368108556864222,
              1.0727982668931104, 1.1093315970881696,
              1.1451847086635345, 1.1787679451151407,
              1.2090714764691635, 1.235163281494783]

    plt.plot(radii, phases_std, 'r', label='theirs')
    plt.plot(radii, phases, 'b', label='ours')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
