import struct

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def load_buf(filename):
    with open(filename, mode='rb') as f:
        content = f.read()
        data = [x for x in struct.iter_unpack('d', content)]
    return np.array(data)


def get_timeseries(index, firststep, laststep, step):
    series = []
    for i in range(firststep, laststep, step):
        filename = 'vtk/hx{:07d}.dat'.format(i)
        hx = load_buf(filename)
        series.append(hx[index])
    return np.squeeze(np.asarray(series))


class FittedSine():
    def __init__(self, tdata, hdata, kz, z, omega):
        self.tdata, self.hdata = tdata, hdata
        self.kz, self.z, self.omega = kz, z, omega
        self.fit()

    def f(self, t, A, phi):
        return A*np.sin(-self.kz*self.z - self.omega*t - phi)

    def jac(self, t, A, phi):
        out = np.empty((t.size, 2))
        arg = -self.kz*self.z - self.omega*t - phi
        out[:,0] = np.sin(arg)
        out[:,1] = -A*np.cos(arg)
        return out

    def fit(self):
        popt, _ = curve_fit(self.f, self.tdata, self.hdata)
        self.A, self.phi = popt

    def __call__(self, t):
        return self.f(t, self.A, self.phi)

    def __repr__(self):
        info = []
        info.append('A = {}'.format(self.A))
        info.append('phi = {}'.format(self.phi))
        return '\n'.join(info)


def get_amplitude_and_phi(index, firststep, laststep, step):
    dt = 8.49e-3
    z = -7.104
    omega = 0.993
    eps = 1.0
    mu = 1.0
    kz = np.sqrt(eps*mu)*omega

    hx = get_timeseries(index, firststep, laststep, step)
    times = dt*np.arange(firststep, laststep, step)
    hxinc = np.sin(-kz*z - omega*times)
    fit = FittedSine(times, hx, kz, z, omega)
    return fit.A, fit.phi


def main():
    firststep = 9000
    laststep = 10575
    step = 25




    x = load_buf('vtk/xcoordinates.dat')
    y = load_buf('vtk/ycoordinates.dat')

    phi = []
    start = x.size//2
    end = x.size//2+1
    for index in range(start, end):
        hx = get_timeseries(index, firststep, laststep, step)
        times = dt*np.arange(firststep, laststep, step)
        hxinc = np.sin(-kz*z - omega*times)
        fit = FittedSine(times, hx, kz, z, omega)
        print("On step {}/{}, phi = {}".format(index - start + 1,
                                               end - start, fit.phi))
        phi.append(fit.phi)
    phi = np.asarray(phi)
    hist, bins = np.histogram(phi)
    width = 0.7*(bins[1] - bins[0])
    center = (bins[:-1] + bins[1:])/2
    plt.bar(center, hist, align='center', width=width)
    plt.show()

    plt.plot(times, hxinc, 'b', label='incident')
    plt.plot(times, hx, 'r', label='scattered')
    plt.plot(times, fit(times), 'g', label='fit')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
