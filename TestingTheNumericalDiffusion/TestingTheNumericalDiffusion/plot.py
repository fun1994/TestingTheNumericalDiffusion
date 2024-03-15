# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 09:06:04 2024

@author: HFKJ059
"""

import numpy as np
from matplotlib import pyplot as plt


def read_1d(path, filename):
    with open("./data/" + path + "/" + filename + ".txt", "r") as file:
        data = file.read()
    data = data.split()
    for i in range(len(data)):
        data[i] = float(data[i])
    data = np.array(data)
    return data

def read_2d(path, filename):
    data = []
    with open("./data/" + path + "/" + filename + ".txt", "r") as file:
        while True:
            line = file.readline()
            if not line:
                break
            data_temp = line.split()
            data.append(data_temp)
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j] = float(data[i][j])
    data = np.array(data)
    return data

def read(path, convection):
    x = read_1d(path, "x")
    y = read_1d(path, "y")
    phi = read_2d(path, "phi_" + convection)
    return x, y, phi

def plot(x, y, phi, title):
    X, Y = np.meshgrid(x, y)
    plt.contourf(X, Y, phi.T, levels=1000, cmap="jet")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)
    plt.colorbar()
    plt.show()

def run(N, convection):
    x, y, phi = read("N=" + N, convection)
    plot(x, y, phi, "N=" + N + ", " + convection)

def main():
    run("20", "UDS")
    run("20", "blend")
    run("40", "UDS")
    run("40", "blend")
    run("80", "UDS")
    run("80", "blend")


main()
