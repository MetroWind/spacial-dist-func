#!/usr/bin/env python3

import sys, os
import json

import matplotlib
import matplotlib.pyplot as Plt

with open("hist.json", 'r') as f:
    Data = json.load(f)

Fig = Plt.figure()                      # Create a `figure' instance
Ax = Fig.add_subplot(111)               # Create a `axes' instance in the figure
Colors = Ax.pcolormesh(Data['x'], Data['y'], Data['c'], cmap="BuPu")
Ax.set_aspect('equal', adjustable='box')
Fig.colorbar(Colors, ax=Ax)

for Atom in Data["specials"]:
    Circle = matplotlib.patches.Circle(Data["specials"][Atom], 0.01)
    Ax.add_patch(Circle)
    Ax.annotate(Atom, Data["specials"][Atom])

Fig.savefig("test.pdf")
