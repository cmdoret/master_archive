import matplotlib.pyplot as plt
import pandas as pd

"""
This script generates plots from custom generated summaries of the BWA
mapping output with different parameters. This is useful to visualize which
combinations of parameters yield the best results.
Cyril Matthey-Doret
02.04.2017
"""
map_sum = pd.read_csv("mapstats.csv",header=1)  # Loading mapping summaries
