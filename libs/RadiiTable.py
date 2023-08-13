import pandas as pd
from os.path import join, dirname

RadiiTable = pd.read_csv(join(dirname(__file__), "..", "assets", "Radii.csv"))

default_radius_table = { element_num: radius for (element_num, radius) in zip(RadiiTable["Element"], RadiiTable["r"])}
