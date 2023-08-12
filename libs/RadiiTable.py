import pandas as pd
from os.path import join, dirname

RadiiTable = pd.read_csv(join(dirname(__file__), "..", "assets", "Radii.csv"))
