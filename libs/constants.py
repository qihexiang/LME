import os

EPS = 1e-8
PRODUCTION = os.getenv("LME_PRODUCTION") in ["Y", "y", "YES", "yes", 1]
