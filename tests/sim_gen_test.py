import os
from c2e2 import *


TOTALV2 = os.path.dirname(os.path.abspath(__file__)) + "/../sample_hyxml/TotalMotionV2.hyxml"
linThermo = os.path.dirname(os.path.abspath(__file__)) + "/../sample_hyxml/linThermo.hyxml"


# autos, props = FileHandler.load_model(TOTALV2) 
autos, props = FileHandler.load_model(linThermo) 

prop = props[0]
sim = Simulator(automata=autos, prop=prop, time_bound=10)


sim.simulate()