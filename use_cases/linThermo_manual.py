from c2e2 import *

# Manually construct a simulator to see if
# it's different from the one loaded from hyxml

# ================= Automaton ======================
# Declare a Automaton
linThermo = Automaton("linThermo")

# Add Variables to the Automaton
x = Variable(name="x", type=REAL, scope=LOCAL)
linThermo.add_variable(x)

# ================== Mode ===========================
# Heating Mode
heating = Mode(name="heating", id=0, isInitial=True)
# Cooling Mode
cooling = Mode(name="cooling", id=1, isInitial=False)

# Dynamics of each mode
h1 = DAI("x_dot = 40-0.5*x")
c1 = DAI("x_dot = 30-0.5*x")

# Invariant
h_inv = Invariant("x<75")
c_inv = Invariant("x>65")

# Initial Set
init_set = RectangleSet("cooling: x>=68 && x<=72")
unsafe_set = RectangleSet("x<=63")

# ================== Transition ======================
# guard 
h2c_guard = Guard("x>=75")
c2h_guard = Guard("x<=65")
heating_to_cooling = Transition(guard=h2c_guard, id=0, src_mode=0, dst_mode=1)
cooling_to_heating = Transition(guard=c2h_guard, id=1, src_mode=1, dst_mode=0)

# ================== Parameters ======================
time_step = 0.01
time_bound = 10

# ================== Add Everything Together =========
# Add the differential equations/invariants
heating.add_dai(h1)
cooling.add_dai(c1)

heating.add_invariant(h_inv)
cooling.add_invariant(c_inv)

# Add the mode/transition to Automaton
linThermo.add_mode(heating)
linThermo.add_mode(cooling)

linThermo.add_transition(heating_to_cooling)
linThermo.add_transition(cooling_to_heating)

# ================== Simulator ========================
# Simulator
sim = Complete_Automata(automata=linThermo,
                        # initial_set=init_set,
                        # unsafe_set=unsafe_set,
                        initial_state=[68,],
                        initial_mode="cooling",
                        time_step=time_step,
                        time_bound=time_bound,
                        name="linear_thermo",
                        unsafe_check=False)

# simualte
data_points = sim.simulate()

print(len(data_points))