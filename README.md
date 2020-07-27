
# Functionalities & Todos
- Overall:
  - [x] Be able to load model from .hyxml file.

- Automaton():
  - [ ] ...

- RectangleSet()
  - [x] RectangleSet(): Can be used to describe Initial Set & Unsafe Set
  - [ ] Have InitialSet/UnsafeSet class inherit from RectangleSet() in order to declare those sets more easily.
  
- DAI Equations()
  - [x] DAI(): Can be used to describe differential equations such as 'x_dot = y + 3'

- Simulation
  - [x] Using C2E2 C++ Backend to do Simualtion
  - [ ] Instead of wrapping all information into a Model(), maybe binding more C++ classes/functions for the python side to use.
  - [ ] Currently output is stored in a file, to speedup the program, maybe return the result directly back to Python.

- Verification
  - [ ] Not yet started.

# Install
```zsh
# Clone this repo.
git clone https://github.com/bznick98/pyC2E2.git
# Go to the project directory
cd pyC2E2
# Install the module into current directory
pip install -e .
# Install all C++ Dependent Libraries (brew is for Mac)
brew install ppl
brew install eigen
brew install glpk
# Make C++ into a library
make
```

# Use Cases
## Loading from a .hyxml file
```python
from c2e2 import *

# Using FileHandler to load a hyxml file (file path may change depending where are you executing the python shell)
automata, properties = FileHandler.load_model("sample_hyxml/TotalMotionV2.hyxml")

# All the automata will be loaded into 'automata'
# All the properties will be loaded in to 'properties'
# You can display info of an automaton by directly printing automaton's name
print(automata)
print(properties)
```

## Operations of each Data Type:
- RectangleSet()
```python
from c2e2 import *

# Empty Set
R1 = RectangleSet()

# Creating Initial Set from string expressions
R2 = RectangleSet("Mode1: x>=1 && x<=3 && y>=-1 && y<=0")
# or Unsafe Set (without specifying mode name)
R3 = RectangleSet("x>10")

# Assign/Update new expressions to the Set
R3("x>5")
# or more explicity R3.set_expression("x>5")

```

- DAI() (Deterministic algebraic inequalities)
```python
from c2e2 import *

# Empty DAI Equation
D1 = DAI()

# Creating DAI from string expressions
D2 = DAI("x = 3 + y^2")
D3 = DAI("x = y**2 / z")

# Assign/Update new expressions to the DAI
D3("x = y + 5")
# or more explicity D3.set_expression("x = y + 5")

```

# Full Example
- See `tests/linThermo.py`, it's an example constructed manually, produce same result as loading linThermo.hyxml then do simulation.

# Some Urgent TODOs
- [ ] Implement Unsafe Set Bound Checkings! (Avoid like x>5 && x<3).

# C++ lib Requirement:
```zsh
brew install ppl
brew install eigen
brew install glpk
```