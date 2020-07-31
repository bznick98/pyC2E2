# c2e2lib-dev
## Functionalities & Todos
- Overall:
  - [x] Be able to load model from .hyxml file.
  - [x] Be able to load model from DryVR Black Box Model.
  - [ ] Support reading DryVR inputFile (mode transition information).
    - Question: It seems like in DryVR's inputFile, mode transition is based on time instead of guards, is it how DryVR is designed?
  - [ ] File Path in the complete.py needs to be revised to support cross-platform

- Automaton():
  - [ ] ...

- RectangleSet()
  - [x] RectangleSet(): Can be used to describe Initial Set & Unsafe Set.
  - [ ] Is it better to have InitialSet/UnsafeSet class inherit from RectangleSet() so that we can declare those sets more explicitly?
  
- DAI Equations()
  - [x] DAI(): Can be used to describe differential equations such as 'x_dot = y + 3'

- Simulation
  - [x] Using C2E2 C++ Backend to do Simualtion
  - [x] Read output from file into a 2d list 
  - [ ] Instead of wrapping all information into a Model(), maybe binding more C++ classes/functions for the python side to use.
  - [ ] Currently output is stored in a file, to speedup the program, maybe return the result directly back to Python.

- Verification
  - [ ] Not yet started.

## Install
```zsh
# Clone this repo.
git clone https://github.com/C2E2-Development-Team/c2e2lib-dev.git
# Go to the project directory
cd c2e2lib-dev
# Install in the current directory (using pip3 if pip is python2)
pip install -e .
# Install all C++ Dependent Libraries (brew is for Mac)
brew install ppl
brew install eigen
brew install glpk
# Make C++ into a library
# TODO: Some path in the Makefile are somehow hardcoded to work on my local computer.
make
```

## Use Cases
### Loading from a .hyxml file
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

### Operations of each Data Type:
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

## Full Example
### Construct by hand
- See `use_cases/linThermo_manual.py`, it's an example constructed manually, produce same result as loading linThermo.hyxml then do simulation.
- Run linThermo example by running: (In order to run all )
  ```
  cd use_cases
  python3 linThermo_test.py
  ```
### Load from hyxml or DryVR file
- See `use_cases/hyxml_examples` or `use_cases/dryVR_examples`, examples will load the model from files and perform simulate.
- Run load model examples  
  ```
  cd use_cases/dryVR_examples/Thermostats
  python3 dryVR_Thermo_test.py
  ```
  or 
  ```
  cd use_cases/hyxml_examples/linThermo
  python3 linThermo_test.py
  ```

## Some Urgent TODOs
- [ ] Add Unsafe Set Boundary Checkings. (Avoid like x>5 && x<3).
- [ ] I used string concatenate instead of os.path.join() to join file paths in source code, so it's likely won't work in Windows. In order to make it cross platform, file path problems need to be addressed.
- [ ] Makefile is a problem, All of my C++ lib is installed by brew, so lots of the INCLUDE path points to `usr/local/Cellar/...`, which might not be the case for all users. Maybe use CMake to build as it's more intended for cross platform apps?

## C++ lib Requirement:
```zsh
ppl
eigen
glpk
```