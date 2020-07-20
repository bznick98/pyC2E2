
# Current Functionalities
- [x] Be able to load model from .hyxml file.
- [x] Some Set Operations.
- [ ] Using C2E2 C++ Backend to do Simualtion/Verification

# Install
```zsh
# Clone this repo.
git clone https://github.com/bznick98/pyC2E2.git
# Go to the project directory
cd pyC2E2
# install into global environment
pip install .
```

# Use Case
- Loading a .hyxml
```python
from pyC2E2 import *

# Using FileHandler to load a hyxml file
automata, properties = FileHandler.load_model("TotalMotionV2.hyxml")

# All the automata will be loaded into 'automata'
# All the properties will be loaded in to 'properties'
# You can display info of an automaton by directly printing automaton's name
print(automata)
print(properties)
```

- Operations related to RectangleSet()
```python
from pyC2E2 import *

# Empty Set
R1 = RectangleSet()

# Creating Initial Set from string expressions
R2 = RectangleSet("Mode1: x>=1 && x<=3 && y>=-1 && y<=0")
# or Unsafe Set
R3 = RectangleSet("x>10")

# Assign/Update new expressions to the Set
R3.set_expression("x>5")
```

# Some Urgent TODOs
- [ ] Implement Unsafe Set Bound Checkings! (Avoid like x>5 && x<3).