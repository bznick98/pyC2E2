# FileHandler.py Class List:
# FileHandler()

import os
import xml.dom.minidom
import xml.etree.ElementTree as ET

from .Automaton import Automaton, Mode, DAI, Invariant, Transition, Guard, Action, Variable, ThinVariable
from .Property import Property
from .Set import RectangleSet
from .SymEq import SymEq
from .constants import *

class FileHandler:
    '''
    This class is used to do basic file operations,
    such as loading a model from hyxml file, etc.
    '''
    
    @staticmethod
    def load_model(file_path):
        '''
        Loading Hybrid Model from the filePath given by the argument.

        Args:
            file_path (str): the path of the model file

        Returns:
            A list of Automaton()
            A list of Property()
            or False if failed to read the file
        '''

        print("Loading file from \"{}\"...".format(file_path))
        # get the file name ("raw_name.ext" = base_name)
        base_name = os.path.basename(file_path)
        raw_name, ext = os.path.splitext(base_name)

        # parsing .hyxml into a Hyxml Tree
        if ext == '.hyxml':
            hyxml_tree = ET.parse(file_path)
            if (hyxml_tree == None):
                return False
            
            hyxml_root = hyxml_tree.getroot()
            hyxml_type = hyxml_root.get('type')

            thinvarprop = ""
            thinvarlist = ""

            if(hyxml_type == 'Model'):
                # read automata and properties respectively
                automata = FileHandler.open_hyxml_model(hyxml_root)
                properties = FileHandler.open_hyxml_properties(hyxml_root)

            else:
                raise Exception("Hyxml type is not 'Model'")
        else:
            raise Exception("Currently only support hyxml")
        
        return automata, properties

    @staticmethod
    def open_hyxml_model(root):
        """
        Loads automata from input xml root

        Args:
            root: XML root

        Returns:
            List of Automaton() objects
        """

        print("Loading hyxml model...", end='')

        automata = []

        for auto in root.iterfind("automaton"):

            name = auto.get("name")
            automaton = Automaton(name)
            print("[NEW ROUND!!!!]")
            print(automaton)
            
            # Load variables
            for var in auto.iterfind("variable"):
                v_name = var.get("name")
                v_scope = var.get("scope")
                v_type = var.get("type")
               
                v = Variable(name=v_name, type=v_type, scope=v_scope)
                # Variable CHECK
                if v not in automaton.variable_list:
                    automaton.add_variable(v)
                else:
                    raise Exception("[AUTOMATON LOAD ERROR]: Variable " + v.name + " Already Exists!")

            # Load thin variables (not sure this, another type of variable?)
            for thinvar in auto.iterfind("thin_variable"):
                v_name = thinvar.get("name")
                v_scope = thinvar.get("scope")
                v_type = thinvar.get("type")
                
                v = ThinVariable(name=v_name, type=v_type, scope=v_scope)
                # ThinVariable CHECK
                if v not in automaton.variable_list:
                    automaton.add_thinvar(v)
                else:
                    raise Exception("[AUTOMATON LOAD ERROR]: ThinVariable " + v.name + " Already Exists!")

            # Load modes
            for mode in auto.iterfind("mode"):
                mode_name = mode.get("name")
                mode_id = int(mode.get("id"))
                mode_init = (mode.get("initial") == "True")
                
                mode_obj = Mode(name=mode_name, id=mode_id, isInitial=mode_init)

                # for each mode, load DAI equations
                for dai in mode.iterfind("dai"):
                    raw_eq = dai.get("equation")                    
                    mode_obj.add_dai(DAI(raw_eq))

                # for each mode, load Invariant
                for inv in mode.iterfind("invariant"):
                    raw_eq = inv.get("equation")
                    # Equation 'cleaning' is needed for inequalities
                    clean_eq = FileHandler.clean_eq(raw_eq)
                    mode_obj.add_invariant(Invariant(clean_eq))
                
                automaton.add_mode(mode_obj, mode_id)
            
            # MODE ID/Name CHECK
            print("BEFORE VERIFYING MODEID!!")
            if not automaton.verify_mode_ids():
                raise Exception("[AUTOMATON LOAD ERROR]: {}'s MODE IDs NOT UNIQUE!\n".format(automaton.name))
            if not automaton.verify_mode_names():
                raise Exception("[AUTOMATON LOAD ERROR]: {}'s MODE NAMES NOT UNIQUE!\n".format(automaton.name))

            # Load transitions
            for tran in auto.iterfind("transition"):
                g = tran.find("guard")
                guard = Guard(FileHandler.clean_eq(g.get("equation")))

                tran_id = int(tran.get("id"))
                tran_src = int(tran.get("source"))
                tran_dest = int(tran.get("destination"))
                
                # Actions
                actions = []
                for act in tran.iterfind("action"):
                    # raw equation
                    raw_eq = act.get("equation")
                    clean_eq = FileHandler.clean_eq(raw_eq)
                    actions.append(Action(clean_eq))

                transition = Transition(guard, actions, tran_id, 
                                        tran_src, tran_dest)
                automaton.add_transition(transition)

            # Transition src/dest CHECK
            if not automaton.verify_transition_src_dest():
                raise Exception("[AUTOMATON LOAD ERROR]: Transition Source/Destination NOT valid!\n\
                    (FROM Automaton: " + automaton.name +")\n")

            automata.append(automaton)

        print("  Done, all automata are loaded. (+Validated)")
        return automata

    @staticmethod
    def open_hyxml_properties(root):
        """ 
        Load properties from hyxml
        
        Args:
            root: XML root

        Returns:
            List of Proerty() objects
        """
        
        print("Loading hyxml properties...", end='')

        prop_list = []
        for prop in root.iterfind('property'):
            # A new Property
            p = Property()
            p.name = prop.get('name')

            p.type = SAFETY 
            p.initial_set_raw = FileHandler.clean_eq(prop.get('initialSet'))
            p.unsafe_set_raw = FileHandler.clean_eq(prop.get('unsafeSet'))

            # Convert RAW Expressions to RectangleSet
            p.initial_set = RectangleSet(p.initial_set_raw)
            p.unsafe_set = RectangleSet(p.unsafe_set_raw)

            # Handle properties parameters
            param = prop.find('parameters')
            if param is not None:
                time_step = param.get('timestep')
                if time_step is None:
                    p.time_step = 0.0
                else:
                    p.time_step = float(time_step)

                time_horizon = param.get('timehorizon')
                if time_horizon is None:
                    p.time_horizon = 0.0
                else:
                    p.time_horizon = float(time_horizon)

                k_value = param.get('kvalue')
                if k_value is None:
                    p.k_value = 0.0
                else:
                    p.k_value = float(k_value)

            # Validate some of the attributes of the Property
            p.validate_initial()
            # Initial Set/Unsafe Set will be validated before doing Verification

            prop_list.append(p)

        print("  Done, all properties are loaded.\n")
        return prop_list

    @staticmethod   
    def clean_eq(eq):
        '''
        given a Raw equation, replace inequality notations
        '''
        r_dict = {'&lt;':'<', '&gt;':'>', '&amp;':'&', ' and ':'&&', ' or ':'||'}
        for term in r_dict:
            eq = eq.replace(term, r_dict[term])
        return eq