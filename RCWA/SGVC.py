import os.path
import inspect
from RCWA.util.load_jsonfile import load_jsonfile

class SGVC(object):
    """Subwavelength Grating Vortex Coronagraph (backbone) class
    
    This class object is responsible for instantiating all objects required 
    to carry out a SGVC simulation.
    
    Args:
        \*\*conf:
            user configuration values
        jsonfile (string):
            JSON script file. If not set, assumes that dictionary has been 
            passed through conf.
            
    Attributes:
        name (type):
            Description
    
    """
    
    def __init__(self, jsonfile=None, **conf):
        """Initializes all modules from a given script file or conf dictionary.
        
        Args: 
            jsonfile (string):
                Path to JSON script file. If not set, assumes that 
                dictionary has been passed through conf.
            conf (dictionary):
                Dictionary containing additional user configuration values and 
                desired module names.
            
        """
        # ensure JSON filename extension
        if jsonfile[-5:].lower() != '.json':
            jsonfile += '.json'
        # ensure JSON file location
        if not os.path.isfile(jsonfile):
            jsonfile = os.path.join(os.path.split(inspect.getfile( \
                    self.__class__))[0], 'script', jsonfile)
        # load JSON script file
        conf = load_jsonfile(jsonfile, **conf)
        
        # load attributes
        self.conf = conf
        
