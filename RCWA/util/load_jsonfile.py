import os.path
import json

def load_jsonfile(jsonfile, **conf):
    
    assert os.path.isfile(jsonfile), "%s is not a file."%jsonfile
    
    try:
        # load JSON script file
        res = json.loads(open(jsonfile).read())
        # use conf keyword arguments to override JSON script
        res.update(conf)
    except IOError:
        raise IOError("'%s' is not a file."%jsonfile)
    except ValueError:
        raise ValueError("Error: script file is formatted incorrectly.")
    
    return res