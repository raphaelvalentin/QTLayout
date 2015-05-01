class Exception(Exception):
    """Base class for exceptions in this module."""
    def __init__(self,msg=None):
        self.msg = msg
    def __str__(self):
        if self.msg:
            return self.msg
        else:
            return ''

# retro compatibility            
Error = Exception            
	    
# raise Exception('mon error')

import warnings as Warnings

# __warnings__ is one of the following ("ignore", "error", "always", "default"...)
try:
    Warnings.simplefilter(__warnings__, Warning)
except:
    Warnings.simplefilter("always", Warning)
   

def formatwarning(message, category, filename, lineno, file=None, line=None):
    return '%s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
Warnings.formatwarning = formatwarning

def warning(test, msg):
    if test:
        Warnings.warn(msg, Warning, stacklevel=2)

################
import sys

class Warning(object):
    Warnings = []
    def __init__(self, warning=''):
        if not warning in Warning.Warnings:
	    Warning.Warnings.append(warning)
	    sys.stderr.write(warning)
	    sys.stderr.flush()

def deprecated(cls):
    old_init = getattr(cls, '__init__')
    # overridding
    def __init__(self, *args, **kwargs):
        Warning('Warning: The class {name} is deprecated.\n'.format(name=cls.__name__))	
        old_init(self, *args, **kwargs)
    old_new = getattr(cls, '__new__')
    # overridding
    def new(cls, *args, **kwargs):
        Warning('Warning: The class {name} is deprecated.\n'.format(name=cls.__name__))
        if 'iteritems' in dir(cls):
            inst = dict.__new__(cls)
        elif 'append' in dir(cls):
            inst = list.__new__(cls)
        else:
            inst = super(cls.__class__, cls).__new__(cls)
        return inst
    cls.__new__ = staticmethod(new)
    setattr(cls, '__init__', __init__)
    return cls
