# load all libraries in globals
from math import *
from syntax import *

class Parameter(float):
    list = list()
    def __new__(cls, name, value, comment=''):
        Parameter.list.append((name, value, comment))
        if name in globals():
            value = globals()[name]
        return float.__new__(cls, value)
#

__globals__ = globals()
__all__= ['RunScript']

class RunScript(object):
    def __init__(self, script, context={}): 
        self.script = script
        self.globals = {k:v for k,v in __globals__.iteritems()}
        self.globals.update(context)
        self.globals.update({"__file__": "<script>", "__name__": "__main__"})
        self.traceback = ''

    def run(self):
        import sys, traceback
        try:
            # execute the script
            code = compile(self.script, '', 'exec')
            exec(code, self.globals)
        except:
            traceback_lines = traceback.format_exc().split('\n')
            # Remove traceback mentioning this file, and a linebreak
            if __name__<>'__main__':
                for i in (1,1):
                    traceback_lines.pop(i)
            self.traceback = "\n".join(traceback_lines)
            sys.stderr.write(self.traceback)
        return self


if __name__=='__main__':
    script = """
def f():
    obj = Primitive()
    obj.append( Point(0,0,0) )
    obj.append( Point(0,1,0) )
    obj.append( Point(1,1,0) )
    obj.append( Point(1,0,0) )
    return obj

obj = f()
print obj
"""
    filename = "D:/Work/dev/python/layerview/test_inductor2.py"
    script = open(filename).read()
    obj = RunScript(script).run()
