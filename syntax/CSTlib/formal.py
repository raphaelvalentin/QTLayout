import math
from exceptions import Exception

class formal(str):
    def __add__(self, obj):
        return formal('({a}+{b})'.format(a=self, b=obj))
    def __radd__(self, obj):
        return formal('({b}+{a})'.format(a=self, b=obj))
    def __mul__(self, obj):
        return formal('({a}*{b})'.format(a=self, b=obj))
    def __rmul__(self, obj):
        return formal('({b}*{a})'.format(a=self, b=obj))
    def __sub__(self, obj):
        return formal('({a}-{b})'.format(a=self, b=obj))
    def __rsub__(self, obj):
        return formal('({b}-{a})'.format(a=self, b=obj))
    def __div__(self, obj):
        return formal('({a}/{b})'.format(a=self, b=obj))
    def __rdiv__(self, obj):
        return formal('({b}/{a})'.format(a=self, b=obj))
    def __pow__(self, obj):
        return formal('({a}**{b})'.format(a=self, b=obj))
    def __rpow__(self, obj):
        return formal('({b}**{a})'.format(a=self, b=obj))
    def __neg__(self):
        return formal('(-{a})'.format(a=self) )       
    def __pos__(self):
        return formal('(+{a})'.format(a=self) )
    def __ge__(self, obj):
        return formal('({a}>={b})'.format(a=self, b=obj) )
    def __le__(self, obj):
        return formal('({a}<={b})'.format(a=self, b=obj) )
    def __gt__(self, obj):
        return formal('({a}>{b})'.format(a=self, b=obj) )
    def __lt__(self, obj):
        return formal('({a}<{b})'.format(a=self, b=obj) )
    def __ne__(self, obj):
        return formal('({a}!={b})'.format(a=self, b=obj) )
    def __eq__(self, obj):
        return formal('({a}=={b})'.format(a=self, b=obj) )
    
    
def sqrt(obj):
    if isinstance(obj, formal):
        return formal('sqrt({a})'.format(a=obj))
    return math.sqrt(obj)

def cos(obj):
    if isinstance(obj, formal):
        return formal('cos({a})'.format(a=obj))
    return math.cos(obj)

def sin(obj):
    if isinstance(obj, formal):
        return formal('sin({a})'.format(a=obj))
    return math.sin(obj)

def tan(obj):
    if isinstance(obj, formal):
        return formal('tan({a})'.format(a=obj))
    return math.tan(obj)

def atan(obj):
    if isinstance(obj, formal):
        return formal('atan({a})'.format(a=obj))
    return math.atan(obj)

def IIf(condexpr, TruePart, FalsePart):
    if isinstance(condexpr, formal):
        return 'IIf({condexpr}, {TruePart}, {FalsePart})'.format(condexpr=condexpr, TruePart=TruePart, FalsePart=FalsePart)
    elif isinstance(condexpr, (int, float)):
        if condexpr:
            return TruePart
        else:
            return FalsePart
    else:
        raise Exception('Wrong type in IIF function')

__min__ = min
def min(*args):
    if len(args)==2:
        if isinstance(args[0], formal) or isinstance(args[1], formal):
            return IIf(args[0]<args[1], args[0], args[1])
    return __min__(*args)

__max__ = max
def max(*args):
    if len(args)==2:
        if isinstance(args[0], formal) or isinstance(args[1], formal):
            return IIf(args[0]<args[1], args[1], args[0])
    return __max__(*args)

#x = formal('x')
#y = formal('y')
#print (x+4.56*2)/2.-y



