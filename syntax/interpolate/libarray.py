import math
    
class array(list):
    def __init__(self, object, dtype=None):
        try:
            for i, element in enumerate(object):
                if isinstance(element, list):
                    object[i] = array(element, dtype)
                elif dtype:
                    object[i] = dtype(element)
            list.__init__(self, object)
        except:
            raise ValueError, 'Setting an array element with a sequence'
 
    def imag(self):
        return array([element.imag for element in self])
    imag = property(fget=imag)
 
    def real(self):
        return array([element.real for element in self])
    real = property(fget=real)

    def shape(self):
        return shape(self)
    shape = property(fget=shape)
    
    def copy(self):
        return array(self)
 
    def __add__(self, object):
        if isinstance(object, (int, float, complex)):
            return array([element+object for element in self])
        elif isinstance(object, list):
            if len(self)==len(object):
                return array([element1+element2 for element1, element2 in zip(self, object)])
            else:
                raise ValueError, 'operands could not be broadcast together'
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __sub__(self, object):
        if isinstance(object, (int, float, complex)):
            return array([element-object for element in self])
        elif isinstance(object, list):
            if len(self)==len(object):
                return array([element1-element2 for element1, element2 in zip(self, object)])
            else:
                raise ValueError, 'operands could not be broadcast together'
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __mul__(self, object):
        if isinstance(object, (int, float, complex)):
            return array([element*object for element in self])
        elif isinstance(object, list):
            if len(self)==len(object):
                return array([element1*element2 for element1, element2 in zip(self, object)])
            else:
                raise ValueError, 'operands could not be broadcast together'
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __div__(self, object):
        if isinstance(object, (int, float, complex)):
            try:
                return array([element/object for element in self])
            except ZeroDivisionError:
                return array([float('nan') for i in xrange(len(self))])
        elif isinstance(object, list):
            if len(object)==len(self):
                result = [None for i in xrange(len(self))]
                for i in xrange(len(self)):
                    try:
                        result[i] = self[i] / object[i]
                    except ZeroDivisionError:
                        result[i] = float('nan')
            else:
                raise ValueError, 'operands could not be broadcast together'
            return array(result)
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __rdiv__(self, object):
        if isinstance(object, (int, float, complex)):
            try:
                return array([object/element for element in self])
            except ZeroDivisionError:
                return array([float('nan') for i in xrange(len(self))])
        elif isinstance(object, list):
            if len(object)==len(self):
                result = [None for i in xrange(len(self))]
                for i in xrange(len(self)):
                    try:
                        result[i] = object[i]/self[i]
                    except ZeroDivisionError:
                        result[i] = float('nan')
            else:
                raise ValueError, 'operands could not be broadcast together'
            return array(result)
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __pow__(self, object):
        if isinstance(object, (int, float, complex)):
            try:
                return array([element**object for element in self])
            except ZeroDivisionError:
                return array([float('nan') for i in xrange(len(self))])
        elif isinstance(object, list):
            if len(object)==len(self):
                result = [None for i in xrange(len(self))]
                for i in xrange(len(self)):
                    try:
                        result[i] = self[i] ** object[i]
                    except ZeroDivisionError:
                        result[i] = float('nan')
            else:
                raise ValueError, 'operands could not be broadcast together'
            return array(result)
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __rpow__(self, object):
        if isinstance(object, (int, float, complex)):
            try:
                return array([object**element for element in self])
            except ZeroDivisionError:
                return array([float('nan') for i in xrange(len(self))])
        elif isinstance(object, list):
            if len(object)==len(self):
                result = [None for i in xrange(len(self))]
                for i in xrange(len(self)):
                    try:
                        result[i] = object[i]**self[i]
                    except ZeroDivisionError:
                        result[i] = float('nan')
            else:
                raise ValueError, 'operands could not be broadcast together'
            return array(result)
        else:
            raise TypeError, 'operands could not be broadcast together'
                   
    def __radd__(self, object):
        return self.__add__(object)
                   
    def __iadd__(self, object):
        return self.__add__(object)
                   
    def __isub__(self, object):
        return self.__sub__(object)
 
    def __rsub__(self, object):
        if isinstance(object, (int, float, complex)):
            return array([object-element for element in self])
        elif isinstance(object, list):
            if len(self)==len(object):
                return array([element2-element1 for element1, element2 in zip(self, object)])
            else:
                raise ValueError, 'operands could not be broadcast together'
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __rmul__(self, object):
        return self.__mul__(object)
 
    def __neg__(self):
        return array([-element for element in self])
 
    def __pos__(self):
        return self
 
    def T(self):
        return transpose(self)
    T = property(fget=T)
   
    def dot(self, obj):
        return dot(self, obj)
       
    def tolist(self):
        l = list(self)
        for i, element in enumerate(self):
            if isinstance(element, array):
               l[i] = element.tolist() 
        return l
 
    def __repr__ (self):
        return "".join(['array('+ str(self.tolist()) +')'])
            
    def astype(self, dtype):
        return array(_astype(self, dtype), dtype=dtype)
            
    def fill(self, value):
        list.__init__(array(_fill(self, value)))
        
def _fill(object, value):
    for i, element in enumerate(object):
        if isinstance(element, list):
            object[i] = _fill(element, value)
        else:
            object[i] = value
    return object

def _astype(object, dtype):
    if isinstance(object, list):
        return [_astype(element, dtype) for element in object]
    else:
        return dtype(object)
        
def shape(obj):
    _shape = []
    if isinstance(obj, list):
        _shape.append(len(obj))
        if isinstance(obj[0], list):
            _shape.extend(shape(obj[0]))
    return tuple(_shape)
       
def dot(arr1, arr2):
    m, n, p = len(arr1), len(arr2[0]), len(arr1[0])
    result = [[0.0 for j in xrange(n)] for i in xrange(m)]
    for i in xrange(m):
        for j in xrange(n):
            for k in xrange(p):
                result[i][j] += arr1[i][k]*arr2[k][j]
    return array(result)
 
def transpose(object):
    try:
        object = list(object)
        return array([[object[i][j] for i in xrange(len(object))] for j in xrange(len(object[0]))])
    except (IndexError, TypeError):
        return array(object)
 
def diag(arr):
    try:
        result = [[0.0 for j in xrange(len(arr))] for i in xrange(len(arr[0]))]
        for i in xrange(len(arr)):
            result[i][i] = arr[i][i]
    except (IndexError, TypeError):
        return array(arr)
    return array(result)
 
def zeros(shape):
    if len(shape)==1:
        return array([0. for i in xrange(shape[0])])
    if len(shape)==2:
        return array([[0.0 for j in xrange(shape[1])] for i in xrange(shape[0])])
    if len(shape)==3:
        return array([[[0.0 for j in xrange(shape[2])] for i in xrange(shape[1])] for k in xrange(shape[0])])
 
def identity(n):
    result = [[0.0 for j in xrange(n)] for i in xrange(n)]
    for i in xrange(n):
        result[i][i] = 1.0
    return array(result)
 
def cos(obj):
    if isinstance(obj, list):
        return array([cos(subobj) for subobj in obj])
    else:
        return math.cos(obj)
 
def sin(obj):
    if isinstance(obj, list):
        return array([sin(subobj) for subobj in obj])
    else:
        return math.sin(obj)
 
def tan(obj):
    if isinstance(obj, list):
        return array([tan(subobj) for subobj in obj])
    else:
        try:
            return math.tan(obj)
        except:
            return float('nan')
 
def exp(obj):
    if isinstance(obj, list):
        return array([exp(subobj) for subobj in obj])
    elif obj<709.7:
        return math.exp(obj)
    return float('nan')
 
def log(obj):
    if isinstance(obj, list):
        return array([log(subobj) for subobj in obj])
    else:
        try:
            return math.log(obj)
        except:
            return float('nan')
           
def log10(obj):
    if isinstance(obj, list):
        return array([log10(subobj) for subobj in obj])
    else:
        try:
            return math.log10(obj)
        except:
            return float('nan')
 
def sqrt(obj):
    if isinstance(obj, list):
        return array([sqrt(subobj) for subobj in obj])
    else:
        try:
            return math.sqrt(obj)
        except:
            return float('nan')
 
def norm(M):
    err = 0.
    for i in M:
        if isinstance(i, list):
            for j in i:
                err+=j*j
        else :
            err+=i*i
    return math.sqrt(err)
 
def atan(obj):
    if isinstance(obj, list):
        return array([atan(subobj) for subobj in obj])
    else:
        try:
            return math.atan(obj)
        except:
            return float('nan')
       
def asin(obj):
    if isinstance(obj, list):
        return array([asin(subobj) for subobj in obj])
    else:
        return math.asin(obj)
 
def acos(obj):
    if isinstance(obj, list):
        return array([acos(subobj) for subobj in obj])
    else:
        return math.acos(obj)
 
def fabs(obj):
    if isinstance(obj, list):
        return array([fabs(subobj) for subobj in obj])
    else:
        return math.fabs(obj)
 
def sinh(obj):
    if isinstance(obj, list):
        return array([sinh(subobj) for subobj in obj])
    else:
        return math.sinh(obj)
 
def cosh(obj):
    if isinstance(obj, list):
        return array([cosh(subobj) for subobj in obj])
    else:
        return math.cosh(obj)
 
__abs__ = abs
def abs(obj):
    if isinstance(obj, list):
        return array([abs(subobj) for subobj in obj])
    else:
        return __abs__(obj)

def DTypeFromList(object):
    res = None
    for element in object:
        t = type(element)
        if isinstance(element, list):
            t = DTypeFromObject(element)
        if t==complex:
            res = complex
            break
        elif t==float and res<>complex:
            res = float
        elif t==int and res<>complex and res<>float:
            res = int
    return res
 
 
if __name__ == '__main__':
    x, y = 5., 3.
    a, b = array([5.]), array([3.])
    assert (a+b)[0]==x+y, ''
    assert (b+a)[0]==y+x, ''
    assert (a-b)[0]==x-y, ''
    assert (b-a)[0]==y-x, ''
    assert (a/b)[0]==x/y, ''
    assert (b/a)[0]==y/x, ''
    assert (a*b)[0]==x*y, ''
    assert (b*a)[0]==y*x, ''
 
    assert (x+b)[0]==x+y, ''
    assert (b+x)[0]==y+x, ''
    assert (x-b)[0]==x-y, ''
    assert (b-x)[0]==y-x, ''
    assert (x/b)[0]==x/y, ''
    assert (b/x)[0]==y/x, ''
    assert (x*b)[0]==x*y, ''
    assert (b*x)[0]==y*x, ''
    assert (b**x)[0]==y**x, ''
    assert (sin(b)*sin(x))[0]==math.sin(y)*math.sin(x), ''
 
    x, y = 5, 3
    a, b = array([5]), array([3])
   
    assert (a+b)[0]==x+y, ''
    assert (b+a)[0]==y+x, ''
    assert (a-b)[0]==x-y, ''
    assert (b-a)[0]==y-x, ''
    assert (a/b)[0]==x/y, ''
    assert (b/a)[0]==y/x, ''
    assert (a*b)[0]==x*y, ''
    assert (b*a)[0]==y*x, ''
 
    assert (x+b)[0]==x+y, ''
    assert (b+x)[0]==y+x, ''
    assert (x-b)[0]==x-y, ''
    assert (b-x)[0]==y-x, ''
    assert (x/b)[0]==x/y, ''
    assert (b/x)[0]==y/x, ''
    assert (x*b)[0]==x*y, ''
    assert (b*x)[0]==y*x, ''
