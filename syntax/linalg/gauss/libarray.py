import math
import pprint

__all__ = ['DTypeFromList', 'abs', 'acos', 'array', 'asin', 'atan', 'atan2', 'copy', 
           'cos', 'cosh', 'diag', 'dot', 'exp', 'fabs', 'flatten', 'identity', 'inf', 
	   'inline', 'log', 'log10', 'matrix', 'nan', 'shape', 'sin', 'sinh', 
	   'sqrt', 'tan', 'transpose', 'zeros']
	   
nan = float('nan')
inf = float('inf')

class array(list):
    def __init__(self, object, dtype=None, _parse=True):
        try:
	    if _parse:
                for i, element in enumerate(object):
                    if isinstance(element, (list,tuple)):
                        object[i] = array(element, dtype)
                    elif dtype:
                        object[i] = dtype(element)
            list.__init__(self, object)
        except:
            raise TypeError, 'Setting an array element with a sequence'
 
    @property
    def imag(self):
        return array([element.imag for element in self], _parse=False)
 
    @property
    def real(self):
        return array([element.real for element in self], _parse=False)

    @property
    def shape(self):
        return shape(self)
    
    def copy(self):
        return copy(self)
 
    def __add__(self, object):
        if isinstance(object, (int, float, complex)):
            return array([element+object for element in self], _parse=False)
        elif isinstance(object, list):
            if len(self)==len(object):
                return array([element1+element2 for element1, element2 in zip(self, object)], _parse=False)
            else:
                raise ValueError, 'operands could not be broadcast together'
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __sub__(self, object):
        if isinstance(object, (int, float, complex)):
            return array([element-object for element in self], _parse=False)
        elif isinstance(object, list):
            if len(self)==len(object):
                return array([element1-element2 for element1, element2 in zip(self, object)], _parse=False)
            else:
                raise ValueError, 'operands could not be broadcast together'
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __mul__(self, object):
        if isinstance(object, (int, float, complex)):
            return array([element*object for element in self], _parse=False)
        elif isinstance(object, list):
            if len(self)==len(object):
                return array([element1*element2 for element1, element2 in zip(self, object)], _parse=False)
            else:
                raise ValueError, 'operands could not be broadcast together'
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __div__(self, object):
        if isinstance(object, (int, float, complex)):
            try:
                return array([element/object for element in self], _parse=False)
            except ZeroDivisionError:
	        print 'Warning: a ZeroDivisionError exception occurs'
                return array([nan for i in xrange(len(self))], _parse=False)
        elif isinstance(object, list):
            if len(object)==len(self):
                result = [0 for i in xrange(len(self))]
                for i in xrange(len(self)):
                    try:
                        result[i] = self[i] / object[i]
                    except ZeroDivisionError:
		        print 'Warning: a ZeroDivisionError exception occurs'
                        result[i] = nan
            else:
                raise ValueError, 'operands could not be broadcast together'
            return array(result, _parse=False)
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __rdiv__(self, object):
        if isinstance(object, (int, float, complex)):
            try:
                return array([object/element for element in self], _parse=False)
            except ZeroDivisionError:
	        print 'Warning: a ZeroDivisionError exception occurs'
                return array([nan for i in xrange(len(self))], _parse=False)
        elif isinstance(object, list):
            if len(object)==len(self):
                result = [None for i in xrange(len(self))]
                for i in xrange(len(self)):
                    try:
                        result[i] = object[i]/self[i]
                    except ZeroDivisionError:
                        result[i] = nan
            else:
                raise ValueError, 'operands could not be broadcast together'
            return array(result, _parse=False)
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __pow__(self, object):
        if isinstance(object, (int, float, complex)):
            try:
                return array([element**object for element in self], _parse=False)
            except ZeroDivisionError:
	        print 'Warning: a ZeroDivisionError exception occurs'
                return array([nan for i in xrange(len(self))], _parse=False)
        elif isinstance(object, list):
            if len(object)==len(self):
                result = [None for i in xrange(len(self))]
                for i in xrange(len(self)):
                    try:
                        result[i] = self[i] ** object[i]
                    except ZeroDivisionError:
		        print 'Warning: a ZeroDivisionError exception occurs'
                        result[i] = nan
            else:
                raise ValueError, 'operands could not be broadcast together'
            return array(result, _parse=False)
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __rpow__(self, object):
        if isinstance(object, (int, float, complex)):
            try:
                return array([object**element for element in self], _parse=False)
            except ZeroDivisionError:
	        print 'Warning: a ZeroDivisionError exception occurs'
                return array([nan for i in xrange(len(self))])
        elif isinstance(object, list):
            if len(object)==len(self):
                result = [None for i in xrange(len(self))]
                for i in xrange(len(self)):
                    try:
                        result[i] = object[i]**self[i]
                    except ZeroDivisionError:
		        print 'Warning: a ZeroDivisionError exception occurs'
                        result[i] = nan
            else:
                raise ValueError, 'operands could not be broadcast together'
            return array(result, _parse=False)
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
            return array([object-element for element in self], _parse=False)
        elif isinstance(object, list):
            if len(self)==len(object):
                return array([element2-element1 for element1, element2 in zip(self, object)], _parse=False)
            else:
                raise ValueError, 'operands could not be broadcast together'
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __rmul__(self, object):
        if isinstance(object, (int, float, complex)):
            return array([element*object for element in self], _parse=False)
        elif isinstance(object, list):
            if len(self)==len(object):
                return array([element1*element2 for element1, element2 in zip(self, object)], _parse=False)
            else:
                raise ValueError, 'operands could not be broadcast together'
        else:
            raise TypeError, 'operands could not be broadcast together'
 
    def __neg__(self):
        return array([-element for element in self], _parse=False)
 
    def __pos__(self):
        return self
 
    @property
    def T(self):
        return transpose(self)
   
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

    def flatten(self):
        return array(flatten(self), _parse=False)

    def __str__(self):
        s = pprint.pformat(self.tolist()).split('\n')
	r = ['array(' + s[0]]
	for line in s[1:]:
	    r.append('      %s'%line)
	r[-1] = r[-1] + ')'
        return "\n".join(r)

    def ___repr__(self):
         return self.__str__()        
    

def flatten(sequence):
    """ yield each element of an irregular list (or tuple, dict...->__instance)"""
    for el in sequence:
        if isinstance(el, list):
            for sub in flatten(el):
                yield sub
        else:
            yield el
    
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
	    l = len(obj[0])
	    if False in [l==len(obji) for obji in obj[1:]]:
	        return tuple(_shape)
            _shape.extend(shape(obj[0]))
    return tuple(_shape)

def copy(obj):
    if isinstance(obj, list):
	return type(obj)([ copy(x) for x in obj])
    else:
        return obj
       
def dot(arr1, arr2):
    m, n, p = len(arr1), len(arr2[0]), len(arr1[0])
    result = [[0.0 for j in xrange(n)] for i in xrange(m)]
    for i in xrange(m):
        for j in xrange(n):
            for k in xrange(p):
                result[i][j] += arr1[i][k]*arr2[k][j]
    return array(result)
       
#def dot(arr1, arr2):
#    tpos_b = zip(*arr2)
#    return array( [[ sum(ea*eb for ea, eb in zip(a,b)) for b in tpos_b] for a in arr1] )
 
def transpose(object):
    m, n = len(object), len(object[0])
    return array([[object[i][j] for i in xrange(m)] for j in xrange(n)])
 
def transpose(object):
    return array(zip(*object))
 
def diag(arr):
    return array([arr[i][i] for i in xrange(len(arr))])
 
def zeros(shape):
     if len(shape)==1:
         return array([0. for i in xrange(shape[0])])
     else:
         return array([zeros(shape[1:]) for i in xrange(shape[0])])

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
            return nan
 
def exp(obj):
    if isinstance(obj, list):
        return array([exp(subobj) for subobj in obj])
    elif obj<709.7:
        return math.exp(obj)
    return nan
 
def log(obj):
    if isinstance(obj, list):
        return array([log(subobj) for subobj in obj])
    else:
        try:
            return math.log(obj)
        except:
            return nan
           
def log10(obj):
    if isinstance(obj, list):
        return array([log10(subobj) for subobj in obj])
    else:
        try:
            return math.log10(obj)
        except:
            return nan
 
def sqrt(obj):
    if isinstance(obj, list):
        return array([sqrt(subobj) for subobj in obj])
    else:
        try:
            return math.sqrt(obj)
        except:
            return nan
 
def atan(obj):
    if isinstance(obj, list):
        return array([atan(subobj) for subobj in obj])
    else:
        try:
            return math.atan(obj)
        except:
            return nan
       
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
 
def atan2(obj1, obj2):
    if isinstance(obj1, list) and isinstance(obj2, list) :
        return array([atan2(subobj1, subobj2) for subobj1, subobj2  in zip(obj1, obj2)])
    else:
        return math.atan2(obj1, obj2)
 
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
 
def matrix(*X):
    size = int(sqrt(len(X)))
    return [[X[size*j+i] for i in xrange(size)] for j in xrange(size)]

def inline(X):
    m, n = len(X), len(X[0])
    return [ X[i][j] for i in xrange(m) for j in xrange(n)]



if __name__ == '__main__':
    import unittest
    import random
    # unittest
    class TestCase(unittest.TestCase):
        def test1(self):
	    a = array([1,2,3])
	    b = array([1,2,3])
	    assert a==b
        def test2(self):
	    a = array([1,2,3])
	    b = a.copy()
	    assert a==b
        def test3(self):
	    a = array([1,2,3])
	    b = a.copy()
	    b[0] = 123
	    assert a<>b
        def test4(self):
	    a = array([random.random() for i in xrange(5)])
	    b = array([random.random() for i in xrange(5)])
	    c1 = array([x+y for x, y in zip(a, b)])
	    c2 = a + b
	    assert c1==c2
        def test5(self):
	    a = array([random.random() for i in xrange(5)])
	    b = array([random.random() for i in xrange(5)])
	    c1 = array([x-y for x, y in zip(a, b)])
	    c2 = a - b
	    assert c1==c2
        def test6(self):
	    a = array([random.random() for i in xrange(5)])
	    b = array([random.random() for i in xrange(5)])
	    c1 = array([x*y for x, y in zip(a, b)])
	    c2 = a * b
	    assert c1==c2
        def test7(self):
	    a = array([random.random() for i in xrange(5)])
	    b = array([random.random() for i in xrange(5)])
	    c1 = array([x/y for x, y in zip(a, b)])
	    c2 = a / b
	    assert c1==c2
        def test8(self):
	    a = array([random.random() for i in xrange(5)])
	    b = array([random.random() for i in xrange(5)])
	    c1 = array([x**y for x, y in zip(a, b)])
	    c2 = a ** b
	    assert c1==c2
        def test9(self):
	    a = array([[1,2],[3, 4]])
	    b = a.copy()
	    assert a==b
        def test10(self):
	    a = array([[1,2],[3, 4]])
	    b = a.copy()
	    b[0][0] = 123
	    assert a<>b
        def test11(self):
	    a = array([1,2,3, 4])
	    s = shape(a)
	    assert s == (4,)
        def test12(self):
	    a = array([[1,2],[3, 4]])
	    s = shape(a)
	    assert s == (2,2)
        def test13(self):
	    a = array([[1,2],[3, 4]])
	    t = a.T
	    assert t == array([[1,3],[2, 4]])
        def test14(self):
	    a = array([[random.random(),random.random()] for i in xrange(5)])
	    b = array([[random.random(),random.random()] for i in xrange(5)])
	    c1 = array([[x[0]+y[0], x[1]+y[1]] for x, y in zip(a, b)])
	    c2 = a + b
	    assert c1==c2
        def test15(self):
	    a = array([[random.random(),random.random()] for i in xrange(5)])
	    b = random.random()
	    c1 = array([[x[0]+b, x[1]+b] for x in a])
	    c2 = a + b
	    assert c1==c2
        def test16(self):
	    a = array([[random.random(),random.random()] for i in xrange(5)])
	    b = random.random()
	    c1 = array([[x[0]+b, x[1]+b] for x in a])
	    c2 = b + a
	    assert c1==c2
        def test17(self):
	    a = array([[random.random(),random.random()] for i in xrange(5)])
	    b = random.random()
	    c1 = array([[x[0]-b, x[1]-b] for x in a])
	    c2 = a - b
	    assert c1==c2
        def test18(self):
	    a = array([[random.random(),random.random()] for i in xrange(5)])
	    b = random.random()
	    c1 = array([[b-x[0], b-x[1]] for x in a])
	    c2 = b - a
	    assert c1==c2
        def test19(self):
	    a = array([[random.random(),random.random()] for i in xrange(5)])
	    b = random.random()
	    c1 = array([[x[0]/b, x[1]/b] for x in a])
	    c2 = a / b
	    assert c1==c2
        def test20(self):
	    a = array([[random.random(),random.random()] for i in xrange(5)])
	    b = random.random()
	    c1 = array([[b/x[0], b/x[1]] for x in a])
	    c2 = b / a
	    assert c1==c2
        def test21(self):
	    a = array([[1,2,3],[4,5,6],[7,8,9]])
	    c1 = array([[66,78,90],[78,93,108], [90,108,126]])
	    c2 = a.T.dot(a)
	    assert c1==c2
        def test22(self):
	    a = array([[1,2,3],[4,5,6],[7,8,9]])
	    c1 = array([1, 5, 9])
	    c2 = diag(a)
	    assert c1==c2
        def test23(self):
	    a = array([[1,2,3],[4,5,6],[7,8,9]])
	    c1 = array([[1,0,0],[0,5,0],[0,0,9]])
	    c2 = diag(a) * identity(3)
	    assert c1==c2
        def test24(self):
	    a = array([[1,2,3],[4,5,6],[7,8,9]])
	    c1 = array([[-1,-2,-3],[-4,-5,-6],[-7,-8,-9]])
	    c2 = -a
	    assert c1==c2
        def test25(self):
	    a = array([[random.random() for i in xrange(5)] for j in xrange(5)])
	    b = array([[random.random() for i in xrange(5)] for j in xrange(5)])
	    c = array([[random.random() for i in xrange(5)] for j in xrange(5)])
	    d = random.random()
	    c1 = array([[(x+y/z)*d for x, y, z in zip(rowx, rowy, rowz)] for rowx, rowy, rowz in zip(a, b, c)])
	    c2 = (a+b/c)*d
	    assert c1==c2
        def test26(self):
	    x = zeros((2,2,2))	    
	    x[0][0][0] = 123.
	    y = ([[[ 123.,  0.], [ 0.,  0.]], [[ 0.,  0.], [ 0.,  0.]]])
	    assert x==y

    unittest.main()

array([[[ 0.,  0.], [ 0.,  0.]], [[ 0.,  0.], [ 0.,  0.]]])
