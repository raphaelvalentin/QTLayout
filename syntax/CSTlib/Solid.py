__all__ = ['Solid', 'Solids', 'Subtract', 'Insert', 'Intersect']
try:
    from CSTlib import Brick, Extrude, Copy
except:
    pass
from exceptions import *
from formal import *

class Solid(list):
    _names = {'Solid':1}
    def __init__(self, *args, **kwargs):
        name = kwargs.get('name', 'Solid')
        if name in Solid._names:
            self.name = "{name}{i}".format(name=name, i=Solid._names[name])
            Solid._names[name] += 1
        else:
            self.name = "{name}".format(name=name)
            Solid._names[name] = 1
        self.component = kwargs.get('component', None)
        self.material = kwargs.get('material', None)
        if args:
            try:
                if not self.component:
                    self.component = args[0].component
                if not self.material:
                    self.material = args[0].material
            except:
                raise Exception('the first argument is not a Brick or a Solid')
        list.__init__(self, args)
                
    def __str__(self):
        ref = self[0]
        for elt in self:
            if isinstance(elt, (Brick, Solid, Extrude)):
                if elt.component <> ref.component:
                    print "Warning: %s.component is different than %s.component in Solid %s"%(elt.name, ref.name, self.name)
                if elt.material <> ref.material:
                    print "Warning: %s.material is different than %s.material in Solid %s"%(elt.name, ref.name, self.name)
        if len(self)==0:
            return ''
        elif len(self)==1:
            return '\n'.join(flatten([[ str(elt) for elt in self ],
                                      [ str(Rename(self[0], self.name)) ]
                                     ]))
        else:
            Volumes = []
            for v1 in self:
                if isinstance(v1, (Brick, Solid, Extrude)):
                    Volumes.append(v1)
                    for v2 in self:
                        if isinstance(v2, Subtract):
                            if v1 in v2[1:]:
                                Volumes.pop()
                                break
            if len(Volumes)<2:
                return '\n'.join(flatten([[ str( elt ) for elt in self ],
                                          [ str( Rename(self[0], self.name) ) ]
                                         ]))
                
            return '\n'.join(flatten([[ str( elt ) for elt in self ],
                                      [ str( Add(*Volumes) ) ],
                                      [ str( Rename(self[0], self.name) ) ]
                                     ]))

    @staticmethod
    def Add(*args):
        return Add(*args)

    @staticmethod
    def Subtract(*args):
        return Subtract(*args)

    @staticmethod
    def Insert(*args):
        return Insert(*args)

    @staticmethod
    def Intersect(*args):
        return Intersect(*args)

    @staticmethod
    def Rename(*args):
        return Rename(*args)

    def __copy__(self):
        solid = Solid()
        for obj in self:
            if '__copy__' in dir(obj):
                solid.append( obj.__copy__() )
            else:
                solid.append( Copy(obj) )
        return solid


class Solids(list):
    def __str__(self):
        return '\n'.join( str(solid) for solid in self)
            
        


def flatten(sequence):
    """ yield each element of an irregular list (or tuple, dict...->__instance)"""
    for el in sequence:
        if isinstance(el, (list, tuple)):
            for sub in flatten(el):
                yield sub
        else:
            yield el

        
class Add(list):
    def __init__(self, *solids):
        for solid in solids:
            if not isinstance(solid, (Brick, Solid, Extrude)):
                raise Exception('solids are not Brick or Solid')
        if len(solids)<2:
            raise Exception('number of solids is less than two')
        list.__init__(self, solids)
    def __str__(self):
        return '\n'.join(flatten([( 'Solid.Add "{solid1.component}:{solid1.name}", "{solid2.component}:{solid2.name}"'.format(solid1=self[0],
                                                                                                                              solid2=solid),
                                    '',
                                   ) for solid in self[1:]
                                 ]))


class Subtract(list):
    def __init__(self, *solids):
        for solid in solids:
            if not isinstance(solid, (Brick, Solid, Extrude)):
                raise Exception('solids are not Brick or Solid')
        if len(solids)<2:
            raise Exception('number of solids is less than two')
        list.__init__(self, solids)
    def __str__(self):
        return '\n'.join(flatten([( 'With Solid',
                                    '  .Version 9',
                                    '  .Subtract "{solid1.component}:{solid1.name}", "{solid2.component}:{solid2.name}"'.format(solid1=self[0],
                                                                                                                              solid2=solid),
                                    '  .Version 1',
                                    'End With',
                                    '',
                                   ) for solid in self[1:]
                                 ]))


class Insert(list):
    def __init__(self, *solids):
        #for solid in solids:
        #    if not isinstance(solid, (Brick, Solid, Extrude)):
        #        raise Exception('solids are not Brick or Solid')
        if len(solids)<2:
            raise Exception('number of solids is less than two')
        list.__init__(self, solids)
    def __str__(self):
        return '\n'.join(flatten([( 'With Solid',
                                    '  .Version 9',
                                    '  .Insert "{solid1.component}:{solid1.name}", "{solid2.component}:{solid2.name}"'.format(solid1=self[0],
                                                                                                                              solid2=solid),
                                    '  .Version 1',
                                    'End With',
                                    '',
                                   ) for solid in self[1:]
                                 ]))

class Intersect(list):
    def __init__(self, *solids):
        for solid in solids:
            if not isinstance(solid, (Brick, Solid, Extrude)):
                raise Exception('solids are not Brick or Solid')
        if len(solids)<2:
            raise Exception('number of solids is less than two')
        list.__init__(self, solids)
    def __str__(self):
        return '\n'.join(flatten([( 'With Solid',
                                    '  .Version 9',
                                    '  .Intersect "{solid1.component}:{solid1.name}", "{solid2.component}:{solid2.name}"'.format(solid1=self[0],
                                                                                                                              solid2=solid),
                                    '  .Version 1',
                                    'End With',
                                    '',
                                   ) for solid in self[1:]
                                 ]))


class Rename(object):
    def __init__(self, solid, name):
        if not isinstance(solid, (Brick, Solid, Extrude)):
            raise Exception('solid is not a Brick or a Solid')
        self.solid = solid
        self.name = name
    def __str__(self):
        return '\n'.join(flatten([( 'Solid.Rename "{self.solid.component}:{self.solid.name}", "{self.name}"'.format(self=self),
                                    '',
                                   )
                                 ]))



