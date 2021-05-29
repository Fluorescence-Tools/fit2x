"""Use a table to redirect the access of attributes and methods to class members


Example
-------

    class ClassA():

        @property
        def y(self):
            return self._y

        @y.setter
        def y(self, v):
            self._y = v

        def print(self, value):
            print("ClassA:", value)

        def __init__(self):
            self._y = 1

    class ClassB():

        def get(self):
            print("Class A.get")
            return self._b

        def set(self, v):
            print("Class A.set:", v)
            self._b = v

        def __init__(self):
            self._b = "b"


    attributes = {
        "NewNameA": ("obj1", "y"), # works with properties
        "NewNameB": ("obj2", "get", "set") # or getter/setter methods
    }
    methods = {
        "method_renamed": ("obj1", "print"), # works with properties
    }
    pd = ParameterDirector(attributes, methods)
    pd.obj1 = ClassA()
    pd.obj2 = ClassB()

    pd.get("NewNameA") == 1
    pd.set("NewNameA", 2)
    pd.NewNameA = 33
    pd.obj1._y == 33
    pd.method_renamed(42)  # prints "42"


"""
import json
import pathlib
import collections.abc

from typing import Dict, Tuple


class ParameterDirector(object):

    _ignore_attributes = list()
    _attributes = dict()
    _methods = dict()

    def __init__(
            self,
            attributes: Dict[str, Tuple[str, str, str]] = None,
            methods: Dict[str, Tuple[str, str]] = None,
            state_file = None
    ):
        if isinstance(attributes, dict):
            self._attributes = attributes
        if isinstance(methods, dict):
            self._methods = methods
        if state_file:
            with open(state_file, 'r') as fp:
                s = fp.readlines()
                self.from_json(s)

    def resolve(self, origin):
        if len(origin) == 3:
            source, getter, setter = origin
        else:
            setter = None
            source, getter = origin
        source_obj = self
        for s in source.split("."):
            if len(s) > 0:
                source_obj = source_obj.__getattribute__(s)
        return source_obj, getter, setter

    def set(self, key, value):
        if key in self._attributes:
            origin = self._attributes[key]
            source_obj, getter_attr, setter = self.resolve(origin)
            if setter is None:
                try:
                    setattr(source_obj, getter_attr, value)
                except AttributeError:
                    pass
                    #print("Cannot set attribute", key)
            else:
                getattr(source_obj, setter)(value)
            return True
        return False

    def get(self, key):
        target = self._attributes.get(key, None)
        is_attribute = True
        if target is None:
            target = self._methods.get(key, None)
            is_attribute = False
        if target is None:
            #print("Cannt get key", key)
            raise AttributeError
        target_obj, getter, _ = self.resolve(target)
        try:
            result = getattr(target_obj, getter)
            if callable(result) and is_attribute:
                return result()
            return result
        except TypeError:
            return None

    @property
    def parameter(self):
        keys, values = list(), list()
        for key in self._attributes.keys():
            value = self.get(key)
            if value is not None:
                values.append(value)
                keys.append(key)
        return dict(zip(keys, values))

    @parameter.setter
    def parameter(self, v):
        for key, value in v.items():
            try:
                self.set(key, value)
            except AttributeError:
                pass
                #print("Cannot set key:", key)

    def __setattr__(self, key, value):
        if key in self._attributes:
            self.set(key, value)
        else:
            super().__setattr__(key, value)

    def __getattr__(self, item):
        return self.get(item)

    def to_json(self, *args, **kwargs) -> str:
        d = dict()
        for key, value in self.parameter.items():
            if key in self._ignore_attributes:
                continue
            if not isinstance(value, str) and isinstance(value, collections.abc.Iterable):
                d[key] = [v for v in value]
            else:
                d[key] = value
        s = json.dumps(d, *args, **kwargs)
        return s

    @classmethod
    def from_json(cls, s):
        if isinstance(s, pathlib.Path):
            json_str = open(s, 'r').read()
        else:
            json_str = s
        re = cls()
        re.__setstate__(json.loads(json_str))
        return re

    def save(self, filename, file_type='json', **kwargs):
        if file_type == 'json':
            s = self.to_json(**kwargs)
            with open(filename, 'w+') as fp:
                fp.write(s)
        else:
            raise NotImplementedError

    def load(self, filename, file_type='json'):
        if file_type == 'json':
            with open(filename, 'r') as fp:
                s = fp.read()
                d = json.loads(s)
                self.__setstate__(d)
        else:
            raise NotImplementedError

    def __getstate__(self):
        state = dict()
        state.update(self.parameter)
        return state

    def __setstate__(self, state):
        # see https://stackoverflow.com/questions/9310053/how-to-make-my-swig-extension-module-work-with-pickle/
        # Requires empty constructor: __init__()
        self.__init__()
        self.parameter = state
