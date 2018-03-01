#!/usr/bin/env python3

from functools import wraps

def accepts(*types):
    def decorator(f):
        @wraps(f)
        def newf(*args):

            if len(args) != len(types):
                raise AssertionError("In function '%s', expected number of inputs is %d, but %d given" % (f.__name__, len(types), len(args)))

            argtypes = tuple(map(type, args))
            if argtypes != types:
                type_str = ', '.join([str(t).split("'")[1] for t in types])
                raise TypeError("Input types must be %s" % type_str)

            return f(*args)

        return newf
    return decorator

def returns(*types):
    def decorator(f):
        @wraps(f)
        def newf(*args):

            result = f(*args)

            try:
                res_type = tuple([type(x) for x in result])
            except:
                res_type = (type(result),)

            if len(res_type) != len(types):
                raise AssertionError("In function '%s', expected number of outputs is %d, but %d returned " % (f.__name__, len(types), len(res_type)))

            if res_type != types:
                type_str = ', '.join([str(t).split("'")[1] for t in types])
                raise TypeError("Output types must be %s" % type_str)

            return result

        return newf
    return decorator

@accepts(int,int)
@returns(str)
def test_function(a, b):
    return 'a'

print(test_function(4, 5))
print(test_function.__name__)
