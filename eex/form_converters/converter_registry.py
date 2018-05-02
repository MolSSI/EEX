
def sanitize_term_order_name(order):
    if isinstance(order, str):
        order = order.lower()

    if order in [2, "two", "bond", "bonds"]:
        return 2
    elif order in [3, "three", "angle", "angles"]:
        return 3
    elif order in [4, "four", "dihedral", "dihedrals"]:
        return 4
    else:
        raise KeyError("EEX: Term order name '%s' not recognized." % str(order))


REGISTERED_CONVERTERS = {}

for order in (2, 3, 4):
    REGISTERED_CONVERTERS[order] = {}


def register_converter(order=None):
    order = sanitize_term_order_name(order)

    def decorator_function(fn):
        REGISTERED_CONVERTERS[order][fn.__name__] = fn

        def wrapper_function(*args, **kwargs):
            return fn(*args, **kwargs)
        return wrapper_function
    return decorator_function
