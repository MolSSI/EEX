from .. import metadata

REGISTERED_CONVERTERS = {}

for order in (2, 3, 4):
    REGISTERED_CONVERTERS[order] = {}


def register_converter(order=None):
    order = metadata.md_helper.sanitize_term_order_name(order)

    def decorator_function(fn):
        REGISTERED_CONVERTERS[order][fn.__name__.lower()] = fn

        def wrapper_function(*args, **kwargs):
            return fn(*args, **kwargs)

        return wrapper_function

    return decorator_function
