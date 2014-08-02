import functools

_debug_switches = dict()

def define(switch_names, default=True):
    for s in switch_names:
        if s not in _debug_switches:
            _debug_switches[s] = default

def undefine(switch_names):
    for s in switch_names:
        if s in _debug_switches:
            del _debug_switches[s]

def enable(switch_names, value=True):
    for s in switch_names:
        _debug_switches[s] = value

def trace(switch_names):
    define(switch_names, False)
    def trace_decerator(func):
        @functools.wraps(func)
        def wrapper_func(*args, **kwargs):
            if any(_debug_switches[s] for s in switch_names):
                return func(*args, **kwargs)
            else:
                return None
        return wrapper_func
    return trace_decerator

