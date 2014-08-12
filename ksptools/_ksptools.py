from .simulation import runsimulation

def loadsystem(cfgfile):
    from re import match
    from ConfigParser import ConfigParser
    from .body import CelestialBody, System
    
    config_parser = ConfigParser()
    config_parser.read(cfgfile)
    system_dict = dict()
    for section in config_parser.sections():
        if match(r'cbody\.[a-zA-Z_ 0-9]+$', section) is not None:
            cbody = CelestialBody.from_config(config_parser, section, system_dict)
            system_dict[cbody.keyname] = cbody
    return System(system_dict.values())

