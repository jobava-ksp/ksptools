import ConfigParser
import re

from . import body

def loadsystem(cfgfile):
    config_parser = ConfigParser.ConfigParser()
    config_parser.read(cfgfile)
    system_dict = dict()
    for section in config_parser.sections():
        if re.match(r'cbody\.[a-zA-Z_ 0-9]+$', section) is not None:
            key, cbody = body.CelestialBody.from_config(config_parser, section, system_dict)
            system_dict[key] = cbody
    return system_dict

