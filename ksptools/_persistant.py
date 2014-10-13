

class PersistantObject(object):
    def __check_persistant(self):
        if not hasattr(self, '__persistant'):
            self.__persistant = list()
    
    def mapvar(self, varname, storename):
        self.__check_persistant()
        self.__persistant.append((varname, storename))
    
    def __getstate__(self):
        state = dict(self.__dict__)
        self.__check_persistant()
        state['persistant_data_'] = self.__persistant
        for v, s in self.__persistant:
            state[s] = getattr(self, v)
        return state
    
    def __setstate__(self, state):
        self.__persistant = state['persistant_data_']
        del state['persistant_data_']
        for v, s in self.__persistant:
            setattr(self, v, state[s])
            del state[s]
        self.__dict__.update(state)

