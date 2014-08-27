from itertools import count


class TreeNode(object):
    def __init__(self, parent=None):
        self._ident = next(TreeNode._nextid)
        self._parent = parent
        self._children = set()
    
    def addnode(self, node):
        self._children.add(node)
        node._parent = self
    
    def removenode(self, node):
        self._children.remove(node)
        node._parent = None
    
    def __eq__(self, other):
        return self._ident == other._ident
    
    def __hash__(self):
        return hash(self._ident)
    
    _nextid = count()


class ModelNode(TreeNode):
    def __init__(self, models, parent=None):
        TreeNode.__init__(self, parent):
        self._modelset = set(models)
    
    def canapply_to_model(self, model):
        return 'any' in self._modelset or model in self._modelset


