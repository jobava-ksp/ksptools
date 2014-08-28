from itertools import count


class TreeNode(object):
    def __init__(self, parent=None):
        """
        :type parent: TreeNode
        """
        self._ident = next(TreeNode._nextid)
        self._parent = parent
        self._children = set()
    
    def addnode(self, node):
        """
        :type node: TreeNode
        """
        self._children.add(node)
        node._parent = self
    
    def removenode(self, node):
        """
        :type node: TreeNode
        """
        self._children.remove(node)
        node._parent = None
    
    def __eq__(self, other):
        """
        :type other: TreeNode
        """
        return self._ident == other._ident
    
    def __hash__(self):
        return hash(self._ident)
    
    _nextid = count()


class ModelNode(TreeNode):
    def __init__(self, models, parent=None):
        """
        :type models: list
        :type parent: TreeNode
        """
        TreeNode.__init__(self, parent)
        self._modelset = set(models)
    
    def inmodel(self, model):
        """
        :type model: str
        """
        return 'any' in self._modelset or model in self._modelset


