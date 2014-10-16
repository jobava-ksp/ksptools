import itertools

from .._persistant import PersistantObject

class Node(PersistantObject):
    def __init__(self, parent_node=None):
        self._id = next(Node._nextid)
        self.children = set()
        self.parent_node = None
        if parent_node is not None:
            parent_node._add(self)
        self.mapvar('_id', 'node_id_')
    
    def __eq__(self, other):
        return self._id == other._id
    
    def _add(self, new_child):
        if new_child.parent_node is not None:
            new_child.parent_node._remove(new_child)
        new_child.parent_node = self
        self.children.add(new_child)
    
    def _remove(self, old_child):
        self.children.remove(old_child)
        old_child.parent_node = None
    
    @staticmethod
    def _path_to_root(from_node):
        root_node = from_node
        while root_node is not None:
            yield root_node
            root_node = root_node.parent_node
    
    @staticmethod
    def _path(from_node, to_node):
        root_path_from = list(Node._path_to_root(from_node))
        root_path_to = list(Node._path_to_root(to_node))
        while len(root_path_from) > 0 and len(root_path_to) > 0:
            if root_path_from[-1] == root_path_to[-1]:
                root = root_path_from[-1]
                root_path_from = root_path_from[:-1]
                root_path_to = root_path_to[:-1]
            else:
                break
        return root_path_from + [root], reversed(root_path_to)
    
    def walk(self, to_node, x, fdown, fup):
        iterdown, iterup = Node._path(self, to_node)
        prev = self
        for next in list(iterdown)[1:]:
            x = fdown(prev, next, x)
            prev = next
        for next in list(iterup):
            x = fup(prev, next, x)
            prev = next
        return x
    
    def __setstate__(self, state):
        PersistantObject.__setstate__(self, state)
        Node._nextid = itertools.count(max(next(Node._nextid), self._id+1))
    
    _nextid = itertools.count()


