import itertools


class Node(object):
    def __init__(self, parent_node=None):
        self._id = next(Node._nextid)
        self.children = set()
        self.parent_node = None
        if parent_node is not None:
            parent_node._add(self)
    
    def __eq__(self, other):
        self._id == other._id
    
    def __hash__(self):
        return hash(self._id)
    
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
            root_path_from = root_path_from + [root_node]
            root_node = root_node.parent_node
    
    @staticmethod
    def path(from_node, to_node):
        root_path_from = list(Node._path_to_root(from_node))
        root_path_to = list(reversed(Node._path_to_root(to_node)
        while root_path_from[-1] == root_path_to[0]:
            root_path_from = root_path_from[:-1]
            root_path_to = root_path_to[1:]
        path_ft = root_path_from + foot_path_to
        return zip(path_ft[:-1], path_ft[1:])
    
    nextid = itertools.count()


