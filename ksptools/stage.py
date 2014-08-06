import collections

from .part import ResourceSink

class PhysicalStage(object):
    def __init__(self, parts, prev=None):
        self.prev = prev
        self.next = None
        if prev is not None:
            prev.next = self
        self.parts = parts
        partsbyclass = collections.defaultdict(list)
        for part in parts:
            for cls in part.partclass:
                partsbyclass[cls].append(part)
        self.partsbyclass = dict(partsbyclass)
        self.resources = ResourceSink.combine(parts)
        self.feeds_in = []
        self.feeds_out = []
