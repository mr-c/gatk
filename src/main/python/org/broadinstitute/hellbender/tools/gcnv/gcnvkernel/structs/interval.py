import numpy as np
from abc import abstractmethod
from typing import Dict, List
from copy import deepcopy


class Interval:
    """ This class represents a genomic interval along with annotations
    Note:
        Equality test and hashing is based on get_key() which excludes all annotations
    """
    def __init__(self, contig: str, start: int, end: int):
        self.contig: str = str(contig)
        self.start: int = int(start)
        self.end: int = int(end)
        self.annotations = dict()
        self._hash = hash(self.get_key())

    def get_key(self):
        return self.contig, self.start, self.end

    def add_annotation(self, key: str, annotation: 'IntervalAnnotation'):
        self.annotations[key] = annotation

    def get_annotation(self, key: str):
        return self.annotations[key].get_value()

    def get_padded(self, padding: int, keep_annotations=False) -> 'Interval':
        assert padding >= 0, "padding must be >= 0"
        padded_interval = Interval(self.contig, self.start - padding, self.end + padding)
        if keep_annotations:
            padded_interval.annotations = deepcopy(self.annotations)
        return padded_interval

    def overlaps_with(self, other):
        assert isinstance(other, Interval), "the other object is not an interval!"
        if self.contig != other.contig:
            return False
        else:
            if other.start <= self.end <= other.end or other.start <= self.start <= other.end:
                return True
            else:
                return False

    def get_midpoint(self):
        return 0.5 * (self.start + self.end)

    def distance(self, other):
        if self.contig != other.contig:
            return np.inf
        else:
            return np.abs(self.get_midpoint() - other.get_midpoint())

    def __eq__(self, other):
        return self.get_key() == other.get_key()

    def __ne__(self, other):
        return self.get_key() != other.get_key()

    def __lt__(self, other):
        return self.get_key() < other.get_key()

    def __gt__(self, other):
        return self.get_key() > other.get_key()

    def __le__(self, other):
        return self.get_key() <= other.get_key()

    def __ge__(self, other):
        return self.get_key() >= other.get_key()

    def __hash__(self):
        return self._hash

    def __str__(self):
        return str(self.get_key())

    def __repr__(self):
        return self.__str__()


class IntervalAnnotation:
    def __init__(self, raw_value):
        self.parsed_value = self.parse(raw_value)

    def get_value(self):
        return self.parsed_value

    @staticmethod
    @abstractmethod
    def parse(raw_value):
        pass

    @staticmethod
    @abstractmethod
    def get_key() -> str:
        pass

    def __repr__(self):
        return repr(self.get_value())

    def __str__(self):
        return str(self.get_value())


class GCContentAnnotation(IntervalAnnotation):
    """ This class represents GC content annotation that can be added to an Interval
    """
    def __init__(self, gc_content):
        super().__init__(gc_content)

    @staticmethod
    def parse(raw_value):
        gc_content = float(raw_value)
        if not 0.0 <= gc_content <= 1.0:
            raise ValueError("GC content ({0}) must be a float in range [0, 1]".format(gc_content))
        return gc_content

    @staticmethod
    def get_key():
        return "GC_CONTENT"


class NameAnnotation(IntervalAnnotation):
    """ This class represents name annotation that can be added to an Interval
    """
    def __init__(self, name: str):
        super().__init__(name)

    @staticmethod
    def parse(raw_value):
        return str(raw_value)

    @staticmethod
    def get_key():
        return "NAME"

    def __repr__(self):
        return repr(self.get_value())


interval_annotations_dict: Dict[str, IntervalAnnotation] = {
    GCContentAnnotation.get_key(): GCContentAnnotation,
    NameAnnotation.get_key(): NameAnnotation
}

interval_annotations_dtypes: Dict[str, object] = {
    GCContentAnnotation.get_key(): float,
    NameAnnotation.get_key(): str
}


def inherit_interval_annotations(src_interval_list: List[Interval], dest_interval_list: List[Interval]):
    src_interval_list_lookup = {interval.get_key(): interval for interval in src_interval_list}
    for dest_interval in dest_interval_list:
        if dest_interval in src_interval_list_lookup:
            src_interval: Interval = src_interval_list_lookup[dest_interval]
            for key, annotation in src_interval.annotations.items():
                if key not in dest_interval.annotations.keys():
                    dest_interval.add_annotation(key, annotation)