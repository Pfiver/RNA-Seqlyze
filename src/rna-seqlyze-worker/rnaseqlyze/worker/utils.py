"""
Some helpful infrastructure for the WorkerStages class
"""

class OrderedClass(type):
    """
    Adds a ".members" attribute to classes that use this
    class as a metaclass. The members attribute contains
    an ordered(!) list of members of that class.

    Used as a metaclass of rnaseqlyze.worker.stages.WorkerStages --
    all methods defined in that class will be automatically called during
    analysis processing.

    Taken from
    http://docs.python.org/dev/py3k/reference/datamodel.html#metaclass-example
    """

    @classmethod
    def __prepare__(metacls, name, bases, **kwds):
        return collections.OrderedDict()

    def __new__(cls, name, bases, namespace, **kwds):
        result = type.__new__(cls, name, bases, dict(namespace))
        result.members = filter(lambda s: not s.startswith('_'), namespace)
        return result

