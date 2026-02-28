from functools import cached_property


def invalidate_computed_properties():
    """
    Decorator function to invalidate the cached properties of the molecule when a new atom is added

    Args:
        properties (list[str]): list of the names of the properties to be invalidated
    """

    def get_cached_properies(cls: type) -> set[str]:
        return {
            name
            for name, value in cls.__dict__.items()
            if isinstance(value, cached_property)
        }

    def decorator(func):
        def wrapper(self, *args, **kwargs):
            result = func(self, *args, **kwargs)
            for prop in get_cached_properies(self.__class__):
                if hasattr(self, prop):
                    delattr(self, prop)
            return result

        return wrapper

    return decorator
