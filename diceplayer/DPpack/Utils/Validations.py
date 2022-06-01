def NotNull(requiredArgs=[]):
    def _NotNull(function):
        def wrapper(*args, **kwargs):
            for arg in requiredArgs:
                try:
                    assert (
                        kwargs.get(arg) is not None
                    ), "Invalid Config File. Keyword {} is required".format(arg)
                except AssertionError as err:
                    print(err)
            return function(*args, **kwargs)

        return wrapper

    return _NotNull