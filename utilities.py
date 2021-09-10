

class IOWrapper:
    """ wraps iterator with io buffer - like methods
    """

    def __init__(self, obj):
        self._iter = obj

    def read(self):
        pass

    def readline(self):
        try:
            return next(self._iter)
        except StopIteration:
            return ""

    def close(self):
        pass


class ToolExecutionError(Exception):

    def __init__(self, name, error):
        self.name = name
        self.error = error
        super().__init__(f'{self.name} error: {self.error}')
