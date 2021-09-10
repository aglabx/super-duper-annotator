from os.path import join
from os.path import exists
from os.path import abspath
from settings import settings


class Tool:

    def __init__(self, input):
        if not exists(input):
            raise FileExistsError(input)

    def execute(self):
        pass

    @property
    def output(self):
        pass
