from os import mkdir
from os import remove
from os import listdir
from os.path import join
from os.path import exists
from json import dump
from json import load
from logging import getLogger
from contextlib import contextmanager


log = getLogger('sudan.cache')


class CacheState:
    INCOMPLETE = 1
    COMPLETE = 2


class Cache:
    METAFILE_NAME = ".cache_meta"

    def __init__(self, folder):
        self.folder = folder
        self.meta_file = join(self.folder, Cache.METAFILE_NAME)

        if not exists(self.folder):
            mkdir(self.folder)

        if not exists(self.meta_file):
            with open(self.meta_file, 'w') as ch:
                dump({'ver': 1}, ch)

        log.info(f'cache initialized {self.folder}')

        m = self.meta
        for file in listdir(self.folder):
            if not file.startswith('.') and not file in m:
                remove(join(self.folder, file))

        new_meta = {}
        for key in m:
            if exists(join(self.folder, key)):
                new_meta[key] = m[key]
        self.meta = new_meta

    @property
    def meta(self):
        with open(self.meta_file, 'r') as ch:
            return load(ch)

    @meta.setter
    def meta(self, meta_inf):
        with open(self.meta_file, 'w') as ch:
            dump(meta_inf, ch, indent=2)

    def _set_prop(self, name, state):
        m = self.meta
        m[name] = state
        self.meta = m

    @contextmanager
    def register(self, name):
        self._set_prop(name, CacheState.INCOMPLETE)
        with open(join(self.folder, name), mode='w', encoding='utf-8') as fh:
            yield fh

        # cache will be complete iff there is no exception
        self._set_prop(name, CacheState.COMPLETE)

    def cachelines(self, name, iterable):
        meta = self.meta
        if name in meta and meta[name] == CacheState.COMPLETE:
            log.info(f"cache hit {name}")
            with open(join(self.folder, name)) as f:
                yield from f
            return

        with self.register(name) as cache_handle:
            for line in iterable:
                cache_handle.write(line)
                yield line


def main():
    cache = Cache('temp')


if __name__ == '__main__':
    main()
