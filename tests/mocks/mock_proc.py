from typing_extensions import List

import itertools


class MockProc:
    pid_counter = itertools.count()

    def __init__(self, *args, **kwargs):
        self.pid = next(MockProc.pid_counter)

        if "exitcode" in kwargs:
            self.exitcode = kwargs["exitcode"]
        else:
            self.exitcode = 0

        self.sentinel = self.pid

    def __call__(self, *args, **kwargs):
        return self

    def start(self):
        pass

    def terminate(self):
        pass


class MockConnection:
    @staticmethod
    def wait(sentinels: List[int]):
        return sentinels
