from diceplayer.shared.config.dice_dto import DiceDTO
from diceplayer.shared.external import External

from multiprocessing import Process, connection
from setproctitle import setproctitle
import sys


class Dice(External):
    __slots__ = ['config', 'step']

    def __init__(self, data: dict):
        self.config: DiceDTO = self.set_config(data)

    @staticmethod
    def set_config(data: dict) -> DiceDTO:
        return DiceDTO.from_dict(data)

    def configure(self, step: any):
        self.step = step

    def start(self, cycle: int):
        procs = [
            Process(target=self._simulation_process, args=(cycle, proc))
            for proc in range(1, self.config.ncores+1)
        ]

        for proc in procs:
            proc.start()

        connection.wait(p.sentinel for p in procs)

    def reset(self):
        del self.step

    def _simulation_process(self, cycle: int, proc: int):
        setproctitle(f"diceplayer-step{cycle:0d}-p{proc:0d}")

        try:
            self._make_proc_dir(cycle, proc)
            self._make_dice_inputs(cycle, proc)
            self._run_dice(cycle, proc)
        except Exception as err:
            sys.exit(err)

    def _make_proc_dir(self, cycle, proc):
        raise NotImplementedError

    def _make_dice_inputs(self, cycle, proc):
        raise NotImplementedError

    def _run_dice(self, cycle, proc):
        raise NotImplementedError
