import logging


class Logger:
    outfile = None

    _logger = None
    _instance = None

    def __new__(cls, *args, **kwargs):
        if not getattr(cls, '_instance'):
            cls._instance = super(Logger, cls).__new__(cls)
        return cls._instance

    def set_logger(self, logger_name, outfile='run.log', level=logging.INFO):
        self.outfile = outfile

        self._logger = logging.getLogger(logger_name)

        if level is not None:
            self._logger.setLevel(level)

        self._create_handlers()

    def _create_handlers(self):
        handlers = []
        if self.outfile is not None:
            handlers.append(logging.FileHandler(self.outfile, mode='a+'))
        else:
            handlers.append(logging.StreamHandler())

        for handler in handlers:
            handler.setFormatter(logging.Formatter('%(message)s'))
            self._logger.addHandler(handler)

    def close(self):
        for handler in self._logger.handlers:
            handler.close()
            self._logger.removeHandler(handler)
