import logging
from pathlib import Path


def valid_logger(func):
    def wrapper(*args, **kwargs):
        logger = args[0]
        assert logger._was_set, "Logger is not set. Please call set_logger() first."

        return func(*args, **kwargs)

    return wrapper


class Logger:
    outfile = None

    _logger = None

    _was_set = False

    def __init__(self, logger_name):
        if self._logger is None:
            self._logger = logging.getLogger(logger_name)

    def set_logger(self, outfile="run.log", level=logging.INFO, stream=None):
        outfile_path = None
        if outfile is not None and stream is None:
            outfile_path = Path(outfile)
            if outfile_path.exists():
                outfile_path.rename(str(outfile_path) + ".backup")

        if level is not None:
            self._logger.setLevel(level)

        self._create_handlers(outfile_path, stream)

        self._was_set = True

    @valid_logger
    def info(self, message):
        self._logger.info(message)

    @valid_logger
    def debug(self, message):
        self._logger.debug(message)

    @valid_logger
    def warning(self, message):
        self._logger.warning(message)

    @valid_logger
    def error(self, message):
        self._logger.error(message)

    def _create_handlers(self, outfile_path: Path, stream):
        handlers = []
        if outfile_path is not None:
            handlers.append(logging.FileHandler(outfile_path, mode="a+"))
        elif stream is not None:
            handlers.append(logging.StreamHandler(stream))
        else:
            handlers.append(logging.StreamHandler())

        for handler in handlers:
            handler.setFormatter(logging.Formatter("%(message)s"))
            self._logger.addHandler(handler)

    def close(self):
        for handler in self._logger.handlers:
            handler.close()
            self._logger.removeHandler(handler)
