from diceplayer.shared.utils.logger import Logger

from importlib import metadata

VERSION = metadata.version("diceplayer")

logger = Logger(__name__)
