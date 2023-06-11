from typing import Protocol, runtime_checkable


@runtime_checkable
class Dataclass(Protocol):
    __dataclass_fields__: dict
