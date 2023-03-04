from typing import runtime_checkable, Protocol


@runtime_checkable
class Dataclass(Protocol):
    __dataclass_fields__: dict
