"""Current PR #347 entry point: antipodal-only characteristic transport.

Surface fronts interact only at their antipodal focus, or traverse the mouths
through the connected throat. Worldtube histories are inputs; a compatible
field/recoil boundary problem remains unsolved. The superseded continuous
local-port action is isolated in legacy_history_throat.py for reproduction.
"""
from .antipodal_throat_transport import (
    AntipodalTransport, HistoryBoundary, MouthHistory, inertial_mouth,
)

__all__ = ['AntipodalTransport', 'HistoryBoundary', 'MouthHistory', 'inertial_mouth']
