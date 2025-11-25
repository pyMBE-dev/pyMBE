import json
from typing import Dict
from .reaction import Reaction, ReactionParticipant


def load_reaction_set(path: str) -> Dict[str, Reaction]:
    with open(path, "r") as f:
        data = json.load(f)

    reactions = {}
    for name, rdata in data["data"].items():

        participants = [
            ReactionParticipant(**p)
            for p in rdata["participants"]
        ]

        reaction = Reaction(
            name=name,
            participants=participants,
            constant=rdata["constant"],
            reaction_type=rdata.get("reaction_type", "acid_base"),
            metadata=rdata.get("metadata")
        )

        reactions[name] = reaction

    return reactions
