from pydash import py_

class UUIDPair:
    def __init__(self, ids) -> None:
        self.a, self.b = ids

    def has_uuid(self, target_id) -> bool:
        return self.another_atom(target_id) is not None

    def another_atom(self, atom_id):
        if self.a == atom_id:
            return self.b

        if self.b == atom_id:
            return self.a

        return None

    def __eq__(self, another) -> bool:
        return (another.a == self.a and another.b == self.b) or (
            another.a == self.b and another.b == self.a
        )

    def __hash__(self) -> int:
        sorted_id = py_.sort([self.a, self.b])
        return hash(tuple(sorted_id))

    def __repr__(self) -> str:
        return f"{self.a} {self.b}"

    @property
    def export(self):
        return f"{str(self.a)} {str(self.b)}"

