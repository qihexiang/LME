from copy import deepcopy
from libs.constants import PRODUCTION


class StateContainer:
    def __init__(self, init_state=None) -> None:
        self.__state__ = init_state
        self.subscribers = set()

    @property
    def state(self):
        # In prodcution mode, return the __state__ directly
        if PRODUCTION:
            return self.__state__
        # In development mode, return an copied state to avoid change the state from outside.
        return deepcopy(self.__state__)

    def add_subscriber(self, subscriber):
        self.subscribers.add(subscriber)

    def remove_subscriber(self, subscriber):
        self.subscribers.remove(subscriber)

    def __broadcast__(self):
        for subscriber in self.subscribers:
            subscriber(self.state)

    def update(self, updator):
        self.__state__ = updator(self.state)
        self.__broadcast__()
