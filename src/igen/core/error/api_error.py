from typing import Optional


class ApiError(Exception):
    def __init__(self, message: str, status_code: int = 500, from_exception: Optional[Exception] = None):
        super().__init__(message)

        self.message = message
        self.status_code = status_code

        if from_exception is not None:
            self.__cause__ = from_exception

    def __str__(self) -> str:
        full_message = f"{self.message}\nStatus Code: {self.status_code}"
        cause = getattr(self, "__cause__", None)

        while cause:
            full_message += f"\nCaused by: {cause}"
            cause = getattr(cause, "__cause__", None)

        return full_message
