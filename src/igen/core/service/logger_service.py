from __future__ import annotations

import logging
from logging import FileHandler, Formatter, Logger, StreamHandler
from threading import RLock


class LoggerService:
    """Thread-safe singleton loggers configured with stream and file handlers."""

    _instances: dict[str, "LoggerService"] = {}
    _lock = RLock()

    def __new__(
        cls, name: str = "app", level: int = logging.INFO, log_file: str | None = None, root_path: str | None = None
    ):
        with cls._lock:
            if name not in cls._instances:
                instance = super().__new__(cls)
                instance._initialize(name, level, log_file, root_path)
                cls._instances[name] = instance

            return cls._instances[name]

    def _initialize(self, name: str, level: int, log_file: str | None, root_path: str | None = None):
        self.logger: Logger = logging.getLogger(name)
        self.logger.setLevel(level)
        formatter = Formatter("[%(asctime)s] [%(levelname)s] %(name)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

        stream_handler = StreamHandler()
        stream_handler.setFormatter(formatter)
        self.logger.addHandler(stream_handler)

        if root_path and log_file:
            logs_dir = root_path / "logs"
            logs_dir.mkdir(parents=True, exist_ok=True)
            log_path = logs_dir / log_file

            file_handler = FileHandler(filename=log_path, mode="a")
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self) -> Logger:
        return self.logger
