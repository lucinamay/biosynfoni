import os


class ChangeDirectory:
    def __init__(self, new_path):
        self.new_path = new_path
        self.original_path = None

    def __enter__(self):
        self.original_path = os.getcwd()
        os.makedirs(self.new_path, exist_ok=True)
        os.chdir(self.new_path)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.original_path is not None:
            os.chdir(self.original_path)
