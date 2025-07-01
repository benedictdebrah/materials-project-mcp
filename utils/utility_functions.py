from pydantic import Field
from dataclasses import dataclass
from pathlib import Path
from datetime import datetime
from typing import Dict, List
import os



MARKDOWN_FOLDER = Path(os.path.dirname(__file__)).parent / "markdown_folder"
print(MARKDOWN_FOLDER)

@dataclass
class MarkdownFile: 
    name: str 
    path: Path 
    content: str = Field(default_factory=str)


class MarkdownResourceManager:
    def __init__(self, folder_path: Path):
        self.folder_path = folder_path
        self.files: Dict[str, MarkdownFile] = {}


    def load_files(self): 
        """"Load all markdown files from the folder path"""

        self.files.clear()

        if not self.folder_path.exists():
            print(f"Folder {self.folder_path} does not exist.")
            return "Folder does not exist."
        
        for file_path in self.folder_path.glob("*.md"):
            try: 
                content = file_path.read_text(encoding="utf-8")
                last_modified = datetime.fromtimestamp(file_path.stat().st_mtime)
                markdown_file = MarkdownFile(
                    name=file_path.stem,
                    path=file_path,
                    content=content,
                 #   last_modified=last_modified
                )
                self.files[file_path.stem] = markdown_file
            except Exception as e:
                print(f"Error reading file {file_path}: {e}")


if __name__ == "__main__":
    manager = MarkdownResourceManager(MARKDOWN_FOLDER)
    manager.load_files()
    print(manager.files["docsmaterialsproject"].content)
