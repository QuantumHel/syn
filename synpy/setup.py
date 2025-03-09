from pathlib import Path
from setuptools import setup  # type: ignore

synpy_core: str = (Path(__file__).parent / "synpy_core").as_uri()

setup(
    install_requires=[
        f"synpy_core @ {synpy_core}",
    ]
)
