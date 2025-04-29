from setuptools import find_packages
from setuptools import setup

setup(
    name="sveqtl",
    version="1.0.dev",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "sveqtl": ["py.typed"],
    },
)
