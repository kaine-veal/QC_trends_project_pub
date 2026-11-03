from setuptools import setup, find_packages

setup(
    name="qm_trends_project",
     version="0.1.0",
    package_dir={"": "package"},
    packages=find_packages(where="package"),
    python_requires=">=3.8",
)