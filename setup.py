from setuptools import setup, find_packages

setup(
    name="muffin",
    version="0.1.0",
    packages=find_packages(),
    package_dir={"": "."},
    package_data={
        "PythonTools": ["*.py"],
    },
    install_requires=[
        "numpy",
        "matplotlib",
        "scipy",
    ],
    python_requires=">=3.6",
    description="MUFFIN - Multi-fluid plasma simulation code",
    author="Laboratory of Plasma Physics",
    author_email="",
    url="https://github.com/LaboratoryOfPlasmaPhysics/MUFFIN",
)