from setuptools import setup, find_packages

setup(
    name="mcpy",
    version="1.0",
    python_requires=">=3.12",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'mcpy': ['*.so'],
    },
    install_requires=[
        "numpy==2.1.2",
        "scipy==1.14.1",
        "matplotlib==3.9.2",
        "pybind11==2.13.6",
    ],
)