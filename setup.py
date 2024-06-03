from setuptools import find_packages, setup

setup(
    name="chamberr",
    version="0.1.0",
    author="Murphy Angelo",
    author_email="howdy@murfalo.com",
    description="Converts AMBER forcefields to CHARMM-readable formats",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/murfalo/chamberr",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[
        # Requires custom for of ParmEd pending upstream merge (see
        # `README.md`)
        # "parmed",
        # "ambertools",  # must install with conda :(
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    license="MIT",
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "chamberr=chamberr.__main__:main",
        ],
    },
)
