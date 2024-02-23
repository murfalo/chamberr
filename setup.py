from setuptools import setup, find_packages

setup(
    name="chamber",
    version="0.1.0",
    author="Murphy Angelo",
    author_email="howdy@murfalo.com",
    description="Converts AMBER forcefields to CHARMM-readable formats",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/murfalo/chamber",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[
        "parmed",
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
            "chamber=chamber.__main__:main",
        ],
    },
)

