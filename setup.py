import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

exec(open("mtneedlet/version.py").read())

setuptools.setup(
    name="mtneedlet",
    version=__version__,
    author="Javier Carron Duque",
    author_email="javier.carron@roma2.infn.it",
    description="A Python package for needlet filtering and Multiple Testing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    install_requires=["numpy","healpy","matplotlib","scipy","pandas"],
    classifiers=["Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
