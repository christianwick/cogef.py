import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cogef.py",
    author="Christian Wick",
    author_email="christian.wick@fau.de",
    description="Cogef.py - COnstraint Geometries simulate external force",
    long_description=long_description,
    long_description_content_type="text/markdown",

    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Linux",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    use_scm_version = {
        "write_to": "src/cogef/_version.py",
        "write_to_template": "__version__ = '{version}'"},
    setup_requires = ["setuptools_scm"]
)
