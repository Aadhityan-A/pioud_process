from setuptools import setup, find_packages

setup(
    name="pioud_process",
    version="0.1.0",
    packages=find_packages() + ['gui'],  # Explicitly include the gui package
    install_requires=[
        "numpy",
        "ase",
        "scipy",
        "dscribe",
        "skmatter",
    ],
    author="Aadhityan A",
    author_email="aadhityan1995@gmail.com",
    description="A package to process ouput of PIOUD",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/pioud_process",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    entry_points={
        'console_scripts': [
            'pioud-process=gui.run_gui:main',
        ],
    },
    include_package_data=True,
    package_data={
        'gui': ['*.py'],
    },
)