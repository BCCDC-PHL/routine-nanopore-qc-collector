from setuptools import setup, find_namespace_packages


setup(
    name='routine-nanopore-qc-collector',
    version='0.1.0',
    packages=find_namespace_packages(),
    entry_points={
        "console_scripts": [
            "routine-nanopore-qc-collector = routine_nanopore_qc_collector.__main__:main",
        ]
    },
    scripts=[],
    package_data={
    },
    install_requires=[
    ],
    description='Collect Routine Nanopore Sequence QC Data',
    url='https://github.com/BCCDC-PHL/routine-nanopore-qc-collector',
    author='Dan Fornika',
    author_email='dan.fornika@bccdc.ca',
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
