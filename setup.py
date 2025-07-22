from setuptools import setup, find_packages

setup(
    name='torchic',
    version='0.5.25',  # Initial version
    author='Giorgio Alberto Lucia, Roberta Ferioli',
    author_email='giogioalberto@gmail.com',
    description='A Python library for data analysis in python and ROOT (tailored for high energy physics)',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',  # Format of README
    url='https://github.com/GiorgioAlbertoLucia/torchic',
    packages=find_packages(),  # Automatically finds your package
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',  # Replace with your license
        'Operating System :: OS Independent',
    ],
    package_data={'torchic': [
        'physics/BetheBloch.hh',
        'roopdf/RooCustomPdfs/RooGausExp.hh',
        'roopdf/RooCustomPdfs/RooGausExp.cxx',
        'roopdf/RooCustomPdfs/RooSillPdf.hh',
        'physics/simulations/ExponentialDecaySimulation.cxx',
        'physics/simulations/TwoBodyDecaySimulation.cxx',
        'physics/simulations/test.cxx',
        ]},  # Include header files
    python_requires='>=3.6',  # Minimum Python version
    install_requires=[
        # List your package dependencies here
        'numpy<2.0',
        'pandas',
        'scikit-learn',
        'pyYAML',
        'uproot',
    ],
    extras_require={
        'ROOT': ['ROOT'],
    },
)
