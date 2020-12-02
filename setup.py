import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="doubling_agent-lkmartin90",  # Replace with your own username
    version="0.0.1",
    author="Lucy Martin",
    author_email="lucy.martin@igmm.ed.ac.uk",
    description="A simple package to model the growth of a tumour from stem cells",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['numpy', 'pandas', 'scipy', 'matplotlib'],
    scripts=['bin/tumour_sim.py', 'bin/animate_data.py', 'bin/plot_data.py', 'bin/steve_sim.py',
             'bin/plot_steve.py', 'bin/read_experimental_data.py', 'bin/state_sim.py', 'bin/density_sim.py'],
)
