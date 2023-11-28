from setuptools import setup, find_packages
import os

# readme_path = os.path.join(os.path.dirname(__file__), "README.md")
# with open(readme_path, "r") as readme_file:
#     readme = readme_file.read()

# with open('requirements.txt', 'r') as f:
#     required = f.read().splitlines()

# with open("morph_utils/__init__.py") as f:
#     for line in f:
#         if line.startswith("__version__"):
#             version = line.split('"')[1]

setup(
    name = 'primate_uprighting',
    version = '0.0.1', #version,
    description = "Functions for common and not so common morphology operations",
    # long_description=readme,
    author = "Sarah Walling-Bell",
    author_email = "sarah.wallingbell@alleninstitute.org",
    # url = 'https://github.com/MatthewMallory/morph_utils',
    packages = find_packages(),
    # install_requires = required,
    include_package_data=True,
    package_data={"morph_utils": ["tests/*"]},
    setup_requires=['pytest-runner'],
)